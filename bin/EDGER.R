#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library(argparser)
	library(Seurat)
	library(Matrix)
	library(edgeR)
})

parser <- arg_parser("Pseudobulk edgeR")
parser <- add_argument(parser, "--input", help = "Input Seurat RDS", type = "character")
parser <- add_argument(parser, "--output", help = "Output TSV", type = "character")
parser <- add_argument(parser, "--lfc_threshold", help = "LFC threshold", type = "numeric", default = 0.0)
args <- parse_args(parser)

obj <- readRDS(args$input)
if (!("RNA" %in% names(obj@assays))) {
	stop("Seurat object is missing RNA assay")
}

counts <- Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")
cell_meta <- obj@meta.data
sample_col <- "Sample"
condition_col <- "Condition"
if (!(sample_col %in% colnames(cell_meta)) || !(condition_col %in% colnames(cell_meta))) {
	stop("Missing required metadata columns: Sample and Condition")
}

sample_vec <- as.character(cell_meta[[sample_col]])
cond_vec <- as.character(cell_meta[[condition_col]])
batch_exists <- "Batch" %in% colnames(cell_meta)
batch_vec <- if (batch_exists) as.character(cell_meta[["Batch"]]) else rep("Batch1", nrow(cell_meta))

sample_factor <- factor(sample_vec)
agg_mat <- Matrix::sparse.model.matrix(~0 + sample_factor)
colnames(agg_mat) <- levels(sample_factor)
pb_counts <- counts %*% agg_mat
pb_counts <- as.matrix(pb_counts)
mode(pb_counts) <- "integer"

sample_ids <- colnames(pb_counts)
sample_meta <- data.frame(Sample = sample_ids, stringsAsFactors = FALSE)
sample_to_condition <- tapply(cond_vec, sample_vec, function(x) unique(x)[1])
sample_meta$Condition <- as.character(sample_to_condition[sample_ids])
sample_to_batch <- tapply(batch_vec, sample_vec, function(x) unique(x)[1])
sample_meta$Batch <- as.character(sample_to_batch[sample_ids])
rownames(sample_meta) <- sample_meta$Sample
sample_meta$Condition <- factor(sample_meta$Condition)
sample_meta$Batch <- factor(sample_meta$Batch)

batch_n_unique <- if (batch_exists) length(unique(as.character(sample_meta$Batch))) else 0
use_batch <- batch_exists && (batch_n_unique > 1)
design_formula <- if (use_batch) ~ Batch + Condition else ~ Condition

dge <- DGEList(counts = pb_counts)
dge <- calcNormFactors(dge)

sample_meta_edger <- sample_meta[colnames(dge), , drop = FALSE]
design_edger <- model.matrix(design_formula, data = sample_meta_edger)
cond_coef_edger <- grep("^Condition", colnames(design_edger))

dge <- estimateDisp(dge, design_edger)
fit <- glmFit(dge, design_edger)
if (args$lfc_threshold > 0) {
	test_edger <- glmTreat(fit, coef = cond_coef_edger[1], lfc = args$lfc_threshold)
	res <- topTags(test_edger, n = Inf, sort.by = "none")$table
} else {
	lrt <- glmLRT(fit, coef = cond_coef_edger[1])
	res <- topTags(lrt, n = Inf, sort.by = "none")$table
}

out <- data.frame(
	gene = rownames(res),
	logFC = res$logFC,
	p_val = res$PValue,
	stringsAsFactors = FALSE
)
write.table(out, file = args$output, sep = "\t", row.names = FALSE, quote = FALSE)
