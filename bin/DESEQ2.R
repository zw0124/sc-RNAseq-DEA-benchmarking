#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library(argparser)
	library(Seurat)
	library(Matrix)
	library(DESeq2)
})

parser <- arg_parser("Pseudobulk DESeq2")
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
reduced_formula <- if (use_batch) ~ Batch else ~ 1

dds <- DESeqDataSetFromMatrix(countData = pb_counts, colData = sample_meta, design = design_formula)
if (args$lfc_threshold > 0) {
	dds <- DESeq(dds, test = "Wald", fitType = "local")
} else {
	dds <- DESeq(dds, test = "LRT", reduced = reduced_formula, fitType = "local")
}
cond_name <- grep("^Condition", resultsNames(dds), value = TRUE)
if (length(cond_name) == 0) {
	stop("No Condition coefficient found in DESeq2 results names")
}
cond_name <- cond_name[1]

if (args$lfc_threshold > 0) {
	res <- as.data.frame(results(dds, name = cond_name, lfcThreshold = args$lfc_threshold, altHypothesis = "greaterAbs"))
} else {
	res <- as.data.frame(results(dds, name = cond_name))
}

out <- data.frame(
	gene = rownames(res),
	logFC = res$log2FoldChange,
	p_val = res$pvalue,
	stringsAsFactors = FALSE
)
write.table(out, file = args$output, sep = "\t", row.names = FALSE, quote = FALSE)
