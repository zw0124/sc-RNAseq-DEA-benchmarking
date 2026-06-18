#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library(argparser)
	library(Seurat)
	library(Matrix)
})

parser <- arg_parser("Pseudobulk t-test")
parser <- add_argument(parser, "--input", help = "Input Seurat RDS", type = "character")
parser <- add_argument(parser, "--output", help = "Output TSV", type = "character")
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

sample_factor <- factor(sample_vec)
agg_mat <- Matrix::sparse.model.matrix(~0 + sample_factor)
colnames(agg_mat) <- levels(sample_factor)

cond_levels <- unique(cond_vec)
if (length(cond_levels) != 2) {
	stop(sprintf("ttest requires exactly 2 conditions, found: %d", length(cond_levels)))
}

sample_meta <- data.frame(Sample = levels(sample_factor), stringsAsFactors = FALSE)
sample_to_condition <- tapply(cond_vec, sample_vec, function(x) unique(x)[1])
sample_meta$Condition <- as.character(sample_to_condition[sample_meta$Sample])
sample_meta$Condition <- factor(sample_meta$Condition)
rownames(sample_meta) <- sample_meta$Sample

group1 <- rownames(sample_meta)[sample_meta$Condition == levels(sample_meta$Condition)[1]]
group2 <- rownames(sample_meta)[sample_meta$Condition == levels(sample_meta$Condition)[2]]

cell_library_size <- Matrix::colSums(counts)
cell_library_size[cell_library_size == 0] <- 1
cell_scale <- Matrix::Diagonal(x = 1e6 / as.numeric(cell_library_size))
cell_cpm <- counts %*% cell_scale
sample_sum_cpm <- cell_cpm %*% agg_mat
sample_ncells <- Matrix::colSums(agg_mat)
sample_ncells[sample_ncells == 0] <- 1
sample_mean_cpm <- sweep(sample_sum_cpm, 2, sample_ncells, "/")
sample_mean_cpm <- as.matrix(sample_mean_cpm)

mean_cpm_group1 <- rowMeans(sample_mean_cpm[, group1, drop = FALSE])
mean_cpm_group2 <- rowMeans(sample_mean_cpm[, group2, drop = FALSE])
logfc_cpm <- log2((mean_cpm_group2 + 1e-9) / (mean_cpm_group1 + 1e-9))

safe_ttest <- function(x1, x2) {
	if (length(x1) < 2 || length(x2) < 2) {
		return(NA_real_)
	}
	if (var(x1) == 0 && var(x2) == 0) {
		return(1.0)
	}
	tryCatch(t.test(x2, x1)$p.value, error = function(e) NA_real_)
}

p_values <- vapply(
	seq_len(nrow(sample_mean_cpm)),
	function(i) safe_ttest(sample_mean_cpm[i, group1], sample_mean_cpm[i, group2]),
	numeric(1)
)

out <- data.frame(
	gene = rownames(sample_mean_cpm),
	logFC = as.numeric(logfc_cpm),
	p_val = as.numeric(p_values),
	stringsAsFactors = FALSE
)
write.table(out, file = args$output, sep = "\t", row.names = FALSE, quote = FALSE)
