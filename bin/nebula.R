#!/usr/bin/env Rscript

library(nebula)
library(argparser)

p <- arg_parser("Run nebula on a Seurat object")
p <- add_argument(p, "--input", help="Path to input Seurat object (RDS file)", type="character")
p <- add_argument(p, "--scenario", help="Scenario name", type="character", default="")
p <- add_argument(p, "--output", help="Path to output results file (TSV)", type="character")
p <- add_argument(p, "--n_cores", help="Number of cores to use", type="integer", default=4)

argv <- parse_args(p)
obj <- readRDS(argv$input)
output_file <- argv$output

required_cols <- c("Sample", "Condition", "nCount_RNA")
missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required metadata columns: %s", paste(missing_cols, collapse=", ")))
}

has_batch_col <- "Batch" %in% colnames(obj@meta.data)
pred_cols <- if (has_batch_col) c("Condition", "Batch") else "Condition"

seuratdata <- scToNeb(obj = obj,
                      assay = "RNA",
                      id = "Sample",
                      pred = pred_cols,
                      offset = "nCount_RNA")

pred_df <- seuratdata$pred
condition_levels <- unique(as.character(pred_df$Condition))
if (length(condition_levels) != 2) {
  stop(sprintf("nebula requires exactly 2 Condition levels, found: %d", length(condition_levels)))
}
pred_df$Condition <- as.integer(as.character(pred_df$Condition) == condition_levels[2])

batch_n_unique <- if ("Batch" %in% colnames(pred_df)) length(unique(as.character(pred_df$Batch))) else 0
use_batch <- (batch_n_unique > 1)

if (use_batch) {
  cat(sprintf("Using design: ~ Condition + Batch (n_unique_batch=%d)\n", batch_n_unique))
  df <- model.matrix(~Condition + Batch, data = pred_df)
} else {
  cat("Using design: ~ Condition (Batch missing or n_unique_batch <= 1)\n")
  df <- model.matrix(~Condition, data = pred_df)
}

data_g <- group_cell(count = seuratdata$count,
                     id = seuratdata$id,
                     pred = df,
                     offset = seuratdata$offset)

if (is.null(data_g)) {
  cat("Cells are already grouped. Using original data.\n")
  re <- nebula(seuratdata$count, seuratdata$id,
               pred = df, offset = seuratdata$offset, ncore=argv$n_cores)
} else {
  cat("Cells have been regrouped. Using grouped data.\n")
  re <- nebula(data_g$count, data_g$id,
               pred = data_g$pred, offset = data_g$offset, ncore=argv$n_cores)
}

logfc_col <- "logFC_Condition"
p_col <- "p_Condition"
if (!(logfc_col %in% colnames(re$summary)) || !(p_col %in% colnames(re$summary))) {
  stop(sprintf("nebula result is missing expected columns: %s, %s", logfc_col, p_col))
}

results_df <- data.frame(
  gene = re$summary$gene,
  logFC = re$summary[[logfc_col]],
  p_val = re$summary[[p_col]]
)

write.table(results_df, file = output_file, sep = '\t', row.names = FALSE, quote = FALSE)
