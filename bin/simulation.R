#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(splatter)
  library(anndataR)
  library(argparser)
  library(SingleCellExperiment)
})

p <- arg_parser("Splatter simulation")
p <- add_argument(p, "--output", help = "output h5ad", type = "character")
p <- add_argument(p, "--scenario", help = "simulation scenario", type = "character")
p <- add_argument(p, "--n_genes", help = "number of genes", type = "integer")
p <- add_argument(p, "--seed", help = "seed", type = "integer")

argv <- parse_args(p)

scenario_match <- regmatches(argv$scenario, regexec("^dataset_n(10|20|30)_lfc(1|1p5|2)$", argv$scenario))[[1]]
n_samples <- as.integer(scenario_match[2])
lfc_loc <- c("1" = 1.0, "1p5" = 1.5, "2" = 2.0)[scenario_match[3]]
n_batches <- n_samples / 10

params_args <- list(
  nGenes = argv$n_genes,
  eqtl.n = 0,
  batch.size = 10,
  batchCells = rep(300, n_batches),
  condition.prob = c(0.5, 0.5),
  cde.prob = 0.1,
  cde.downProb = 0.5,
  cde.facLoc = unname(lfc_loc) * log(2),
  cde.facScale = 0.3
)

if (n_batches > 1) {
  params_args$batch.facLoc <- seq(0, by = 0.5, length.out = n_batches)
  params_args$batch.facScale <- rep(0.1, n_batches)
}

params <- do.call(newSplatPopParams, params_args)
params <- setParams(params, update = list(similarity.scale = 0.3))
vcf <- mockVCF(n.samples = n_samples)

check_design <- function(sim) {
  design <- unique(as.data.frame(colData(sim)[, c("Sample", "Batch", "Condition")]))
  !any(table(design$Batch, design$Condition) == 0)
}

current_seed <- argv$seed
repeat {
  set.seed(current_seed)
  sim <- splatPopSimulate(vcf = vcf, params = params, verbose = FALSE)
  if (check_design(sim)) {
    break
  }
  current_seed <- current_seed + 1
}

gene_ids <- rownames(sim)
canonical_gene_ids <- gsub("_", "-", gene_ids, fixed = TRUE)
rownames(sim) <- canonical_gene_ids
rowData(sim)$orig_gene_id <- gene_ids
rowData(sim)$Gene <- canonical_gene_ids

rowData(sim)[, c(
  "chromosome", "geneStart", "geneEnd", "geneMiddle",
  "eQTL.group", "eQTL.condition", "eSNP.ID", "eSNP.chromosome",
  "eSNP.loc", "eSNP.MAF", "eQTL.EffectSize"
)] <- NULL

assays(sim)$BCV <- NULL
assays(sim)$BaseCellMeans <- NULL
assays(sim)$BatchCellMeans <- NULL
assays(sim)$CellMeans <- NULL
assays(sim)$TrueCounts <- NULL

metadata(sim) <- list()
gc()

write_h5ad(sim, argv$output, x_mapping = "counts")
