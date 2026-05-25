#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(splatter)
  library(anndataR)
  library(SingleCellExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: sim-runtime.R <output.h5ad> <n_genes> <n_cells> <seed>")
}

output <- args[[1]]
n_genes <- as.integer(args[[2]])
n_cells <- as.integer(args[[3]])
seed <- as.integer(args[[4]])
if (anyNA(c(n_genes, n_cells, seed))) {
  stop("<n_genes>, <n_cells>, and <seed> must be integers")
}

make_simulation_inputs <- function(n_genes, n_cells) {
  list(
    vcf = mockVCF(n.samples = 10),
    params = newSplatPopParams(
      nGenes = n_genes,
      eqtl.ES.shape = 0,
      eqtl.ES.rate = 0,
      out.prob = 0,
      eqtl.n = 0,
      eqtl.group.specific = 0,
      eqtl.condition.specific = 0,
      batch.size = 10,
      batchCells = n_cells / 10,
      condition.prob = c(0.5, 0.5),
      cde.prob = 0.025,
      cde.downProb = 0.5,
      cde.facLoc = 2.5,
      cde.facScale = 0.4
    )
  )
}

get_design <- function(sim) {
  unique(as.data.frame(colData(sim)[, c("Sample", "Batch", "Condition")]))
}

has_balanced_design <- function(sim) {
  design_table <- with(get_design(sim), table(Batch, Condition))
  !any(design_table == 0)
}

simulate_valid_design <- function(inputs, seed) {
  repeat {
    set.seed(seed)
    sim <- splatPopSimulate(vcf = inputs$vcf, params = inputs$params, verbose = FALSE)
    if (has_balanced_design(sim)) {
      return(list(sim = sim, seed = seed))
    }

    message("Seed ", seed, " produced a confounded design. Retrying with next seed...")
    seed <- seed + 1
  }
}

clean_simulation <- function(sim) {
  orig_gene_ids <- rownames(sim)
  canonical_gene_ids <- gsub("_", "-", orig_gene_ids, fixed = TRUE)
  stopifnot(anyDuplicated(canonical_gene_ids) == 0)

  rownames(sim) <- canonical_gene_ids
  rowData(sim)$orig_gene_id <- orig_gene_ids
  if ("Gene" %in% colnames(rowData(sim))) {
    rowData(sim)$Gene <- canonical_gene_ids
  }

  drop_rowdata <- intersect(
    c(
      "chromosome", "geneStart", "geneEnd", "geneMiddle",
      "eQTL.group", "eQTL.condition", "eSNP.ID", "eSNP.chromosome",
      "eSNP.loc", "eSNP.MAF", "eQTL.EffectSize"
    ),
    colnames(rowData(sim))
  )
  if (length(drop_rowdata) > 0) {
    rowData(sim)[, drop_rowdata] <- NULL
  }

  drop_assays <- intersect(
    c("BCV", "BaseCellMeans", "BatchCellMeans", "CellMeans", "TrueCounts"),
    assayNames(sim)
  )
  if (length(drop_assays) > 0) {
    assays(sim) <- assays(sim)[setdiff(assayNames(sim), drop_assays)]
  }

  metadata(sim) <- list()
  sim
}

message("Running runtime simulation with base seed: ", seed)
result <- simulate_valid_design(make_simulation_inputs(n_genes, n_cells), seed)
sim <- clean_simulation(result$sim)

message("Valid experimental design generated successfully with final seed: ", result$seed)
message("Sample distribution across batches and conditions:")
print(with(get_design(sim), table(Batch, Condition)))

gc()
write_h5ad(sim, output, x_mapping = "counts")
