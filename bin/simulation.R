#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(splatter)
  library(anndataR)
  library(argparser)
  library(SingleCellExperiment)
})

p <- arg_parser("Splatter population simulation")
p <- add_argument(p, "--output", help = "Output h5ad", type = "character")
p <- add_argument(p, "--design_output", help = "Output sample design TSV", type = "character")
p <- add_argument(p, "--scenario", help = "Simulation scenario", type = "character")
p <- add_argument(p, "--n_genes", help = "Number of genes", type = "integer")
p <- add_argument(p, "--seed", help = "Requested seed", type = "integer")
argv <- parse_args(p)

lfc_values <- c("1" = 1.0, "1p5" = 1.5, "2" = 2.0)

parse_scenario <- function(scenario) {
  regular <- regmatches(scenario, regexec("^dataset_n(10|20|30)_lfc(1|1p5|2)$", scenario))[[1]]
  if (length(regular) > 0) {
    n_samples <- as.integer(regular[2])
    return(list(
      family = "sample_size", n_samples = n_samples,
      batch_sizes = rep(10L, n_samples / 10L), condition_prob = c(0.5, 0.5),
      cells_per_sample = 300, sample_cells = FALSE,
      target_mean_cells = 300, target_cv_cells = 0,
      lfc = unname(lfc_values[regular[3]]), custom_batches = FALSE
    ))
  }

  sample_size <- regmatches(scenario, regexec("^review_n(18|22)_lfc(1|1p5|2)$", scenario))[[1]]
  if (length(sample_size) > 0) {
    n_samples <- as.integer(sample_size[2])
    return(list(
      family = "sample_size", n_samples = n_samples,
      batch_sizes = rep(as.integer(n_samples / 2), 2), condition_prob = c(0.5, 0.5),
      cells_per_sample = 300, sample_cells = FALSE,
      target_mean_cells = 300, target_cv_cells = 0,
      lfc = unname(lfc_values[sample_size[3]]), custom_batches = FALSE
    ))
  }

  cell_cv <- regmatches(scenario, regexec("^review_cv(0p5|1|1p5)_lfc1p5$", scenario))[[1]]
  if (length(cell_cv) > 0) {
    cv <- c("0p5" = 0.5, "1" = 1.0, "1p5" = 1.5)[cell_cv[2]]
    return(list(
      family = "cell_imbalance", n_samples = 30L,
      batch_sizes = rep(10L, 3), condition_prob = c(0.5, 0.5),
      cells_per_sample = 300, sample_cells = TRUE,
      n_cells_shape = 1 / cv^2, n_cells_rate = (1 / cv^2) / 300,
      target_mean_cells = 300, target_cv_cells = unname(cv),
      lfc = 1.5, custom_batches = FALSE
    ))
  }

  cell_level <- regmatches(scenario, regexec("^review_m(5|10|25|50|100|300)_lfc1p5$", scenario))[[1]]
  if (length(cell_level) > 0) {
    cells <- as.integer(cell_level[2])
    return(list(
      family = "cells_per_subject", n_samples = 20L,
      batch_sizes = rep(10L, 2), condition_prob = c(0.5, 0.5),
      cells_per_sample = cells, sample_cells = FALSE,
      target_mean_cells = cells, target_cv_cells = 0,
      lfc = 1.5, custom_batches = FALSE
    ))
  }

  if (scenario == "review_stress_lfc1p5") {
    return(list(
      family = "combined_imbalance", n_samples = 30L,
      batch_sizes = c(6L, 10L, 14L), condition_prob = c(1 / 3, 2 / 3),
      cells_per_sample = 300, sample_cells = TRUE,
      n_cells_shape = 1, n_cells_rate = 1 / 300,
      target_mean_cells = 300, target_cv_cells = 1,
      lfc = 1.5, custom_batches = TRUE
    ))
  }

  stop(sprintf("Unknown simulation scenario: %s", scenario))
}

build_params <- function(config) {
  n_batches <- length(config$batch_sizes)
  args <- list(
    nGenes = argv$n_genes,
    eqtl.n = 0,
    batch.size = config$batch_sizes[1],
    batchCells = rep(config$cells_per_sample, n_batches),
    condition.prob = config$condition_prob,
    cde.prob = 0.1,
    cde.downProb = 0.5,
    cde.facLoc = config$lfc * log(2),
    cde.facScale = 0.3
  )
  if (config$sample_cells) {
    args$nCells.sample <- TRUE
    args$nCells.shape <- config$n_cells_shape
    args$nCells.rate <- config$n_cells_rate
  }
  if (n_batches > 1) {
    args$batch.facLoc <- seq(0, by = 0.5, length.out = n_batches)
    args$batch.facScale <- rep(0.1, n_batches)
  }
  params <- do.call(newSplatPopParams, args)
  setParams(params, update = list(similarity.scale = 0.3))
}

clean_for_batch_resimulation <- function(sim) {
  for (assay_name in intersect(
    c("counts", "TrueCounts", "CellMeans", "BCV", "BaseCellMeans", "BatchCellMeans"),
    assayNames(sim)
  )) {
    assay(sim, assay_name) <- NULL
  }
  drop_rowdata <- intersect(
    c(
      "chromosome", "geneStart", "geneEnd", "geneMiddle",
      "eQTL.group", "eQTL.condition", "eSNP.ID", "eSNP.chromosome",
      "eSNP.loc", "eSNP.MAF", "eQTL.EffectSize"
    ),
    colnames(rowData(sim))
  )
  if (length(drop_rowdata) > 0) rowData(sim)[, drop_rowdata] <- NULL
  drop_coldata <- intersect(c("Batch", "Group"), colnames(colData(sim)))
  if (length(drop_coldata) > 0) colData(sim)[, drop_coldata] <- NULL
  sim
}

assign_custom_batches <- function(sim, batch_sizes) {
  samples <- unique(as.character(colData(sim)$Sample))
  stopifnot(sum(batch_sizes) == length(samples))
  shuffled <- sample(samples, length(samples), replace = FALSE)
  ends <- cumsum(batch_sizes)
  starts <- c(1L, head(ends, -1L) + 1L)
  sample_batches <- unlist(lapply(seq_along(batch_sizes), function(i) {
    setNames(rep(paste0("Batch", i), batch_sizes[i]), shuffled[starts[i]:ends[i]])
  }))
  colData(sim)$Batch <- unname(sample_batches[as.character(colData(sim)$Sample)])
  sim
}

get_lnorm_factors <- function(n_factors, selection_prob, negative_prob, location, scale) {
  selected <- as.logical(rbinom(n_factors, 1, selection_prob))
  directions <- (-1)^rbinom(sum(selected), 1, negative_prob)
  selected_factors <- rlnorm(sum(selected), location, scale)
  directions[selected_factors < 1] <- -directions[selected_factors < 1]
  factors <- rep(1, n_factors)
  factors[selected] <- selected_factors^directions
  factors
}

simulate_batch_effects <- function(sim, params) {
  n_genes <- nrow(sim)
  batch_fac_loc <- getParam(params, "batch.facLoc")
  batch_fac_scale <- getParam(params, "batch.facScale")
  samples <- unique(as.character(colData(sim)$Sample))

  for (sample_id in samples) {
    batch <- unique(as.character(colData(sim)$Batch[colData(sim)$Sample == sample_id]))
    batch_number <- as.integer(sub("^Batch", "", batch))
    rowData(sim)[[paste(sample_id, "BatchFac", sep = "_")]] <- get_lnorm_factors(
      n_genes, 1, 0.5, batch_fac_loc[batch_number], batch_fac_scale[batch_number]
    )
  }

  batch_means <- lapply(samples, function(sample_id) {
    gene_means <- metadata(sim)$Simulated_Means$Group1[, sample_id]
    n_cells <- sum(colData(sim)$Sample == sample_id)
    factors <- as.matrix(rowData(sim)[, paste(sample_id, "BatchFac", sep = "_")])
    as.matrix(factors[, rep(1, n_cells)]) * gene_means
  })
  batch_means <- do.call(cbind, batch_means)
  rownames(batch_means) <- rownames(sim)
  colnames(batch_means) <- colnames(sim)
  assays(sim)$BatchCellMeans <- batch_means

  cell_props <- t(t(batch_means) / colSums(batch_means))
  assays(sim)$BaseCellMeans <- t(t(cell_props) * colData(sim)$ExpLibSize)

  bcv_common <- getParam(params, "bcv.common")
  bcv_df <- getParam(params, "bcv.df")
  cell_means <- list()
  bcv_values <- list()
  for (sample_id in samples) {
    sample_cells <- colData(sim)$Sample == sample_id
    base_sample <- assays(sim[, sample_cells])$BaseCellMeans
    sample_bcv <- (bcv_common + (1 / sqrt(base_sample))) *
      sqrt(bcv_df / rchisq(n_genes, df = bcv_df))
    bcv_values[[sample_id]] <- sample_bcv
    cell_means[[sample_id]] <- matrix(
      rgamma(
        as.numeric(n_genes) * sum(sample_cells),
        shape = 1 / (sample_bcv^2),
        scale = base_sample * (sample_bcv^2)
      ),
      nrow = n_genes,
      ncol = sum(sample_cells)
    )
  }
  cell_means <- do.call(cbind, cell_means)
  rownames(cell_means) <- rownames(sim)
  colnames(cell_means) <- colnames(sim)
  assays(sim)$BCV <- do.call(cbind, bcv_values)
  assays(sim)$CellMeans <- cell_means

  true_counts <- matrix(
    rpois(as.numeric(n_genes) * ncol(sim), lambda = cell_means),
    nrow = n_genes, ncol = ncol(sim), dimnames = dimnames(sim)
  )
  assays(sim)$TrueCounts <- true_counts
  assays(sim)$counts <- true_counts
  sim
}

simulate_once <- function(config, params) {
  vcf <- mockVCF(n.samples = config$n_samples)
  if (!config$custom_batches) {
    return(splatPopSimulate(vcf = vcf, params = params, verbose = FALSE))
  }

  internal <- params
  internal <- setParam(internal, "batch.size", 1)
  internal <- setParam(internal, "batch.facLoc", 0)
  internal <- setParam(internal, "batch.facScale", 0)
  internal <- setParam(internal, "batchCells", mean(getParam(params, "batchCells")))
  sim <- splatPopSimulate(vcf = vcf, params = internal, verbose = FALSE)
  sim <- clean_for_batch_resimulation(sim)
  sim <- assign_custom_batches(sim, config$batch_sizes)
  simulate_batch_effects(sim, params)
}

get_sample_design <- function(sim) {
  design <- unique(as.data.frame(colData(sim)[, c("Sample", "Batch", "Condition")]))
  design$Sample <- as.character(design$Sample)
  design$Batch <- as.character(design$Batch)
  design$Condition <- as.character(design$Condition)
  design
}

has_valid_design <- function(sim) {
  tab <- with(get_sample_design(sim), table(Batch, Condition))
  ncol(tab) == 2 && !any(tab == 0)
}

write_design <- function(sim, config, requested_seed, used_seed, output) {
  design <- get_sample_design(sim)
  sample_cells <- as.data.frame(table(as.character(colData(sim)$Sample)), stringsAsFactors = FALSE)
  colnames(sample_cells) <- c("Sample", "sample_cells")
  design <- merge(design, sample_cells, by = "Sample", all.x = TRUE, sort = FALSE)

  condition_levels <- sort(unique(design$Condition))
  condition_counts <- table(factor(design$Condition, levels = condition_levels))
  batch_condition <- as.data.frame(
    table(design$Batch, factor(design$Condition, levels = condition_levels)),
    stringsAsFactors = FALSE
  )
  colnames(batch_condition) <- c("Batch", "Condition", "n")
  batch_total <- aggregate(n ~ Batch, batch_condition, sum)
  colnames(batch_total)[2] <- "batch_n"
  condition_count_for_batch <- function(batch, condition) {
    value <- batch_condition$n[
      batch_condition$Batch == batch & batch_condition$Condition == condition
    ]
    if (length(value) == 0) 0L else as.integer(value[1])
  }

  actual_mean <- mean(design$sample_cells)
  actual_cv <- if (actual_mean > 0) sd(design$sample_cells) / actual_mean else NA_real_
  design$scenario <- argv$scenario
  design$requested_seed <- requested_seed
  design$used_seed <- used_seed
  design$n_total <- nrow(design)
  design$n_case <- as.integer(condition_counts[1])
  design$n_control <- as.integer(condition_counts[2])
  design$batch_n <- batch_total$batch_n[match(design$Batch, batch_total$Batch)]
  design$batch_case <- vapply(
    design$Batch, condition_count_for_batch, integer(1), condition = condition_levels[1]
  )
  design$batch_control <- vapply(
    design$Batch, condition_count_for_batch, integer(1), condition = condition_levels[2]
  )
  design$target_case_prob <- config$condition_prob[1]
  design$target_control_prob <- config$condition_prob[2]
  design$target_mean_cells <- config$target_mean_cells
  design$target_cv_cells <- config$target_cv_cells
  design$actual_mean_cells <- actual_mean
  design$actual_cv_cells <- actual_cv
  design$total_cells <- sum(design$sample_cells)
  design$lfc <- config$lfc

  design <- design[, c(
    "scenario", "requested_seed", "used_seed", "Sample", "Batch", "Condition",
    "sample_cells", "n_total", "n_case", "n_control", "batch_n",
    "batch_case", "batch_control", "target_case_prob",
    "target_control_prob", "target_mean_cells", "target_cv_cells",
    "actual_mean_cells", "actual_cv_cells", "total_cells", "lfc"
  )]
  write.table(design, output, sep = "\t", quote = FALSE, row.names = FALSE)
}

config <- parse_scenario(argv$scenario)
params <- build_params(config)
current_seed <- argv$seed
repeat {
  set.seed(current_seed)
  sim <- simulate_once(config, params)
  if (has_valid_design(sim)) break
  message("Seed ", current_seed, " produced a confounded design; retrying")
  current_seed <- current_seed + 1L
}

gene_ids <- rownames(sim)
canonical_gene_ids <- gsub("_", "-", gene_ids, fixed = TRUE)
stopifnot(anyDuplicated(canonical_gene_ids) == 0)
rownames(sim) <- canonical_gene_ids
rowData(sim)$orig_gene_id <- gene_ids
rowData(sim)$Gene <- canonical_gene_ids

drop_rowdata <- intersect(
  c(
    "chromosome", "geneStart", "geneEnd", "geneMiddle",
    "eQTL.group", "eQTL.condition", "eSNP.ID", "eSNP.chromosome",
    "eSNP.loc", "eSNP.MAF", "eQTL.EffectSize"
  ),
  colnames(rowData(sim))
)
if (length(drop_rowdata) > 0) rowData(sim)[, drop_rowdata] <- NULL
for (assay_name in intersect(
  c("BCV", "BaseCellMeans", "BatchCellMeans", "CellMeans", "TrueCounts"),
  assayNames(sim)
)) {
  assay(sim, assay_name) <- NULL
}

write_design(sim, config, argv$seed, current_seed, argv$design_output)
metadata(sim) <- list()
gc()
write_h5ad(sim, argv$output, x_mapping = "counts")
