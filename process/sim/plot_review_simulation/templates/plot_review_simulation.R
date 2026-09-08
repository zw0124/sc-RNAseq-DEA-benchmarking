#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) stop(sprintf("Missing argument: %s", flag))
  args[index + 1]
}

metrics_dir <- get_arg("--metrics-dir")
design_dir <- get_arg("--design-dir")
fpr_dir <- get_arg("--fpr-dir")
target_thresholds <- c(0.01, 0.05, 0.10)

method_order <- c(
  "ttest", "wilcox", "edger", "limma", "deseq2",
  "nebula", "memento", "BOOTSTRAP3", "BOOTSTRAP2"
)
method_labels <- c(
  ttest = "T-test", wilcox = "Wilcoxon", edger = "edgeR", limma = "limma",
  deseq2 = "DESeq2", nebula = "NEBULA", memento = "memento",
  BOOTSTRAP3 = "BOOTSTRAP3", BOOTSTRAP2 = "BOOTSTRAP2"
)

method_groups <- list(
  Wilcoxon = c("ttest", "wilcox"),
  Pseudobulk = c("edger", "limma", "deseq2"),
  Mixed = "nebula",
  Deconvolution = "memento",
  Bootstrapping = c("BOOTSTRAP3", "BOOTSTRAP2")
)
base_colors <- c(
  Wilcoxon = "#000000", Pseudobulk = "#ff0000", Mixed = "#fdaa05",
  Deconvolution = "#0080ff", Bootstrapping = "#00aa55"
)

build_method_colors <- function(methods) {
  colors <- character(0)
  for (group in names(method_groups)) {
    selected <- intersect(method_groups[[group]], methods)
    if (length(selected) == 0) next
    if (length(selected) == 1) {
      values <- base_colors[group]
    } else {
      light <- colorRampPalette(c("white", base_colors[group]))(10)[4]
      dark <- colorRampPalette(c(base_colors[group], "black"))(10)[4]
      values <- colorRampPalette(c(light, dark))(length(selected))
    }
    names(values) <- selected
    colors <- c(colors, values)
  }
  colors
}

parse_metric_file <- function(path) {
  stem <- sub("\\.tsv$", "", basename(path))
  parts <- strsplit(stem, "_", fixed = TRUE)[[1]]
  if (length(parts) < 5 || !identical(parts[1:2], c("fdr", "power"))) return(NULL)
  method <- parts[length(parts)]
  run <- suppressWarnings(as.integer(parts[length(parts) - 1]))
  scenario <- paste(parts[3:(length(parts) - 2)], collapse = "_")
  dt <- tryCatch(fread(path), error = function(e) NULL)
  required <- c("p_type", "threshold", "observed_fdr", "power")
  if (is.null(dt) || !all(required %in% names(dt))) return(NULL)
  dt <- dt[p_type == "adj_target_grid"]
  if (nrow(dt) == 0) return(NULL)
  dt[, `:=`(scenario = scenario, run = run, method = method)]
  dt[, .(scenario, run, method, threshold, observed_fdr, power, discoveries)]
}

describe_scenario <- function(scenario) {
  regular <- regmatches(scenario, regexec("^dataset_n(10|20|30)_lfc(1|1p5|2)$", scenario))[[1]]
  if (length(regular) > 0) {
    n <- as.integer(regular[2])
    lfc <- c("1" = 1, "1p5" = 1.5, "2" = 2)[regular[3]]
    rows <- data.table(family = "sample_size", setting = paste0("N=", n), setting_value = n, lfc = lfc)
    if (n == 30 && lfc == 1.5) {
      rows <- rbind(
        rows,
        data.table(family = "cell_imbalance", setting = "CV=0", setting_value = 0, lfc = lfc)
      )
    }
    return(rows)
  }

  sample_size <- regmatches(scenario, regexec("^review_n(18|22)_lfc(1|1p5|2)$", scenario))[[1]]
  if (length(sample_size) > 0) {
    n <- as.integer(sample_size[2])
    lfc <- c("1" = 1, "1p5" = 1.5, "2" = 2)[sample_size[3]]
    return(data.table(family = "sample_size", setting = paste0("N=", n), setting_value = n, lfc = lfc))
  }

  cell_cv <- regmatches(scenario, regexec("^review_cv(0p5|1|1p5)_lfc1p5$", scenario))[[1]]
  if (length(cell_cv) > 0) {
    cv <- c("0p5" = 0.5, "1" = 1, "1p5" = 1.5)[cell_cv[2]]
    return(data.table(family = "cell_imbalance", setting = paste0("CV=", cv), setting_value = cv, lfc = 1.5))
  }

  cell_level <- regmatches(scenario, regexec("^review_m(5|10|25|50|100|300)_lfc1p5$", scenario))[[1]]
  if (length(cell_level) > 0) {
    cells <- as.integer(cell_level[2])
    return(data.table(
      family = "cells_per_subject", setting = paste0("M=", cells),
      setting_value = cells, lfc = 1.5
    ))
  }

  if (scenario == "review_stress_lfc1p5") {
    return(data.table(
      family = "combined_imbalance", setting = "Combined imbalance",
      setting_value = 1, lfc = 1.5
    ))
  }
  NULL
}

metric_files <- list.files(metrics_dir, pattern = "^fdr_power_.*\\.tsv$", full.names = TRUE)
metric_rows <- lapply(metric_files, parse_metric_file)
metric_rows <- metric_rows[!vapply(metric_rows, is.null, logical(1))]
metrics <- if (length(metric_rows) > 0) rbindlist(metric_rows, fill = TRUE) else data.table()
if (nrow(metrics) == 0) stop("No usable FDR/power files were found")

descriptions <- lapply(unique(metrics$scenario), function(scenario) {
  result <- describe_scenario(scenario)
  if (is.null(result)) return(NULL)
  result[, scenario := scenario]
  result
})
descriptions <- descriptions[!vapply(descriptions, is.null, logical(1))]
description_dt <- rbindlist(descriptions, fill = TRUE)
runs <- merge(metrics, description_dt, by = "scenario", allow.cartesian = TRUE)
runs[, threshold := round(as.numeric(threshold), 2)]
runs[, method := factor(method, levels = method_order)]
runs <- runs[!is.na(method)]
setorder(runs, family, setting_value, lfc, threshold, method, run)
fwrite(runs, "review_simulation_runs.tsv", sep = "\t")
family_outputs <- c(
  sample_size = "review_sample_size_runs.tsv",
  cell_imbalance = "review_cell_imbalance_runs.tsv",
  cells_per_subject = "review_cells_per_subject_runs.tsv",
  combined_imbalance = "review_combined_imbalance_runs.tsv"
)
for (family_name in names(family_outputs)) {
  fwrite(runs[family == family_name], family_outputs[[family_name]], sep = "\t")
}

long_runs <- melt(
  runs,
  id.vars = c("scenario", "run", "method", "threshold", "family", "setting", "setting_value", "lfc"),
  measure.vars = c("observed_fdr", "power"),
  variable.name = "metric",
  value.name = "value"
)
summary_dt <- long_runs[, .(
  mean = mean(value, na.rm = TRUE),
  sd = if (sum(!is.na(value)) > 1) sd(value, na.rm = TRUE) else NA_real_,
  n_runs = sum(!is.na(value))
), by = .(family, setting, setting_value, lfc, threshold, method, metric)]
setorder(summary_dt, family, setting_value, lfc, threshold, method, metric)
fwrite(summary_dt, "review_simulation_summary.tsv", sep = "\t")

design_files <- list.files(design_dir, pattern = "^design_.*\\.tsv$", full.names = TRUE)
design_rows <- lapply(design_files, function(path) {
  dt <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  dt
})
design_rows <- design_rows[!vapply(design_rows, is.null, logical(1))]
design_dt <- if (length(design_rows) > 0) rbindlist(design_rows, fill = TRUE) else data.table()
fwrite(design_dt, "review_simulation_design.tsv", sep = "\t")

parse_fpr_file <- function(path) {
  stem <- sub("\\.tsv$", "", basename(path))
  parts <- strsplit(stem, "_", fixed = TRUE)[[1]]
  if (length(parts) < 4 || parts[1] != "fpr") return(NULL)
  method <- parts[length(parts)]
  run <- suppressWarnings(as.integer(parts[length(parts) - 1]))
  scenario <- paste(parts[2:(length(parts) - 2)], collapse = "_")
  if (!(method %in% method_order) || is.na(run)) return(NULL)
  dt <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(dt) || !all(c("raw_p_value", "fpr") %in% names(dt))) return(NULL)
  dt[, `:=`(scenario = scenario, run = run, method = method)]
  dt[, .(scenario, run, method, raw_p_value, fpr)]
}

fpr_files <- list.files(fpr_dir, pattern = "^fpr_.*\\.tsv$", full.names = TRUE)
fpr_rows <- lapply(fpr_files, parse_fpr_file)
fpr_rows <- fpr_rows[!vapply(fpr_rows, is.null, logical(1))]
fpr_runs <- if (length(fpr_rows) > 0) rbindlist(fpr_rows, fill = TRUE) else data.table(
  scenario = character(), run = integer(), method = character(),
  raw_p_value = numeric(), fpr = numeric()
)
if (nrow(fpr_runs) > 0) {
  fpr_runs <- merge(fpr_runs, description_dt, by = "scenario", allow.cartesian = TRUE)
  fpr_runs[, method := factor(method, levels = method_order)]
  fpr_runs <- fpr_runs[!is.na(method)]
  setorder(fpr_runs, family, setting_value, lfc, raw_p_value, method, run)
  fpr_summary <- fpr_runs[, .(
    mean_fpr = mean(fpr, na.rm = TRUE),
    sd_fpr = if (sum(!is.na(fpr)) > 1) sd(fpr, na.rm = TRUE) else NA_real_,
    n_runs = uniqueN(run[!is.na(fpr)])
  ), by = .(family, setting, setting_value, lfc, raw_p_value, method)]
} else {
  fpr_summary <- data.table(
    family = character(), setting = character(), setting_value = numeric(),
    lfc = numeric(), raw_p_value = numeric(), method = character(),
    mean_fpr = numeric(), sd_fpr = numeric(), n_runs = integer()
  )
}
fwrite(fpr_runs, "review_simulation_fpr_runs.tsv", sep = "\t")
fwrite(fpr_summary, "review_simulation_fpr_summary.tsv", sep = "\t")

available_methods <- intersect(method_order, as.character(unique(na.omit(runs$method))))
method_colors <- build_method_colors(available_methods)
bootstrap_colors <- c(BOOTSTRAP3 = "#7B2CBF", BOOTSTRAP2 = "#008837")
bootstrap_methods <- intersect(names(bootstrap_colors), names(method_colors))
method_colors[bootstrap_methods] <- bootstrap_colors[bootstrap_methods]
method_shapes <- c(
  ttest = 16, wilcox = 17, edger = 15, limma = 18,
  deseq2 = 3, nebula = 8, memento = 1,
  BOOTSTRAP3 = 2, BOOTSTRAP2 = 0
)
method_linetypes <- c(
  ttest = "solid", wilcox = "dashed", edger = "dotted", limma = "dotdash",
  deseq2 = "longdash", nebula = "solid", memento = "dashed",
  BOOTSTRAP3 = "dotted", BOOTSTRAP2 = "dotdash"
)

empty_plot <- function(path, label) {
  plot <- ggplot() + annotate("text", x = 0, y = 0, label = label, size = 5) + theme_void()
  ggsave(path, plot, width = 8, height = 5, dpi = 300, bg = "white")
}

save_distribution <- function(
  family_name, output, x_label, width, height,
  thresholds = target_thresholds
) {
  dt <- copy(long_runs[family == family_name & threshold %in% thresholds])
  if (nrow(dt) == 0) {
    empty_plot(output, paste("No data for", family_name))
    return(invisible(NULL))
  }
  plot_methods <- intersect(available_methods, as.character(unique(dt$method)))
  setting_order <- unique(dt[order(setting_value)]$setting)
  dt[, setting := factor(setting, levels = setting_order)]
  dt[, method := factor(as.character(method), levels = plot_methods)]
  dt[, metric := factor(metric, levels = c("observed_fdr", "power"), labels = c("Observed FDR", "Power"))]
  dt[, threshold_label := factor(
    sprintf("Nominal FDR = %.2f", threshold),
    levels = sprintf("Nominal FDR = %.2f", thresholds)
  )]

  if (family_name == "sample_size") {
    dt[, lfc_label := factor(paste0("LFC = ", lfc), levels = paste0("LFC = ", c(1, 1.5, 2)))]
    facets <- facet_grid(metric + lfc_label ~ threshold_label, scales = "free_y")
  } else {
    facets <- facet_grid(metric ~ threshold_label, scales = "free_y")
  }

  plot <- ggplot(dt, aes(x = setting, y = value, color = method, fill = method)) +
    geom_boxplot(
      position = position_dodge(width = 0.82), width = 0.72,
      outlier.shape = NA, alpha = 0.18, linewidth = 0.35
    ) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.82),
      size = 0.8, alpha = 0.65
    ) +
    facets +
    scale_color_manual(values = method_colors[plot_methods], labels = method_labels[plot_methods], drop = FALSE) +
    scale_fill_manual(values = method_colors[plot_methods], labels = method_labels[plot_methods], drop = FALSE) +
    labs(x = x_label, y = NULL, color = "Method", fill = "Method") +
    theme_cowplot(12) +
    theme(
      legend.position = "right",
      strip.background = element_rect(fill = "grey15", color = "grey30"),
      strip.text = element_text(color = "white", face = "bold"),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )

  if (family_name %in% c("sample_size", "cell_imbalance", "cells_per_subject", "combined_imbalance")) {
    reference_lines <- data.table(
      metric = factor("Observed FDR", levels = levels(dt$metric)),
      threshold_label = factor(
        sprintf("Nominal FDR = %.2f", thresholds),
        levels = levels(dt$threshold_label)
      ),
      nominal_fdr = thresholds
    )
    plot <- plot +
      geom_hline(
        data = reference_lines,
        aes(yintercept = nominal_fdr),
        inherit.aes = FALSE,
        color = "#d62728",
        linetype = "dashed",
        linewidth = 0.6
      )
  }

  if (family_name == "combined_imbalance") {
    plot <- plot +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
  ggsave(output, plot, width = width, height = height, dpi = 300, bg = "white")
}

figure1_runs <- runs[
  family == "sample_size" &
    setting_value %in% c(10, 18, 20, 22, 30) &
    threshold %in% c(0.01, 0.05, 0.10)
]
if (nrow(figure1_runs) == 0) {
  empty_plot("reviewer_figure1.png", "No Figure 1 data")
  fwrite(data.table(
    method = character(), lfc = numeric(), setting_value = numeric(), threshold = numeric(),
    observed_fdr = numeric(), power = numeric(), observed_fdr_sd = numeric(),
    power_sd = numeric(), n_runs = integer(), lfc_label = character(), n_label = character()
  ), "reviewer_figure1_summary.tsv", sep = "\t")
} else {
  reviewer_agg <- figure1_runs[, .(
    observed_fdr = mean(observed_fdr, na.rm = TRUE),
    power = mean(power, na.rm = TRUE),
    observed_fdr_sd = if (sum(!is.na(observed_fdr)) > 1) sd(observed_fdr, na.rm = TRUE) else NA_real_,
    power_sd = if (sum(!is.na(power)) > 1) sd(power, na.rm = TRUE) else NA_real_,
    n_runs = uniqueN(run)
  ), by = .(method, lfc, setting_value, threshold)]
  reviewer_agg[, method := factor(method, levels = available_methods)]
  reviewer_agg[, lfc_label := factor(
    paste0("LFC=", lfc),
    levels = paste0("LFC=", c(1, 1.5, 2))
  )]
  reviewer_agg[, n_label := factor(
    paste0("n=", setting_value),
    levels = paste0("n=", c(10, 18, 20, 22, 30))
  )]
  setorder(reviewer_agg, lfc_label, n_label, method, threshold)
  fwrite(reviewer_agg, "reviewer_figure1_summary.tsv", sep = "\t")

  plot <- ggplot(
    reviewer_agg,
    aes(observed_fdr, power, color = method, shape = method, linetype = method)
  ) +
    geom_vline(xintercept = 0.1, linetype = "dashed", color = "#b2182b", linewidth = 0.6) +
    geom_path(aes(group = method), linewidth = 0.75, alpha = 0.82) +
    geom_point(size = 2.0, stroke = 0.8) +
    facet_grid(rows = vars(lfc_label), cols = vars(n_label), scales = "free") +
    scale_color_manual(values = method_colors, labels = method_labels[available_methods], drop = FALSE) +
    scale_shape_manual(values = method_shapes[available_methods], labels = method_labels[available_methods], drop = FALSE) +
    scale_linetype_manual(values = method_linetypes[available_methods], labels = method_labels[available_methods], drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.10))) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.10))) +
    labs(
      x = "Observed FDR", y = "Power",
      color = "Method", shape = "Method", linetype = "Method"
    ) +
    theme_cowplot(13) +
    theme(
      legend.position = "right",
      strip.background = element_rect(fill = "grey15", color = "grey30"),
      strip.text = element_text(color = "white", face = "bold", size = 11),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 9)
    )
  ggsave("reviewer_figure1.png", plot, width = 18, height = 10, dpi = 300, bg = "white")
}

for (threshold in target_thresholds) {
  threshold_id <- sub("\\.", "p", sprintf("%.2f", threshold))
  save_distribution(
    "sample_size",
    paste0("review_sample_size_fdr", threshold_id, ".png"),
    "Total subjects",
    12,
    16,
    thresholds = threshold
  )
}
save_distribution("cell_imbalance", "review_cell_imbalance.png", "Cell-count CV", 12, 8)
save_distribution("cells_per_subject", "review_cells_per_subject.png", "Cells per subject", 12, 8)
save_distribution("combined_imbalance", "review_combined_imbalance.png", NULL, 11, 8)
