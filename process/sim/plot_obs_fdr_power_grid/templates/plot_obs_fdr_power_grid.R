#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(cowplot)

show_error_bars <- FALSE
adj_target_thresholds <- c(0.01, 0.05, 0.10)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || args[1] != "--input_dir") {
  stop("Usage: plot_obs_fdr_power_grid.R --input_dir <directory_containing_tsv_files>")
}

input_dir <- args[2]
if (!dir.exists(input_dir)) {
  stop(paste("Input directory does not exist:", input_dir))
}

fdr_power_files <- list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE)


paths_df <- data.frame(scenario = character(0), run = character(0), method = character(0), path = character(0))
for (p in fdr_power_files) {
  filename <- basename(p)
  filename <- gsub(".tsv$", "", filename)
  parts <- strsplit(filename, "_")[[1]]

  len <- length(parts)
  if (len < 5 || parts[1] != "fdr" || parts[2] != "power") {
    next
  }

  method <- parts[len]
  run <- parts[len - 1]
  scenario <- paste(parts[3:(len - 2)], collapse = "_")

  paths_df <- rbind(paths_df, data.frame(scenario = scenario, run = run, method = method, path = p))
}

extract_design <- function(scenario) {
  matched <- regmatches(scenario, regexec("^dataset_n(10|20|30)_lfc(1|1p5|2)$", scenario))[[1]]
  if (length(matched) == 0) {
    return(list(valid = FALSE, n = NA_character_, lfc = NA_character_))
  }

  lfc_label <- switch(
    matched[3],
    "1" = "1",
    "1p5" = "1.5",
    "2" = "2",
    NA_character_
  )

  list(valid = TRUE, n = matched[2], lfc = lfc_label)
}


if (nrow(paths_df) == 0) {
  png("Fig_ObsFDR_Power_Grid_empty.png", width = 1800, height = 1400)
  plot.new()
  text(0.5, 0.5, "No fdr_power files found", cex = 1.3)
  dev.off()
  quit(save = "no")
}

data_all <- data.table(
  p_type = character(),
  threshold = numeric(),
  observed_fdr = numeric(),
  power = numeric(),
  method = character(),
  scenario = character(),
  run = character()
)

for (i in seq_len(nrow(paths_df))) {
  dt <- fread(paths_df[i, ]$path)

  if (!all(c("p_type", "threshold", "observed_fdr", "power") %in% names(dt))) {
    next
  }

  dt <- dt[p_type == "adj_target_grid"]
  dt <- dt[round(threshold, 2) %in% round(adj_target_thresholds, 2)]
  if (nrow(dt) == 0) {
    next
  }

  dt$method <- paths_df[i, ]$method
  dt$scenario <- paths_df[i, ]$scenario
  dt$run <- paths_df[i, ]$run

  keep_cols <- c("p_type", "threshold", "observed_fdr", "power", "method", "scenario", "run")
  data_all <- rbind(data_all, dt[, ..keep_cols], fill = TRUE)
}

if (nrow(data_all) == 0) {
  png("Fig_ObsFDR_Power_Grid_empty.png", width = 1800, height = 1400)
  plot.new()
  text(0.5, 0.5, "No adj_target_grid rows found", cex = 1.3)
  dev.off()
  quit(save = "no")
}

design_info <- lapply(unique(data_all$scenario), extract_design)
design_dt <- data.table(
  scenario = unique(data_all$scenario),
  valid = vapply(design_info, function(x) x$valid, logical(1)),
  n = vapply(design_info, function(x) x$n, character(1)),
  lfc = vapply(design_info, function(x) x$lfc, character(1))
)

design_dt <- design_dt[valid == TRUE]
if (nrow(design_dt) == 0) {
  png("Fig_ObsFDR_Power_Grid_empty.png", width = 1800, height = 1400)
  plot.new()
  text(0.5, 0.5, "No scenario matched dataset_n{10|20|30}_lfc{1|1p5|2}", cex = 1.3)
  dev.off()
  quit(save = "no")
}

data_all <- merge(data_all, design_dt[, .(scenario, n, lfc)], by = "scenario", all.x = FALSE, all.y = FALSE)

agg <- data_all[, .(
  observed_fdr_mean = mean(observed_fdr, na.rm = TRUE),
  power_mean = mean(power, na.rm = TRUE),
  observed_fdr_sd = if (.N > 1) sd(observed_fdr, na.rm = TRUE) else NA_real_,
  power_sd = if (.N > 1) sd(power, na.rm = TRUE) else NA_real_,
  n_runs = .N
), by = .(scenario, n, lfc, method, threshold)]

agg <- agg[is.finite(observed_fdr_mean) & is.finite(power_mean)]
if (nrow(agg) == 0) {
  png("Fig_ObsFDR_Power_Grid_empty.png", width = 1800, height = 1400)
  plot.new()
  text(0.5, 0.5, "No finite points after aggregation", cex = 1.3)
  dev.off()
  quit(save = "no")
}

METHOD_ORDER <- c(
  "ttest", "wilcox", "edger", "limma", "deseq2", "nebula", "memento", "BOOTSTRAP3", "BOOTSTRAP2"
)
METHOD_TO_LEGEND <- c(
  "nebula" = "NEBULA", "ttest" = "T-test", "wilcox" = "Wilcoxon",
  "edger" = "edgeR", "limma" = "limma", "deseq2" = "DESeq2", "memento" = "memento", "BOOTSTRAP3" = "BOOTSTRAP3", "BOOTSTRAP2" = "BOOTSTRAP2"
)

available_methods <- intersect(METHOD_ORDER, unique(agg$method))
agg <- agg[method %in% available_methods]
agg$method <- factor(agg$method, levels = available_methods)

if (nrow(agg) == 0) {
  png("Fig_ObsFDR_Power_Grid_empty.png", width = 1800, height = 1400)
  plot.new()
  text(0.5, 0.5, "No methods matched plotting order", cex = 1.3)
  dev.off()
  quit(save = "no")
}

groups <- list(
  "Wilcoxon" = c("ttest", "wilcox"),
  "Pseudobulk" = c("edger", "limma", "deseq2"),
  "Mixed" = c("nebula"),
  "Deconvolution" = c("memento"),
  "Bootstrapping" = c("BOOTSTRAP3", "BOOTSTRAP2")
)
base_colors <- c(
  "Wilcoxon" = "#000000ff",
  "Pseudobulk" = "#ff0000ff",
  "Mixed" = "#fdaa05ff",
  "Deconvolution" = "#0080ffff",
  "Bootstrapping" = "#00aa55ff"
)

method_colors <- c()
for (g in names(groups)) {
  methods_in_group <- intersect(groups[[g]], available_methods)
  if (length(methods_in_group) > 0) {
    base <- base_colors[g]
    n <- length(methods_in_group)
    if (n == 1) {
      cols <- base
    } else {
      c_light <- colorRampPalette(c("white", base))(10)[4]
      c_dark <- colorRampPalette(c(base, "black"))(10)[4]
      cols <- colorRampPalette(c(c_light, c_dark))(n)
    }
    names(cols) <- methods_in_group
    method_colors <- c(method_colors, cols)
  }
}

METHOD_SHAPES <- c(
  "ttest" = 16, "wilcox" = 17, "edger" = 15, "limma" = 18,
  "deseq2" = 3, "nebula" = 8, "memento" = 1,
  "BOOTSTRAP3" = 2, "BOOTSTRAP2" = 0
)
METHOD_LINETYPES <- c(
  "ttest" = "solid", "wilcox" = "dashed", "edger" = "dotted", "limma" = "dotdash",
  "deseq2" = "longdash", "nebula" = "solid", "memento" = "dashed",
  "BOOTSTRAP3" = "dotted", "BOOTSTRAP2" = "dotdash"
)

agg$n <- factor(agg$n, levels = c("10", "20", "30"))
agg$lfc <- factor(agg$lfc, levels = c("1", "1.5", "2"))
agg <- agg[order(method, scenario, threshold)]

max_x <- max(agg$observed_fdr_mean + ifelse(is.na(agg$observed_fdr_sd), 0, agg$observed_fdr_sd), na.rm = TRUE)
if (!is.finite(max_x)) {
  max_x <- 0.2
}
x_upper <- min(1, max(0.1, max_x * 1.05))

p <- ggplot(agg, aes(x = observed_fdr_mean, y = power_mean, color = method, linetype = method, shape = method)) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "red", alpha = 0.7, linewidth = 0.8) +
  geom_path(aes(group = interaction(method, scenario)), linewidth = 0.7, alpha = 0.72) +
  geom_point(size = 1.4, alpha = 0.82) +
  facet_grid(rows = vars(lfc), cols = vars(n), labeller = labeller(
    n = function(x) paste0("n=", x),
    lfc = function(x) paste0("LFC=", x)
  )) +
  scale_color_manual(values = method_colors, labels = METHOD_TO_LEGEND[available_methods]) +
  scale_linetype_manual(values = METHOD_LINETYPES[available_methods], labels = METHOD_TO_LEGEND[available_methods]) +
  scale_shape_manual(values = METHOD_SHAPES[available_methods], labels = METHOD_TO_LEGEND[available_methods]) +
  coord_cartesian(xlim = c(0, 0.3), ylim = c(0, 0.8)) +
  labs(
    x = "Observed FDR",
    y = "Power",
    color = "Method",
    linetype = "Method",
    shape = "Method"
  ) +
  theme_cowplot(14) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey15", color = "grey30"),
    strip.text = element_text(color = "white", face = "bold"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.25)
  )

if (show_error_bars) {
  p <- p +
    geom_errorbar(aes(ymin = pmax(0, power_mean - power_sd), ymax = pmin(1, power_mean + power_sd)),
                  linewidth = 0.3, alpha = 0.35) +
    geom_errorbarh(aes(xmin = pmax(0, observed_fdr_mean - observed_fdr_sd), xmax = pmin(1, observed_fdr_mean + observed_fdr_sd)),
                   linewidth = 0.3, alpha = 0.35)
}

output_filename <- "Fig_ObsFDR_Power_Grid_mean_runs.png"
ggsave(output_filename, p, width = 2400, height = 1800, units = "px", dpi = 300, bg = "white")

for (method_id in available_methods) {
  method_label <- METHOD_TO_LEGEND[[method_id]]
  method_dt <- agg[method == method_id]
  p_method <- ggplot(method_dt, aes(x = observed_fdr_mean, y = power_mean)) +
    geom_vline(xintercept = 0.1, linetype = "dashed", color = "red", alpha = 0.7, linewidth = 0.8) +
    geom_path(aes(group = scenario), color = method_colors[[method_id]], linewidth = 0.8, alpha = 0.9) +
    geom_point(color = method_colors[[method_id]], shape = METHOD_SHAPES[[method_id]], size = 1.6, alpha = 0.9) +
    facet_grid(rows = vars(lfc), cols = vars(n), labeller = labeller(
      n = function(x) paste0("n=", x),
      lfc = function(x) paste0("LFC=", x)
    )) +
    coord_cartesian(xlim = c(0, 0.3), ylim = c(0, 0.8)) +
    labs(
      title = method_label,
      x = "Observed FDR",
      y = "Power"
    ) +
    theme_cowplot(14) +
    theme(
      strip.background = element_rect(fill = "grey15", color = "grey30"),
      strip.text = element_text(color = "white", face = "bold"),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey92", linewidth = 0.25)
    )

  method_file_id <- gsub("[^A-Za-z0-9]+", "_", method_id)
  ggsave(
    paste0("Fig_ObsFDR_Power_Grid_mean_runs_", method_file_id, ".png"),
    p_method,
    width = 2000,
    height = 1600,
    units = "px",
    dpi = 300,
    bg = "white"
  )
}
