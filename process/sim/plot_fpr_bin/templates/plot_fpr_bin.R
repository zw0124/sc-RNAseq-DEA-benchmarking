#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparser)
  library(ggplot2)
  library(dplyr)
})

p <- arg_parser("Plot FPR boxplots faceted by negative-control sample size")
p <- add_argument(p, "--input-dir", help = "Directory containing fpr_*.tsv files", type = "character", default = ".")
p <- add_argument(p, "--p-val", help = "P-value threshold used to read FPR", type = "numeric", default = 0.05)
p <- add_argument(p, "--output", help = "Output PNG file", type = "character")
argv <- parse_args(p)

input_dir <- argv$input_dir
p_val <- argv$p_val
output_file <- argv$output

if (is.na(p_val)) {
  p_val <- 0.05
}

files <- list.files(input_dir, pattern = "^fpr_.*[.]tsv", full.names = TRUE)
if (length(files) == 0) {
  stop("No FPR files found in input directory.")
}

read_one <- function(f) {
  base <- basename(f)
  no_ext <- if (endsWith(base, ".tsv")) substr(base, 1, nchar(base) - 4) else base
  parts <- strsplit(no_ext, "_")[[1]]

  if (length(parts) < 4 || parts[1] != "fpr") {
    return(NULL)
  }

  method <- parts[length(parts)]
  run <- parts[length(parts) - 1]
  scenario <- paste(parts[2:(length(parts) - 2)], collapse = "_")
  sample_size <- strsplit(tail(strsplit(scenario, "_")[[1]], 1), "v", fixed = TRUE)[[1]][1]

  df <- read.table(f, header = TRUE, sep = "\t")
  if (!("raw_p_value" %in% names(df)) || !("fpr" %in% names(df))) {
    return(NULL)
  }

  df <- df[order(df[["raw_p_value"]]), ]
  idx <- which.min(abs(df[["raw_p_value"]] - p_val))

  data.frame(
    scenario = scenario,
    run = run,
    sample_size = sample_size,
    method = method,
    p_threshold_used = df[["raw_p_value"]][idx],
    fpr = df[["fpr"]][idx],
    stringsAsFactors = FALSE
  )
}

rows <- lapply(files, read_one)
rows <- rows[!vapply(rows, is.null, logical(1))]
if (length(rows) == 0) {
  stop("No valid FPR tables found.")
}

plot_df <- bind_rows(rows) %>%
  mutate(
    sample_size = ifelse(sample_size %in% c("5", "10", "15"), sample_size, scenario),
    sample_size = factor(sample_size, levels = c("5", "10", "15"), labels = c("n=10", "n=20", "n=30"))
  )

METHOD_ORDER <- c(
  "ttest", "wilcox", "edger", "limma", "deseq2", "nebula",
  "memento", "BOOTSTRAP3", "BOOTSTRAP2"
)
available_methods <- intersect(METHOD_ORDER, unique(plot_df[["method"]]))
if (length(available_methods) == 0) {
  available_methods <- sort(unique(plot_df[["method"]]))
}
plot_df <- plot_df %>%
  mutate(method = factor(method, levels = available_methods))

y_upper <- max(plot_df[["fpr"]], na.rm = TRUE)
if (!is.finite(y_upper)) {
  y_upper <- 0.1
}
y_upper <- min(1, max(0.1, y_upper * 1.05))

p <- ggplot(plot_df, aes(x = method, y = fpr, fill = method)) +
  geom_boxplot(width = 0.62, alpha = 0.85, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.6) +
  geom_hline(yintercept = p_val, linetype = "dashed", color = "firebrick", linewidth = 0.7) +
  facet_wrap(vars(sample_size), nrow = 1) +
  labs(
    x = "Method",
    y = "FPR"
  ) +
  scale_y_continuous(breaks = function(x) pretty(x, n = 6)) +
  coord_cartesian(ylim = c(0, y_upper)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.background = element_rect(fill = "grey15", color = "grey30"),
    strip.text = element_text(color = "white", face = "bold"),
    panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.6),
    panel.spacing.x = grid::unit(0.8, "lines")
  )

ggsave(output_file, p, width = 13, height = 6, dpi = 300)
