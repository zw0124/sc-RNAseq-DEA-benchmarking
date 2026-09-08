#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) stop(sprintf("Missing argument: %s", flag))
  args[index + 1]
}

input <- get_arg("--input")
output <- get_arg("--output")
dt <- fread(input)

method_order <- c(
  "ttest", "wilcox", "edger", "limma", "deseq2",
  "nebula", "memento", "bootstrap3", "bootstrap2"
)
method_labels <- c(
  ttest = "T-test", wilcox = "Wilcoxon", edger = "edgeR", limma = "limma",
  deseq2 = "DESeq2", nebula = "NEBULA", memento = "memento",
  bootstrap3 = "BOOTSTRAP3", bootstrap2 = "BOOTSTRAP2"
)
groups <- list(
  Wilcoxon = c("ttest", "wilcox"),
  Pseudobulk = c("edger", "limma", "deseq2"),
  Mixed = "nebula",
  Deconvolution = "memento",
  Bootstrapping = c("bootstrap3", "bootstrap2")
)
base_colors <- c(
  Wilcoxon = "#000000", Pseudobulk = "#ff0000", Mixed = "#fdaa05",
  Deconvolution = "#0080ff", Bootstrapping = "#00aa55"
)

if (nrow(dt) == 0) {
  plot <- ggplot() + annotate("text", x = 0, y = 0, label = "No peak RSS records", size = 5) + theme_void()
  ggsave(output, plot, width = 10, height = 5, dpi = 300, bg = "white")
  quit(save = "no")
}

available <- intersect(method_order, unique(dt$method))
colors <- character(0)
for (group in names(groups)) {
  selected <- intersect(groups[[group]], available)
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

dt <- dt[method %in% available]
dt[, method := factor(method, levels = available)]

fixed_cells <- dt[scenario == "fixed_cells"]
fixed_genes <- dt[scenario == "fixed_genes"]

common_theme <- theme_cowplot(14) +
  theme(
    legend.title = element_text(face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA)
  )

p1 <- ggplot(fixed_cells, aes(x = n_genes, y = mean_peak_rss_gib, color = method)) +
  geom_line(linewidth = 0.8) +
  labs(x = "Genes", y = "Peak RSS (GiB)", color = "Method") +
  scale_color_manual(
    values = colors,
    breaks = available,
    labels = method_labels[available],
    drop = FALSE
  ) +
  common_theme +
  background_grid(major = "xy", minor = "none", size.major = 0.2, size.minor = 0.2)

p2 <- ggplot(fixed_genes, aes(x = n_cells, y = mean_peak_rss_gib, color = method)) +
  geom_line(linewidth = 0.8) +
  labs(x = "Cells", y = "Peak RSS (GiB)", color = "Method") +
  scale_color_manual(
    values = colors,
    breaks = available,
    labels = method_labels[available],
    drop = FALSE
  ) +
  common_theme +
  background_grid(major = "xy", minor = "none", size.major = 0.2, size.minor = 0.2)

legend <- get_legend(p1 + theme(legend.margin = margin(0, 0, 20, 20)))
panels <- plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"),
  labels = "AUTO",
  label_size = 20,
  ncol = 2,
  label_x = 0,
  label_y = 1,
  hjust = -0.1,
  vjust = 1.1,
  align = "h"
)
memory_benchmark <- plot_grid(panels, legend, ncol = 2, rel_widths = c(2.8, 0.8))

ggsave(
  memory_benchmark,
  width = 3194,
  height = 1600,
  units = "px",
  dpi = 300,
  filename = output,
  bg = "white"
)
