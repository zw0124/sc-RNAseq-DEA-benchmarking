#!/usr/bin/env Rscript

suppressPackageStartupMessages({
     library(data.table)
     library(ggplot2)
     library(dplyr)
     library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
     pos <- match(flag, args)
     if (is.na(pos) || pos == length(args)) {
          stop(sprintf("Missing argument: %s", flag))
     }
     args[pos + 1]
}

fixed_cells <- get_arg("--fixed-cells")
fixed_genes <- get_arg("--fixed-genes")
output_png <- get_arg("--output")

results.cells.fixed <- fread(fixed_cells)
results.genes.fixed <- fread(fixed_genes)

metric_cols <- setdiff(colnames(results.cells.fixed), c("n_genes", "n_cells"))
metric_cols <- metric_cols[metric_cols %in% colnames(results.genes.fixed)]

method_order <- c(
     "ttest",
     "wilcox",
     "edger",
     "limma",
     "deseq2",
     "nebula",
     "memento",
     "bootstrap3",
     "bootstrap2"
)

method_to_legend <- c(
     ttest = "T-test",
     wilcox = "Wilcoxon",
     edger = "edgeR",
     limma = "limma",
     deseq2 = "DESeq2",
     nebula = "NEBULA",
     memento = "memento",
     bootstrap3 = "BOOTSTRAP3",
     bootstrap2 = "BOOTSTRAP2"
)

metric_cols <- intersect(method_order, metric_cols)

groups <- list(
     Wilcoxon = c("ttest", "wilcox"),
     Pseudobulk = c("edger", "limma", "deseq2"),
     Mixed = c("nebula"),
     Deconvolution = c("memento"),
     Bootstrapping = c("bootstrap3", "bootstrap2")
)

base_colors <- c(
     Wilcoxon = "#000000ff",
     Pseudobulk = "#ff0000ff",
     Mixed = "#fdaa05ff",
     Deconvolution = "#0080ffff",
     Bootstrapping = "#00aa55ff"
)

method_colors <- c()
for (g in names(groups)) {
     methods_in_group <- intersect(groups[[g]], metric_cols)
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

results.cells.long <- results.cells.fixed %>%
     melt(
          id.vars = c("n_genes", "n_cells"),
          measure.vars = metric_cols,
          value.name = "time",
          variable.name = "method"
     ) %>%
     filter(!is.na(time))

results.genes.long <- results.genes.fixed %>%
     melt(
          id.vars = c("n_genes", "n_cells"),
          measure.vars = metric_cols,
          value.name = "time",
          variable.name = "method"
     ) %>%
     filter(!is.na(time))

results.cells.long[["method"]] <- factor(results.cells.long[["method"]], levels = metric_cols)
results.genes.long[["method"]] <- factor(results.genes.long[["method"]], levels = metric_cols)

p1 <- ggplot(results.cells.long, aes(x = n_genes, y = time, color = method)) +
     geom_line(linewidth = 0.8) +
     labs(x = "Genes", y = "Time (s)", color = "Method") +
     scale_color_manual(values = method_colors, breaks = metric_cols, labels = method_to_legend[metric_cols], drop = FALSE) +
     theme_cowplot(14) +
     theme(
          legend.title = element_text(face = "bold"),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          legend.background = element_rect(fill = "white", color = NA),
          legend.key = element_rect(fill = "white", color = NA)
     ) +
     background_grid(major = "xy", minor = "none", size.major = 0.2, size.minor = 0.2)

p2 <- ggplot(results.genes.long, aes(x = n_cells, y = time, color = method)) +
     geom_line(linewidth = 0.8) +
     labs(x = "Cells", y = "Time (s)", color = "Method") +
     scale_color_manual(values = method_colors, breaks = metric_cols, labels = method_to_legend[metric_cols], drop = FALSE) +
     theme_cowplot(14) +
     theme(
          legend.title = element_text(face = "bold"),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          legend.background = element_rect(fill = "white", color = NA),
          legend.key = element_rect(fill = "white", color = NA)
     ) +
     background_grid(major = "xy", minor = "none", size.major = 0.2, size.minor = 0.2)

legend.double.plot <- get_legend(p1 + theme(legend.margin = margin(0, 0, 20, 20)))

double.plot <- plot_grid(
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

runtime_benchmark <- plot_grid(double.plot, legend.double.plot, ncol = 2, rel_widths = c(2.8, 0.8))

ggsave(runtime_benchmark, width = 3194, height = 1600, units = "px", dpi = 300, filename = output_png, bg = "white")
