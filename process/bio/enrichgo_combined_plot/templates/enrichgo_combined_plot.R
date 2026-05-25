#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(argparser)
})

p <- arg_parser("Generate combined bubble plot from enrichment results")
p <- add_argument(p, "--color-cap", help = "Upper cap for -log10(q) color scale", type = "numeric", default = 8)
argv <- parse_args(p)

top_n <- 5
y_label <- "Enriched Reactome pathway"
color_cap <- argv$color_cap

bubble_file <- "combined_reactome_bubble.png"
terms_file <- "combined_reactome_terms.tsv"

strategy_levels <- c("lfc0_filter0", "lfc0_filter0p5", "lfc0p5_filter0p5")
strategy_labels <- c(
  "A",
  "B",
  "C"
)
names(strategy_labels) <- strategy_levels

pathway_group_levels <- c("IFN-related pathways", "Other pathways")
pathway_block_levels <- pathway_group_levels

parse_ratio <- function(x) {
  ifelse(
    is.na(x),
    NA_real_,
    vapply(strsplit(as.character(x), "/", fixed = TRUE), function(parts) {
      if (length(parts) != 2) return(NA_real_)
      num <- suppressWarnings(as.numeric(parts[1]))
      den <- suppressWarnings(as.numeric(parts[2]))
      if (is.na(num) || is.na(den) || den == 0) return(NA_real_)
      num / den
    }, numeric(1))
  )
}

parse_ratio_den <- function(x) {
  ifelse(
    is.na(x),
    NA_real_,
    vapply(strsplit(as.character(x), "/", fixed = TRUE), function(parts) {
      if (length(parts) != 2) return(NA_real_)
      den <- suppressWarnings(as.numeric(parts[2]))
      if (is.na(den) || den == 0) return(NA_real_)
      den
    }, numeric(1))
  )
}

classify_pathway <- function(terms) {
  vapply(tolower(as.character(terms)), function(term) {
    if (grepl("interferon|\\bifn\\b|antiviral|virus|viral", term)) {
      "IFN-related pathways"
    } else {
      "Other pathways"
    }
  }, character(1))
}


classify_pathway_block <- function(groups) {
  as.character(groups)
}


read_enrich_file <- function(f) {
  base <- basename(f)
  strategy <- NA_character_
  for (s in strategy_levels) {
    if (grepl(paste0("^bio_enrich_reactome_", s, "_"), base)) {
      strategy <- s
      break
    }
  }
  if (is.na(strategy)) return(NULL)

  dt <- fread(f)
  if (nrow(dt) == 0) return(NULL)

  req_cols <- c("method", "Description", "p.adjust", "GeneRatio", "Count")
  if (!all(req_cols %in% names(dt))) return(NULL)

  dt <- dt[, ..req_cols]
  dt[, strategy := strategy]
  dt[]
}

tsv_files <- list.files(".", pattern = "^bio_enrich_reactome_.*\\.tsv$", full.names = TRUE)
enrich_rows <- lapply(tsv_files, read_enrich_file)
enrich_rows <- enrich_rows[!vapply(enrich_rows, is.null, logical(1))]

all_enrich_dt <- if (length(enrich_rows) > 0) {
  rbindlist(enrich_rows, fill = TRUE)
} else {
  data.table(strategy = character(0), method = character(0), Description = character(0),
             p.adjust = numeric(0), GeneRatio = character(0), Count = integer(0))
}

if (nrow(all_enrich_dt) == 0) {
  fwrite(all_enrich_dt, terms_file, sep = "\t")
  p <- ggplot() +
    theme_void()
  ggsave(bubble_file, p, width = 16, height = 8, dpi = 300)
  quit(save = "no")
}

all_enrich_dt <- all_enrich_dt[strategy %in% strategy_levels]
all_enrich_dt[, strategy := factor(strategy, levels = strategy_levels)]

top_rows <- list()
for (s in strategy_levels) {
  for (m in unique(all_enrich_dt[strategy == s]$method)) {
    sub <- all_enrich_dt[strategy == s & method == m]
    setorder(sub, p.adjust)
    top_rows[[length(top_rows) + 1]] <- sub[1:min(top_n, .N)]
  }
}

top_union_dt <- rbindlist(top_rows, fill = TRUE)
term_union <- unique(top_union_dt$Description)

strategy_methods <- unique(all_enrich_dt[, .(strategy, method)])
grid_dt <- strategy_methods[, .(Description = term_union), by = .(strategy, method)]

plot_dt <- merge(
  grid_dt,
  all_enrich_dt,
  by = c("strategy", "method", "Description"),
  all.x = TRUE
)

plot_dt[, GeneRatioNum := parse_ratio(GeneRatio)]
plot_dt[is.na(GeneRatioNum), GeneRatioNum := 0]
plot_dt[, NegLogAdjP := fifelse(!is.na(p.adjust) & p.adjust > 0, -log10(p.adjust), NA_real_)]
plot_dt[, NegLogAdjPColor := pmin(NegLogAdjP, color_cap)]
plot_dt[, DrawPoint := is.na(p.adjust) | GeneRatioNum > 0]
plot_dt[, PathwayGroup := factor(classify_pathway(Description), levels = pathway_group_levels)]
plot_dt[, PathwayBlock := factor(classify_pathway_block(PathwayGroup), levels = pathway_block_levels)]

method_n <- copy(all_enrich_dt)
method_n[, InputGeneN := parse_ratio_den(GeneRatio)]
method_n <- method_n[!is.na(InputGeneN), .(InputGeneN = max(InputGeneN)), by = .(strategy, method)]
method_n[, method_label := paste0(method, "\nn=", InputGeneN)]
plot_dt <- merge(plot_dt, method_n, by = c("strategy", "method"), all.x = TRUE)
plot_dt[is.na(method_label), method_label := method]

term_order <- top_union_dt[, .(
  best_p = min(p.adjust, na.rm = TRUE),
  PathwayGroup = factor(classify_pathway(Description)[1], levels = pathway_group_levels)
), by = Description]
term_order[, PathwayBlock := factor(classify_pathway_block(PathwayGroup), levels = pathway_block_levels)]
setorder(term_order, PathwayBlock, PathwayGroup, best_p)

term_order[, y_top_index := seq_len(.N)]
term_order[, y_pos := nrow(term_order) - y_top_index + 1]
plot_dt <- merge(
  plot_dt,
  term_order[, .(Description, y_pos)],
  by = "Description",
  all.x = TRUE
)
plot_dt[, Description := factor(Description, levels = rev(term_order$Description))]

plot_dt[, strategy_label := strategy_labels[as.character(strategy)]]
plot_dt[, strategy_label := factor(strategy_label, levels = strategy_labels)]

method_label_levels <- unique(plot_dt[order(strategy, method)]$method_label)
strip_x_label <- " "
plot_dt[, method_label := factor(method_label, levels = c(strip_x_label, method_label_levels))]

fwrite(plot_dt, terms_file, sep = "\t")

color_breaks <- pretty(c(0, color_cap), n = 5)
color_breaks <- color_breaks[color_breaks >= 0 & color_breaks <= color_cap]
if (!color_cap %in% color_breaks) {
  color_breaks <- sort(unique(c(color_breaks, color_cap)))
}
color_labels <- as.character(color_breaks)
color_labels[length(color_labels)] <- paste0(">=", color_cap)

block_bands <- term_order[
  ,
  .(
    start_index = min(y_top_index),
    end_index = max(y_top_index)
  ),
  by = PathwayBlock
]
block_bands[, `:=`(
  xmin = -Inf,
  xmax = Inf,
  ymin = nrow(term_order) - end_index + 0.5,
  ymax = nrow(term_order) - start_index + 1.5
)]
block_palette <- c("#e31a8a", "#63d2c6", "#bdbdbd")
block_bands[, BlockFill := block_palette[seq_len(.N)]]
block_label_map <- c(
  "IFN-related pathways" = "IFN-related pathways",
  "Other pathways" = "Other pathways"
)
block_bands[, BlockLabel := block_label_map[as.character(PathwayBlock)]]
block_bands[, ymid := (ymin + ymax) / 2]
block_bands[, method_label := factor(strip_x_label, levels = c(strip_x_label, method_label_levels))]
block_bands[, strategy_label := factor(strategy_labels[1], levels = strategy_labels)]

group_breaks <- copy(block_bands[end_index < nrow(term_order)])
group_breaks[, strategy_label := NULL]
group_breaks[, yintercept := nrow(term_order) - end_index + 0.5]

p <- ggplot(plot_dt, aes(x = method_label, y = y_pos)) +
  geom_tile(
    data = block_bands,
    aes(x = method_label, y = ymid, height = ymax - ymin, fill = BlockFill),
    inherit.aes = FALSE,
    width = 0.28,
    alpha = 0.95
  ) +
  geom_text(
    data = block_bands,
    aes(x = method_label, y = ymid, label = BlockLabel),
    inherit.aes = FALSE,
    angle = 90,
    size = 3.4,
    color = "#1a1a1a"
  ) +
  geom_hline(
    data = group_breaks,
    aes(yintercept = yintercept),
    color = "#9e9e9e",
    linewidth = 0.35
  ) +
  geom_point(
    data = plot_dt[DrawPoint == TRUE],
    aes(size = GeneRatioNum, color = NegLogAdjPColor),
    alpha = 0.9,
    na.rm = FALSE
  ) +
  facet_grid(. ~ strategy_label, scales = "free_x", space = "free_x") +
  scale_size_continuous(name = "GeneRatio", range = c(1.2, 9)) +
  scale_color_gradient(
    name = expression(-log[10]("q")),
    low = "#4575b4",
    high = "#d73027",
    na.value = "#d9d9d9",
    limits = c(0, color_cap),
    breaks = color_breaks,
    labels = color_labels
  ) +
  scale_fill_identity() +
  scale_y_continuous(
    breaks = term_order$y_pos,
    labels = term_order$Description,
    limits = c(0.5, nrow(term_order) + 0.5),
    expand = expansion(mult = 0)
  ) +
  labs(x = "Algorithm", y = y_label) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9),
    panel.grid.major = element_line(color = "#efefef"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#f2f2f2", color = "#bdbdbd"),
    strip.text.x = element_text(face = "bold", size = 11),
    panel.spacing.x = unit(0.8, "lines")
  )

plot_height_px <- min(16000, max(2600, 450 + length(term_union) * 78))
plot_width_px <- 7200

ggsave(
  bubble_file,
  p,
  width = plot_width_px,
  height = plot_height_px,
  units = "px",
  dpi = 300
)
