#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(argparser)
})

p <- arg_parser("Generate dataset-specific bubble plots from enrichment results")
p <- add_argument(p, "--color-cap", help = "Upper cap for -log10(q) color scale", type = "numeric", default = 8)
argv <- parse_args(p)

top_n <- 5
color_cap <- argv$color_cap
terms_file <- "reactome_terms.tsv"

strategy_levels <- c("lfc0_filter0", "lfc0_filter0p5", "lfc0p5_filter0p5")
strategy_labels <- c("A", "B", "C")
names(strategy_labels) <- strategy_levels

method_order <- c(
  "deseq2",
  "limma",
  "edger",
  "wilcox",
  "ttest",
  "memento",
  "nebula",
  "BOOTSTRAP3",
  "BOOTSTRAP2"
)

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

empty_enrich_dt <- function() {
  data.table(
    dataset = character(0),
    strategy = character(0),
    method = character(0),
    Description = character(0),
    p.adjust = numeric(0),
    GeneRatio = character(0),
    Count = integer(0)
  )
}

parse_enrich_filename <- function(file) {
  base <- basename(file)
  pattern <- paste0(
    "^bio_enrich_reactome_(.+)_(",
    paste(strategy_levels, collapse = "|"),
    ")_(.+)\\.tsv$"
  )
  match <- regmatches(base, regexec(pattern, base))[[1]]
  if (length(match) != 4) return(NULL)
  list(dataset = match[2], strategy = match[3], method = match[4])
}

read_enrich_file <- function(file) {
  info <- parse_enrich_filename(file)
  if (is.null(info)) return(NULL)

  file_method <- data.table(
    dataset = info$dataset,
    strategy = info$strategy,
    method = info$method
  )

  dt <- fread(file)
  req_cols <- c("method", "Description", "p.adjust", "GeneRatio", "Count")
  if (nrow(dt) == 0 || !all(req_cols %in% names(dt))) {
    return(list(methods = file_method, rows = empty_enrich_dt()))
  }

  dt <- dt[, ..req_cols]
  dt[is.na(method) | method == "", method := info$method]
  dt[, dataset := info$dataset]
  dt[, strategy := info$strategy]
  dt <- dt[, .(dataset, strategy, method, Description, p.adjust, GeneRatio, Count)]
  list(methods = file_method, rows = dt)
}

method_sort_key <- function(methods) {
  match_idx <- match(methods, method_order)
  fifelse(is.na(match_idx), length(method_order) + frank(methods, ties.method = "dense"), match_idx)
}

write_empty_plot <- function(dataset) {
  output_file <- sprintf("%s_reactome_bubble.png", dataset)
  p <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No Reactome enrichment results", size = 5) +
    theme_void()
  ggsave(output_file, p, width = 10, height = 6, dpi = 300, bg = "white")
}

build_dataset_plot <- function(dataset, all_methods_dt, all_enrich_dt) {
  target_dataset <- dataset
  dataset_methods <- unique(all_methods_dt[dataset == target_dataset & strategy %in% strategy_levels])
  dataset_enrich <- all_enrich_dt[dataset == target_dataset & strategy %in% strategy_levels]
  output_file <- sprintf("%s_reactome_bubble.png", dataset)

  if (nrow(dataset_methods) == 0 || nrow(dataset_enrich) == 0) {
    write_empty_plot(dataset)
    return(data.table())
  }

  dataset_methods[, strategy := factor(strategy, levels = strategy_levels)]
  dataset_enrich[, strategy := factor(strategy, levels = strategy_levels)]

  top_rows <- list()
  for (s in strategy_levels) {
    methods_s <- dataset_methods[as.character(strategy) == s, method]
    methods_s <- methods_s[order(method_sort_key(methods_s))]
    for (m in methods_s) {
      sub <- dataset_enrich[as.character(strategy) == s & method == m]
      if (nrow(sub) == 0) next
      setorder(sub, p.adjust)
      top_rows[[length(top_rows) + 1]] <- sub[1:min(top_n, .N)]
    }
  }

  if (length(top_rows) == 0) {
    write_empty_plot(dataset)
    return(data.table())
  }

  top_union_dt <- rbindlist(top_rows, fill = TRUE)
  term_union <- unique(top_union_dt$Description)

  grid_dt <- dataset_methods[, .(Description = term_union), by = .(dataset, strategy, method)]
  plot_dt <- merge(
    grid_dt,
    dataset_enrich,
    by = c("dataset", "strategy", "method", "Description"),
    all.x = TRUE
  )

  plot_dt[, GeneRatioNum := parse_ratio(GeneRatio)]
  plot_dt[is.na(GeneRatioNum), GeneRatioNum := 0]
  plot_dt[, NegLogAdjP := fifelse(!is.na(p.adjust) & p.adjust > 0, -log10(p.adjust), NA_real_)]
  plot_dt[, NegLogAdjPColor := pmin(NegLogAdjP, color_cap)]
  plot_dt[, DrawPoint := GeneRatioNum > 0]
  plot_dt[, PathwayGroup := factor(classify_pathway(Description), levels = pathway_group_levels)]
  plot_dt[, PathwayBlock := factor(as.character(PathwayGroup), levels = pathway_block_levels)]

  method_n <- copy(dataset_enrich)
  method_n[, InputGeneN := parse_ratio_den(GeneRatio)]
  method_n <- method_n[!is.na(InputGeneN), .(InputGeneN = max(InputGeneN)), by = .(dataset, strategy, method)]
  method_n[, method_label := paste0(method, "\nn=", InputGeneN)]
  plot_dt <- merge(plot_dt, method_n, by = c("dataset", "strategy", "method"), all.x = TRUE)
  plot_dt[is.na(method_label), method_label := method]

  term_order <- top_union_dt[, .(
    best_p = min(p.adjust, na.rm = TRUE),
    PathwayGroup = factor(classify_pathway(Description)[1], levels = pathway_group_levels)
  ), by = Description]
  term_order[, PathwayBlock := factor(as.character(PathwayGroup), levels = pathway_block_levels)]
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
  plot_dt[, method_sort := method_sort_key(method)]
  method_axis_dt <- unique(plot_dt[, .(strategy, method, method_label, method_sort)])
  setorder(method_axis_dt, strategy, method_sort, method)
  method_axis_dt[, method_axis := paste(as.character(strategy), method, sep = "::")]
  strip_x_key <- "__pathway_block__"
  method_axis_levels <- c(strip_x_key, method_axis_dt$method_axis)
  method_axis_labels <- c(
    setNames("", strip_x_key),
    setNames(method_axis_dt$method_label, method_axis_dt$method_axis)
  )
  plot_dt[, method_axis := factor(
    paste(as.character(strategy), method, sep = "::"),
    levels = method_axis_levels
  )]

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
  block_palette <- c(
    "IFN-related pathways" = "#e31a8a",
    "Other pathways" = "#63d2c6"
  )
  block_bands[, BlockFill := block_palette[as.character(PathwayBlock)]]
  block_label_map <- c(
    "IFN-related pathways" = "IFN-related pathways",
    "Other pathways" = "Other pathways"
  )
  block_bands[, BlockLabel := block_label_map[as.character(PathwayBlock)]]
  block_bands[, ymid := (ymin + ymax) / 2]
  block_bands[, method_axis := factor(strip_x_key, levels = method_axis_levels)]
  block_bands[, strategy_label := factor(strategy_labels[1], levels = strategy_labels)]

  scale_dt <- rbindlist(list(
    plot_dt[, .(method_axis, y_pos, strategy_label)],
    block_bands[, .(method_axis, y_pos = ymid, strategy_label)]
  ))

  group_breaks <- copy(block_bands[end_index < nrow(term_order)])
  group_breaks[, strategy_label := NULL]
  group_breaks[, yintercept := nrow(term_order) - end_index + 0.5]

  p <- ggplot(scale_dt, aes(x = method_axis, y = y_pos)) +
    geom_blank() +
    geom_tile(
      data = block_bands,
      aes(x = method_axis, y = ymid, height = ymax - ymin, fill = BlockFill),
      inherit.aes = FALSE,
      width = 0.28,
      alpha = 0.95
    ) +
    geom_text(
      data = block_bands,
      aes(x = method_axis, y = ymid, label = BlockLabel),
      inherit.aes = FALSE,
      angle = 90,
      size = 4.8,
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
    scale_x_discrete(labels = method_axis_labels) +
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
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 22) +
    theme(
      axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 22),
      axis.text.y = element_text(size = 21),
      panel.grid.major = element_line(color = "#efefef"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey15", color = "grey30"),
      strip.text.x = element_text(
        color = "white",
        face = "bold",
        size = 26,
        margin = margin(t = 10, r = 6, b = 10, l = 6)
      ),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 18),
      panel.spacing.x = unit(1, "lines")
    )

  plot_height_px <- min(16000, max(4000, 900 + length(term_union) * 125))
  plot_width_px <- 10800

  ggsave(
    output_file,
    p,
    width = plot_width_px,
    height = plot_height_px,
    units = "px",
    dpi = 300
  )

  plot_dt[]
}

tsv_files <- list.files(".", pattern = "^bio_enrich_reactome_.*\\.tsv$", full.names = TRUE)
parsed <- lapply(tsv_files, read_enrich_file)
parsed <- parsed[!vapply(parsed, is.null, logical(1))]

all_methods_dt <- if (length(parsed) > 0) {
  unique(rbindlist(lapply(parsed, `[[`, "methods"), fill = TRUE))
} else {
  data.table(dataset = character(0), strategy = character(0), method = character(0))
}

all_enrich_dt <- if (length(parsed) > 0) {
  rbindlist(lapply(parsed, `[[`, "rows"), fill = TRUE)
} else {
  empty_enrich_dt()
}

detected_datasets <- sort(unique(all_methods_dt$dataset))
preferred_datasets <- c("kang", "perez")
dataset_levels <- c(
  intersect(preferred_datasets, detected_datasets),
  setdiff(detected_datasets, preferred_datasets)
)

plot_terms <- list()
for (dataset_name in dataset_levels) {
  plot_terms[[dataset_name]] <- build_dataset_plot(dataset_name, all_methods_dt, all_enrich_dt)
}

terms_dt <- rbindlist(plot_terms, fill = TRUE)
fwrite(terms_dt, terms_file, sep = "\t")
