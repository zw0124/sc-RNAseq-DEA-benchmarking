#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(purrr)

files <- list.files(pattern = "^fpr_.*[.]tsv")
if (length(files) == 0) {
  stop("No FPR files found in work directory.")
}

message(sprintf("[PLOT_FPR] files found: %d", length(files)))

plot_title <- "$plot_title"
output_file <- "$output_file"

read_data <- function(f) {
  
  base <- basename(f)
  base_no_ext <- if (endsWith(base, ".tsv")) substr(base, 1, nchar(base) - 4) else base
  parts <- strsplit(base_no_ext, "_")[[1]]
  method <- parts[length(parts)]
  run <- parts[length(parts) - 1]
  scenario <- paste(parts[2:(length(parts) - 2)], collapse = "_")
  sample_size <- strsplit(tail(strsplit(scenario, "_")[[1]], 1), "v", fixed = TRUE)[[1]][1]
  
  df <- read.table(f, header = TRUE, sep = "\t")
  df[["method"]] <- method
  df[["run"]] <- run
  df[["scenario"]] <- scenario
  df[["sample_size"]] <- sample_size
  df[["file_id"]] <- f
  return(df)
}

data_list <- lapply(files, read_data)
data <- bind_rows(data_list)

message(sprintf("[PLOT_FPR] rows loaded: %d", nrow(data)))
message(sprintf("[PLOT_FPR] methods detected: %s", paste(sort(unique(data[["method"]])), collapse = ", ")))

if (!all(c("raw_p_value", "fpr", "method", "file_id") %in% names(data))) {
  stop("FPR tables are missing required columns: raw_p_value/fpr.")
}

data <- data %>%
  mutate(
    sample_size = ifelse(sample_size %in% c("5", "10", "15"), sample_size, scenario),
    sample_size = factor(sample_size, levels = c("5", "10", "15"), labels = c("n=10", "n=20", "n=30"))
  )

x_grid <- seq(0.01, 0.10, by = 0.01)

interp_func <- function(df) {
  df <- df[is.finite(df[["raw_p_value"]]) & is.finite(df[["fpr"]]), ]
  df <- df[order(df[["raw_p_value"]]), ]
  if (nrow(df) < 2 || length(unique(df[["raw_p_value"]])) < 2) {
    message(sprintf("[PLOT_FPR] skipping %s: fewer than two finite raw_p_value points", df[["file_id"]][1]))
    return(NULL)
  }
  res <- approx(df[["raw_p_value"]], df[["fpr"]], xout = x_grid, rule = 2)
  data.frame(
    raw_p_value = res[["x"]],
    fpr = res[["y"]],
    method = df[["method"]][1],
    sample_size = df[["sample_size"]][1],
    file_id = df[["file_id"]][1]
  )
}

interpolated_list <- data %>%
  split(.[["file_id"]]) %>%
  map(interp_func)

interpolated <- bind_rows(interpolated_list)
if (nrow(interpolated) == 0) {
  stop("No FPR tables had at least two finite points for interpolation.")
}

summary_df <- interpolated %>%
  group_by(sample_size, method, raw_p_value) %>%
  summarise(
    mean_fpr = mean(fpr, na.rm = TRUE),
    sd_fpr = sd(fpr, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    se_fpr = sd_fpr / sqrt(n),
    lower = pmax(0, mean_fpr - se_fpr),
    upper = pmin(1, mean_fpr + se_fpr)
  )

summary_file <- paste0(tools::file_path_sans_ext(output_file), ".summary.tsv")
write.table(summary_df, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
message(sprintf("[PLOT_FPR] summary saved: %s", summary_file))

p <- ggplot(summary_df, aes(x = raw_p_value, y = mean_fpr, color = method, fill = method)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  facet_wrap(vars(sample_size), nrow = 1) +
  labs(
    x = "p-value threshold",
    y = "False Positive Rate (FPR)",
    title = plot_title
  ) +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 0.10), ylim = c(0, 0.10)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray")

ggsave(output_file, p, width = 13, height = 6, dpi = 300)
