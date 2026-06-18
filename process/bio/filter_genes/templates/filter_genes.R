#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(argparser)
})

p <- arg_parser("Filter DE genes by LFC and adjusted p-value thresholds")
p <- add_argument(p, "--method",        help = "DE method name",              type = "character")
p <- add_argument(p, "--de-result",     help = "Path to DE result TSV",       type = "character")
p <- add_argument(p, "--lfc-threshold", help = "Log fold change threshold",   type = "numeric", default = 0)
p <- add_argument(p, "--adj-p-cutoff",  help = "Adjusted p-value cutoff",     type = "numeric", default = 0.05)
p <- add_argument(p, "--output-file",   help = "Output filtered genes TSV",   type = "character")
argv <- parse_args(p)

method        <- argv$method
de_result     <- argv$de_result
lfc_threshold <- argv$lfc_threshold
adj_p_cutoff  <- argv$adj_p_cutoff
output_file   <- argv$output_file

pick_first <- function(cols, candidates) {
  hit <- intersect(candidates, cols)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

dt <- fread(de_result)

if (nrow(dt) > 0) {
  gene_col <- pick_first(names(dt), c("gene", "Gene", "genes", "symbol", "SYMBOL"))
  if (is.na(gene_col)) {
    gene_col <- names(dt)[1]
  }
  setnames(dt, gene_col, "gene")

  adj_col <- pick_first(names(dt), c("p_val_adj", "padj", "adj.P.Val", "adj_p", "FDR", "q_val", "qvalue", "p.adjust"))
  raw_col <- pick_first(names(dt), c("p_val", "pvalue", "P.Value", "p", "de_pval"))
  lfc_col <- pick_first(names(dt), c("log2FC", "lfc", "logFC", "avg_log2FC", "avg_logFC", "de_coef"))

  if (is.na(adj_col) && !is.na(raw_col)) {
    dt[, p_adj_bh := p.adjust(get(raw_col), method = "BH")]
    adj_col <- "p_adj_bh"
  }
  
  if (!is.na(adj_col) && !is.na(lfc_col)) {
    dt[, gene := as.character(gene)]
    dt <- dt[!is.na(gene) & gene != ""]
    dt[, adj_value := suppressWarnings(as.numeric(get(adj_col)))]
    dt[, lfc_value := suppressWarnings(as.numeric(get(lfc_col)))]

    sig <- dt[!is.na(adj_value) & !is.na(lfc_value) & adj_value <= adj_p_cutoff & abs(lfc_value) >= lfc_threshold]
    
    if (nrow(sig) > 0) {
      sig_genes <- unique(sig$gene)
      out_dt <- data.table(method = method, gene = sig_genes)
      fwrite(out_dt, output_file, sep = "\t")
      quit(save = "no")
    }
  }
}

fwrite(data.table(method = character(0), gene = character(0)), output_file, sep = "\t")
