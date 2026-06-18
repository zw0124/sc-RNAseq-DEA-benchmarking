#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(argparser)
})

p <- arg_parser("Run Reactome pathway enrichment analysis on filtered genes")
p <- add_argument(p, "--method",         help = "DE method name",                     type = "character")
p <- add_argument(p, "--filtered-genes", help = "Path to filtered genes TSV",         type = "character")
p <- add_argument(p, "--universe-file",  help = "Path to universe genes TSV",         type = "character")
p <- add_argument(p, "--output-file",    help = "Output enrichment result TSV",       type = "character")
argv <- parse_args(p)

method              <- argv$method
filtered_genes_file <- argv$filtered_genes
universe_file       <- argv$universe_file
output_file         <- argv$output_file

empty_out <- data.table(method = character(0), Description = character(0),
                        p.adjust = numeric(0), GeneRatio = character(0), Count = integer(0))

write_empty <- function() {
  fwrite(empty_out, output_file, sep = "\t")
  quit(save = "no")
}

read_gene_column <- function(file) {
  dt <- fread(file)
  if (nrow(dt) == 0) return(character(0))
  if ("gene" %in% names(dt)) {
    genes <- dt$gene
  } else {
    genes <- dt[[1]]
  }
  genes <- unique(as.character(genes))
  genes[!is.na(genes) & genes != ""]
}

map_ensembl_to_entrez <- function(genes) {
  genes <- unique(genes[!is.na(genes) & genes != ""])
  if (length(genes) == 0) return(character(0))

  mapped <- tryCatch(
    suppressMessages(
      AnnotationDbi::select(
        org.Hs.eg.db,
        keys = genes,
        keytype = "ENSEMBL",
        columns = c("ENTREZID")
      )
    ),
    error = function(e) NULL
  )

  if (is.null(mapped) || nrow(mapped) == 0 || !("ENTREZID" %in% names(mapped))) {
    return(character(0))
  }

  entrez <- unique(as.character(mapped$ENTREZID))
  entrez[!is.na(entrez) & entrez != ""]
}

sig_genes <- read_gene_column(filtered_genes_file)
if (length(sig_genes) == 0) write_empty()

universe_genes <- read_gene_column(universe_file)
if (length(universe_genes) == 0) write_empty()

sig_entrez <- map_ensembl_to_entrez(sig_genes)
universe_entrez <- map_ensembl_to_entrez(universe_genes)
sig_entrez <- intersect(sig_entrez, universe_entrez)

if (length(sig_entrez) == 0 || length(universe_entrez) == 0) write_empty()

enrich_res <- tryCatch(
  ReactomePA::enrichPathway(
    gene = sig_entrez,
    universe = universe_entrez,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  ),
  error = function(e) NULL
)

if (is.null(enrich_res)) write_empty()

res <- as.data.frame(enrich_res)
if (nrow(res) == 0 || !("p.adjust" %in% names(res))) write_empty()

res_dt <- as.data.table(res)
res_dt[, method := method]
out_dt <- res_dt[, .(method, Description, p.adjust, GeneRatio, Count)]
fwrite(out_dt, output_file, sep = "\t")
