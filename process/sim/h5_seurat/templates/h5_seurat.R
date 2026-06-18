#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(anndataR)
  library(SingleCellExperiment)
  library(Seurat)
  library(Matrix)
})

sce <- read_h5ad("${input_anndata}", as = "SingleCellExperiment", x_mapping = "counts")

cell_meta <- as.data.frame(colData(sce))
stopifnot("counts" %in% assayNames(sce))
stopifnot("Condition" %in% colnames(cell_meta))
stopifnot("Sample" %in% colnames(cell_meta))

mainExpName(sce) <- "RNA"
assay(sce, "counts") <- as(assay(sce, "counts"), "dgCMatrix")

seu <- as.Seurat(sce, counts = "counts", data = "counts")
DefaultAssay(seu) <- "RNA"

rna_counts <- Seurat::GetAssayData(seu, assay = "RNA", layer = "counts")
seu[["nCount_RNA"]] <- Matrix::colSums(rna_counts)
seu[["nFeature_RNA"]] <- Matrix::colSums(rna_counts > 0)

saveRDS(seu, file = "${meta.scenario}_${meta.run}.seurat.rds")

colnames(sce) <- as.character(colnames(sce))

write_h5ad(
  sce,
  "${meta.scenario}_${meta.run}.h5_seurat.h5ad",
  x_mapping = "counts"
)
