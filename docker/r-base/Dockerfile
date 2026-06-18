FROM bioconductor/bioconductor_docker:RELEASE_3_23

RUN R -e "install.packages(c('argparser','data.table','dplyr','purrr','ggplot2','cowplot','scales','pracma','gtools','Seurat'))"
RUN R -e "BiocManager::install(c('SingleCellExperiment','anndataR','rhdf5','ReactomePA','org.Hs.eg.db','AnnotationDbi'), update=FALSE, ask=FALSE)"
