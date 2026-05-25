FROM bioconductor/bioconductor_docker:RELEASE_3_23

RUN R -e "install.packages(c('argparser','Seurat'))"
RUN R -e "BiocManager::install('DESeq2', update=FALSE, ask=FALSE)"
