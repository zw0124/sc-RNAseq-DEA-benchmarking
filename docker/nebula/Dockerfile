FROM bioconductor/bioconductor_docker:RELEASE_3_16

RUN R -e "install.packages(c('argparser','remotes','fs'))"
RUN R -e "BiocManager::install('SingleCellExperiment', update=FALSE, ask=FALSE)"
RUN R -e "remotes::install_github('lhe17/nebula')"
