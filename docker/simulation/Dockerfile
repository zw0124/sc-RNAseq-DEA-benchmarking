FROM bioconductor/bioconductor_docker:RELEASE_3_23

RUN R -e "install.packages('argparser')"
RUN R -e "BiocManager::install(c('splatter','SingleCellExperiment','anndataR','rhdf5','VariantAnnotation'), update=FALSE, ask=FALSE)"
RUN git clone https://github.com/bmbolstad/preprocessCore.git \
    && R CMD INSTALL --configure-args='--disable-threading' preprocessCore \
    && rm -rf preprocessCore
