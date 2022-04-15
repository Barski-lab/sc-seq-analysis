##########################################################################################################
# Dockerfile
#
# Software:         set of R scripts for single-cell data analysis
# Software Version: v0.0.4
# Description:      Dockerized version of the set of R scripts for
#                   single-cell data analysis
# Website:          https://github.com/Barski-lab/workflows
# Provides:         set of R scripts for single-cell data analysis
# Base Image:       satijalab/seurat:4.0.6
# Build Cmd:        docker build --no-cache --rm -t biowardrobe2/sc-tools:v0.0.4 -f sc-tools-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/sc-tools:v0.0.4
# Run Cmd:          docker run --rm -ti biowardrobe2/sc-tools:v0.0.4 /bin/bash
##########################################################################################################
#
# v0.0.4
# - add/update scripts
#   sc_rna_filter.R
#   sc_multiome_filter.R
#   sc_rna_reduce.R
#   sc_atac_reduce.R
#   sc_rna_cluster.R
#   sc_atac_cluster.R
#   sc_wnn_cluster.R
# - install MACS2 after numba as it fails on numpy
#
# v0.0.3
# - add sc_gex_reduce.R script
# - add sc_gex_cluster.R
# Notes:
# - pip3 install -U numba is required. This will reinstall numpy
#   to the version that is needed by numba thus making possible to
#   use umap in R through reticulate (use umap-learn in runUMAP).
#
# v0.0.2
# - add sc_gex_filter.R script
#
# v0.0.1
# - initial version of sc_multiome_filter.R script
#
##########################################################################################################


### Base Image
FROM satijalab/seurat:4.0.6
LABEL maintainer="misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive

WORKDIR /tmp

ENV R_MAX_VSIZE=200000000000
ENV CB_VERSION "1.1.1"
ENV MACS2_VERSION "2.2.7.1"

COPY ./scripts/sc_tools/sc_rna_filter.R /usr/local/bin/sc_rna_filter.R
COPY ./scripts/sc_tools/sc_multiome_filter.R /usr/local/bin/sc_multiome_filter.R
COPY ./scripts/sc_tools/sc_rna_reduce.R /usr/local/bin/sc_rna_reduce.R
COPY ./scripts/sc_tools/sc_atac_reduce.R /usr/local/bin/sc_atac_reduce.R
COPY ./scripts/sc_tools/sc_rna_cluster.R /usr/local/bin/sc_rna_cluster.R
COPY ./scripts/sc_tools/sc_atac_cluster.R /usr/local/bin/sc_atac_cluster.R
COPY ./scripts/sc_tools/sc_wnn_cluster.R /usr/local/bin/sc_wnn_cluster.R
COPY ./scripts/sc_tools/modules/*.R /usr/local/bin/modules/

### Installing dependencies
RUN apt-get update && \
    apt-get install libgcc-10-dev python3-dev python3-pip libxml2-dev libcurl4-openssl-dev && \
                    libssl-dev pandoc libudunits2-dev libgdal-dev libcairo2-dev libharfbuzz-dev && \
                    libfribidi-dev libbz2-dev -y && \
    pip3 install scipy numpy && \
    pip3 install cellbrowser==${CB_VERSION} && \
    pip3 install -U numba && \
    pip3 install MACS2==${MACS2_VERSION} && \
    R -e 'install.packages("argparse", repo = "https://cloud.r-project.org/")' && \
    R -e 'install.packages("tidyverse", repo = "https://cloud.r-project.org/")' && \
    R -e 'install.packages("modules", repo = "https://cloud.r-project.org/")' && \
    R -e "BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats','limma', 'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment', 'batchelor', 'Matrix.utils', 'DESeq2', 'Rsamtools', 'rtracklayer', 'glmGamPoi'))" && \
    R -e 'install.packages("Signac", repo = "https://cloud.r-project.org/")' && \
    R -e 'install.packages("bestNormalize", repo = "https://cloud.r-project.org/")' && \
### Installing R scripts
    chmod +x /usr/local/bin/sc_rna_filter.R && \
    chmod +x /usr/local/bin/sc_multiome_filter.R && \
    chmod +x /usr/local/bin/sc_rna_reduce.R && \
    chmod +x /usr/local/bin/sc_atac_reduce.R && \
    chmod +x /usr/local/bin/sc_rna_cluster.R && \
    chmod +x /usr/local/bin/sc_atac_cluster.R && \
    chmod +x /usr/local/bin/sc_wnn_cluster.R && \
### Cleaning
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/* && \
    strip /usr/local/bin/*; true
