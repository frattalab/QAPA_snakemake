# Dockerfile to in
##### BASE IMAGE #####
FROM quay.io/biocontainers/bioconductor-dexseq:1.44.0--r42hdfd78af_0

##### METADATA #####


# Install required r packages
RUN Rscript -e 'install.packages(c("tidyverse", "optparse", "BiocManager", "glue"), dependencies=TRUE, repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install("tximport", dependencies=TRUE, version = "3.16")'
