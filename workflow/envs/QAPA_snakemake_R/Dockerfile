# Dockerfile to in
##### BASE IMAGE #####
FROM rocker/tidyverse:4.2.2

##### METADATA #####


# Install required r packages
RUN Rscript -e 'install.packages(c("optparse", "BiocManager", "glue"), dependencies=TRUE, repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install(c("tximport", "satuRn", "limma", "DESeq2", "BioCParallel", "SummarizedExperiment"), dependencies=TRUE)'
