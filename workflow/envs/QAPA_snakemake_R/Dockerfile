# Dockerfile to in
##### BASE IMAGE #####
FROM rocker/tidyverse:4.2.2

##### METADATA #####


# Install required r packages
RUN Rscript -e 'install.packages(c("optparse", "BiocManager", "glue", "devtools"), dependencies=TRUE, repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'devtools::install_github("SamBryce-Smith/satuRn", ref = "9c346e8")'
RUN Rscript -e 'BiocManager::install(c("tximport", "limma", "DESeq2", "BioCParallel", "SummarizedExperiment"), dependencies=TRUE)'
