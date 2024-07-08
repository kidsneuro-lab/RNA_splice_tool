FROM r-base:latest

# UNIX
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git

# R
RUN R -e "install.packages(c('BiocManager','data.table', 'DBI', 'pool', 'ggplot2', 'magrittr', 'RSQLite', 'openxlsx', 'formattable', 'devtools'), repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('GenomicRanges','GenomicFeatures','GenomicAlignments','BSgenome.Hsapiens.UCSC.hg38','BSgenome.Hsapiens.UCSC.hg19','BSgenome.Hsapiens.1000genomes.hs37d5'))"
