#install BioconductoR dependencies
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
BiocManager::install("GenomicFeatures")
BiocManager::install("GenomicAlignments")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")

#install other R dependencies
install.packages("data.table")
install.packages("DBI")
install.packages("pool")
install.packages("ggplot2")
install.packages("magrittr")
install.packages("openxlsx")
install.packages("formattable")