#!/usr/bin/env bash
# install_all_packages.sh
# ---------------------------------------
# Installs CRAN and Bioconductor dependencies for your R package:
# Depends: data.table, ggplot2
# Imports: data.table, GenomicFeatures, GenomicAlignments, GenomicRanges,
#          BSgenome.Hsapiens.UCSC.hg38, BSgenome.Hsapiens.UCSC.hg19,
#          BSgenome.Hsapiens.1000genomes.hs37d5, magrittr, openxlsx, formattable

set -euo pipefail

# 1. Ensure Rscript is available
if ! command -v Rscript &>/dev/null; then
  echo "ERROR: Rscript not found. Please install R first." >&2
  exit 1
fi

# 2. Run a single Rscript block to install CRAN + Bioconductor packages
Rscript -e "
# --- CRAN packages ---
cran_pkgs <- c(
  'data.table',
  'DBI',
  'pool',
  'ggplot2',
  'magrittr',
  'openxlsx',
  'formattable',
  'optparse'
)

# check and install missing CRAN pkgs
installed_cran <- rownames(installed.packages())
to_install_cran <- setdiff(cran_pkgs, installed_cran)
if (length(to_install_cran)) {
  message('Installing CRAN packages: ', paste(to_install_cran, collapse = ', '))
  install.packages(to_install_cran)
} else {
  message('All CRAN packages are already installed.')
}

# --- Bioconductor packages ---
if (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager')
}
bioc_pkgs <- c(
  'GenomicFeatures',
  'GenomicAlignments',
  'GenomicRanges',
  'BSgenome.Hsapiens.UCSC.hg38',
  'BSgenome.Hsapiens.UCSC.hg19',
  'BSgenome.Hsapiens.1000genomes.hs37d5'
)

installed_all <- rownames(installed.packages())
to_install_bioc <- setdiff(bioc_pkgs, installed_all)
if (length(to_install_bioc)) {
  message('Installing Bioconductor packages: ', paste(to_install_bioc, collapse = ', '))
  BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)
} else {
  message('All Bioconductor packages are already installed.')
}
"

echo "✅ All requested packages are installed."