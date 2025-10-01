#!/usr/bin/env bash
# cortar.sh — shell wrapper for the R cortar() function using Rscript + optparse

set -euo pipefail

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Error: Rscript not found in PATH." >&2
  exit 127
fi

Rscript --vanilla - "$@" <<'RSCRIPT'
suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) {
    stop("The 'optparse' package is required. Install with install.packages('optparse').")
  }
  if (!requireNamespace("cortar", quietly = TRUE)) {
    stop("The 'cortar' package is required. Install and ensure it's on .libPaths().")
  }
  library(optparse)
  library(cortar)
})

option_list <- list(
  make_option(c("-f", "--file"),       type="character", help="Input file path [required]"),
  make_option(c("-m", "--mode"),       type="character", default="default", help="Mode [default: %default]"),
  make_option(c("-a", "--assembly"),   type="character", default="hg38", help="Genome assembly [default: %default]"),
  make_option(c("-n", "--annotation"), type="character", default="UCSC", help="Annotation source [default: %default]"),
  make_option(c("-t", "--input_type"), type="character", default="bamfile", help="Input type [default: %default]"),
  make_option(c("-p", "--paired"),     type="logical", default=TRUE, help="Paired-end TRUE/FALSE [default: %default]"),
  make_option(c("-s", "--stranded"),   type="integer", default=2, help="Strandedness (0,1,2) [default: %default]"),
  make_option(c("-u", "--subset"),     type="character", default=NULL, help="Comma-separated subset list [default: %default]"),
  make_option(c("-o", "--output_dir"), type="character", default="~", help="Output directory [default: %default]"),
  make_option(c("-g", "--genelist"),   type="character", default=NULL, help="Path to genelist [default: %default]"),
  make_option(c("-x", "--prefix"),     type="character", default="", help="Output prefix [default: %default]"),
  make_option(c("-d", "--debug"),      type="logical", default=FALSE, help="Debug TRUE/FALSE [default: %default]"),
  make_option(c("-r", "--ria"),        type="logical", default=FALSE, help="Reads in absentia TRUE/FALSE [default: %default]")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
opts <- parse_args(parser)

if (is.null(opts$file) || opts$file == "") {
  print_help(parser)
  stop("Error: --file is required.", call.=FALSE)
}

# Normalize
to_null <- function(x) if (is.null(x) || (is.character(x) && x == "")) NULL else x
opts$genelist <- to_null(opts$genelist)

if (!is.null(opts$subset) && nzchar(opts$subset)) {
  items <- trimws(unlist(strsplit(opts$subset, ",")))
  opts$subset <- if (length(items)) items else NULL
} else {
  opts$subset <- NULL
}

opts$paired <- as.logical(opts$paired)
opts$debug  <- as.logical(opts$debug)
opts$ria    <- as.logical(opts$ria)
opts$stranded <- as.integer(opts$stranded)

# ---- Pretty-print parameters ----
cat("\n========================================\n")
cat(" Running cortar() with parameters\n")
cat("========================================\n")

# Create a named list of key → string
params <- list(
  file       = opts$file,
  mode       = opts$mode,
  assembly   = opts$assembly,
  annotation = opts$annotation,
  input_type = opts$input_type,
  paired     = opts$paired,
  stranded   = opts$stranded,
  subset     = if (is.null(opts$subset)) "NULL" else paste(opts$subset, collapse = ","),
  output_dir = opts$output_dir,
  genelist   = if (is.null(opts$genelist)) "NULL" else opts$genelist,
  prefix     = opts$prefix,
  debug      = opts$debug,
  ria        = opts$ria
)

# Align keys nicely
max_key <- max(nchar(names(params)))
for (nm in names(params)) {
  cat(sprintf("%-*s : %s\n", max_key, nm, params[[nm]]))
}
cat("========================================\n\n")

# ---- Run cortar ----
res <- tryCatch({
  cortar::cortar(
    file        = opts$file,
    mode        = opts$mode,
    assembly    = opts$assembly,
    annotation  = opts$annotation,
    input_type  = opts$input_type,
    paired      = opts$paired,
    stranded    = opts$stranded,
    subset      = opts$subset,
    output_dir  = opts$output_dir,
    genelist    = opts$genelist,
    prefix      = opts$prefix,
    debug       = opts$debug,
    ria         = opts$ria
  )
}, error=function(e) {
  message("cortar() failed: ", conditionMessage(e))
  quit(status=1)
})

if (!is.null(res)) {
  message("cortar() completed.")
  cls <- paste(class(res), collapse=",")
  nm  <- tryCatch(names(res), error=function(...) NULL)
  cat("Result class:", cls, "\n")
  if (!is.null(nm)) cat("Result fields:", paste(nm, collapse=", "), "\n")
}
RSCRIPT