#!/usr/bin/env Rscript

# Load necessary library
suppressPackageStartupMessages(library(cortar))
suppressPackageStartupMessages(library(optparse))

# Define command-line options
option_list <- list(
  make_option(c("-g", "--genes"), type = "character", default = NULL,
              help = "Comma-separated list of genes (e.g., EMD,DMD)", metavar = "genes"),
  make_option(c("-H", "--hg"), type = "integer", default = NULL,
              help = "Genome version (e.g., 38)", metavar = "hg"),
  make_option(c("-o", "--overhang"), type = "integer", default = 1000,
              help = "Overhang value [default %default]", metavar = "overhang")
)

# Create option parser
parser <- OptionParser(option_list = option_list,
                       usage = "Usage: subset_bam.R --genes EMD,DMD --hg 38 [--overhang 1000]")

# Parse arguments
opts <- parse_args(parser)

# Validate required arguments
if (is.null(opts$genes)) {
  print_help(parser)
  stop("Error: --genes parameter is required.", call. = FALSE)
}

if (is.null(opts$hg)) {
  print_help(parser)
  stop("Error: --hg parameter is required.", call. = FALSE)
}

# Split genes into a vector
genes <- unlist(strsplit(opts$genes, split = ","))

# Call the subsetBamfiles function
regions <- capture.output(subsetBamfiles(genes = genes, hg = opts$hg, overhang = opts$overhang))

# Reformat regions into a format accepted by samtools
reformatted_regions <- list()

for (region in regions) {
  reformatted_regions[[region]] <- gsub(".+(chr\\w+):(\\d+)-(\\d+).+", "\\1:\\2-\\3", region)
}

# Format each region with double single quotes
output <- paste(lapply(reformatted_regions, function(x) { return(paste0("", x, "")) }), 
                collapse = '\n')

# Print the output
cat(output)

