#!/usr/bin/env Rscript

# Load necessary library
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(cortar))

# Define command-line options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input folder", metavar = "input"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output folder", metavar = "output")
)

# Create option parser
parser <- OptionParser(option_list = option_list,
                       usage = "Usage: process_cortar.R --input <input folder> --output <output folder>")

# Parse arguments
opts <- parse_args(parser)

# -------------------------- Error Checking --------------------------
# Check if input parameter is provided
if (is.null(opts$input)) {
  cat("Error: --input not provided.\n")
  print_help(parser)
  quit(status = 1)
}

# Check if output parameter is provided
if (is.null(opts$output)) {
  cat("Error: --output not provided.\n")
  print_help(parser)
  quit(status = 1)
}

# Check if input directory exists
if (!file.exists(opts$input)) {
  cat("Error: Samples file does not exist:", opts$input, "\n")
  quit(status = 1)
}

# Check if output directory exists; if not, attempt to create it
if (!dir.exists(opts$output)) {
  cat("Output directory does not exist. Attempting to create:", opts$output, "\n")
  dir.create(opts$output, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(opts$output)) {
    cat("Error: Failed to create output directory:", opts$output, "\n")
    quit(status = 1)
  } else {
    cat("Successfully created output directory:", opts$output, "\n")
  }
}

# Check if cortar function is available
if (!exists("cortar")) {
  cat("Error: 'cortar' function is not available. Please ensure 'cortar' package or function is loaded.\n")
  quit(status = 1)
}

# -------------------------- Runtime Output --------------------------
cat("Starting cortar processing...\n")
cat("Input directory:", opts$input, "\n")
cat("Output directory:", opts$output, "\n")
cat("Debug mode: TRUE\n")
cat("RIA mode: TRUE\n")

# -------------------------- Run cortar --------------------------
# Wrap in a tryCatch to handle runtime errors gracefully
result <- tryCatch({
  cortar(
    file = opts$input,
    output_dir = opts$output,
    ria = TRUE
  )
}, error = function(e) {
  cat("Error running cortar:\n", e$message, "\n")
  quit(status = 1)
})

cat("cortar processing completed successfully.\n")