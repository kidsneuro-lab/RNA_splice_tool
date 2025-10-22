`%nin%` <- Negate(`%in%`)

# Constants for coverage thresholds
DEFAULT_COVERAGE_HET <- 60
DEFAULT_COVERAGE_HOM_HEMI <- 30
DEFAULT_COVERAGE_NONE <- 0

# Constants for junction boundaries
INTRON_JUNCTION_UPSTREAM <- 4
INTRON_JUNCTION_DOWNSTREAM <- 3

# Cache for genome assemblies (improves performance when running multiple analyses)
.genome_cache <- new.env(parent = emptyenv())

#' Check if debug mode is enabled
#'
#' @param debug Debug parameter (can be a path string or FALSE)
#' @return TRUE if debug mode is enabled, FALSE otherwise
#' @keywords internal
is_debug_enabled <- function(debug) {
  return(is.character(debug) && nchar(debug) > 0)
}

#' Get Coverage Threshold
#'
#' Internal helper to determine the appropriate coverage threshold based on
#' the coverage type (heterozygous, homozygous, hemizygous, or custom value).
#'
#' @param coverage_type Character string indicating coverage type ("het", "hom",
#'   "hemi", "") or a numeric string
#'
#' @return Numeric coverage threshold value
#' @keywords internal
get_coverage_threshold <- function(coverage_type) {
  if (coverage_type == "het") {
    return(DEFAULT_COVERAGE_HET)
  } else if (coverage_type %in% c("hom", "hemi")) {
    return(DEFAULT_COVERAGE_HOM_HEMI)
  } else if (coverage_type == "") {
    return(DEFAULT_COVERAGE_NONE)
  } else {
    return(as.numeric(coverage_type))
  }
}

#' Create a GRanges Object from Data Table
#'
#' Internal helper function to create a GRanges object from genomic coordinates.
#' This reduces code duplication across multiple functions that create GRanges.
#'
#' @param seqnames Character vector or column of chromosome names
#' @param starts Numeric vector of start positions
#' @param ends Numeric vector of end positions
#' @param strands Character vector of strand information ("+", "-", or "*")
#' @param annotation Annotation type ("UCSC" or "1000genomes")
#' @param metadata Optional list of metadata columns to add to the GRanges object
#'
#' @return A GRanges object with the specified coordinates and optional metadata
#' @keywords internal
create_granges <- function(seqnames, starts, ends, strands, annotation, metadata = NULL) {
  gr <- GenomicRanges::GRanges(
    seqnames = seqnames,
    IRanges::IRanges(start = starts, end = ends),
    strand = strands
  )
  
  # Set seqlevels style based on annotation
  if (annotation == "UCSC") {
    GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
  }
  
  # Add metadata if provided
  if (!is.null(metadata)) {
    for (col_name in names(metadata)) {
      GenomicRanges::mcols(gr)[col_name] <- metadata[[col_name]]
    }
  }
  
  return(gr)
}

#' Validate Required Parameters
#'
#' Internal helper to validate that required parameters are present and valid.
#' Throws informative errors if validation fails.
#'
#' @param param_name Name of the parameter being validated
#' @param param_value Value of the parameter
#' @param allowed_values Optional character vector of allowed values
#' @param check_file_exists If TRUE, checks that the parameter is a valid file path
#' @param check_dir_exists If TRUE, checks that the parameter is a valid directory path
#'
#' @return NULL if validation passes, throws error otherwise
#' @keywords internal
validate_parameter <- function(param_name, param_value, 
                               allowed_values = NULL, 
                               check_file_exists = FALSE,
                               check_dir_exists = FALSE) {
  # Check for NULL or missing
  if (is.null(param_value) || (is.character(param_value) && nchar(param_value) == 0)) {
    stop(sprintf("Parameter '%s' is required but was not provided or is empty.", param_name))
  }
  
  # Check allowed values
  if (!is.null(allowed_values) && !(param_value %in% allowed_values)) {
    stop(sprintf(
      "Parameter '%s' has invalid value '%s'. Allowed values are: %s",
      param_name,
      param_value,
      paste(allowed_values, collapse = ", ")
    ))
  }
  
  # Check file exists
  if (check_file_exists && !file.exists(param_value)) {
    stop(sprintf(
      "Parameter '%s': File '%s' does not exist or is not readable. Current working directory: '%s'",
      param_name,
      param_value,
      getwd()
    ))
  }
  
  # Check directory exists
  if (check_dir_exists && !dir.exists(param_value)) {
    stop(sprintf(
      "Parameter '%s': Directory '%s' does not exist. Current working directory: '%s'",
      param_name,
      param_value,
      getwd()
    ))
  }
  
  return(invisible(NULL))
}

#' Validate BAM Files
#'
#' Internal helper to validate that BAM files exist and are readable.
#' Provides clear error messages indicating which files are problematic.
#'
#' @param bamfiles Character vector of BAM file paths
#'
#' @return NULL if validation passes, throws error otherwise
#' @keywords internal
validate_bam_files <- function(bamfiles) {
  missing_files <- character(0)
  
  for (bamfile in bamfiles) {
    if (!file.exists(bamfile)) {
      missing_files <- c(missing_files, bamfile)
    }
  }
  
  if (length(missing_files) > 0) {
    stop(sprintf(
      "The following BAM files do not exist or are not readable:\n%s\nCurrent working directory: '%s'",
      paste("  -", missing_files, collapse = "\n"),
      getwd()
    ))
  }
  
  return(invisible(NULL))
}

#' Get Genome Assembly (with caching)
#'
#' Internal helper to retrieve genome assembly objects with caching for
#' improved performance when running multiple analyses.
#'
#' @param assembly Genome assembly version ("hg38" or "hg19")
#' @param annotation Annotation type ("UCSC", "1000genomes", or "NCBI")
#'
#' @return BSgenome object for the specified assembly and annotation
#' @keywords internal
get_genome_assembly <- function(assembly, annotation) {
  cache_key <- paste(assembly, annotation, sep = "_")
  
  # Check if already cached
  if (exists(cache_key, envir = .genome_cache)) {
    return(get(cache_key, envir = .genome_cache))
  }
  
  # Load genome assembly
  genome <- NULL
  if (assembly == "hg19") {
    if (annotation == "UCSC") {
      genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    } else if (annotation == "1000genomes") {
      genome <- BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
    }
  } else if (assembly == "hg38") {
    if (annotation == "UCSC") {
      genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    } else if (annotation == "NCBI") {
      genome <- BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
    }
  }
  
  if (is.null(genome)) {
    stop(sprintf(
      "Invalid assembly/annotation combination: %s/%s",
      assembly, annotation
    ))
  }
  
  # Cache and return
  assign(cache_key, genome, envir = .genome_cache)
  return(genome)
}

#' Clear Genome Assembly Cache
#'
#' Internal helper to clear the genome assembly cache. Useful for memory
#' management when processing large numbers of samples.
#'
#' @return NULL (invisibly)
#' @keywords internal
clear_genome_cache <- function() {
  rm(list = ls(envir = .genome_cache), envir = .genome_cache)
  gc()  # Force garbage collection
  return(invisible(NULL))
}
