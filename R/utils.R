`%nin%` <- Negate(`%in%`)

# Constants for coverage thresholds
DEFAULT_COVERAGE_HET <- 60
DEFAULT_COVERAGE_HOM_HEMI <- 30
DEFAULT_COVERAGE_NONE <- 0

# Constants for junction boundaries
INTRON_JUNCTION_UPSTREAM <- 4
INTRON_JUNCTION_DOWNSTREAM <- 3

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
