`%nin%` <- Negate(`%in%`)

#' Check if debug mode is enabled
#'
#' @param debug Debug parameter (can be a path string or FALSE)
#' @return TRUE if debug mode is enabled, FALSE otherwise
#' @keywords internal
is_debug_enabled <- function(debug) {
  return(is.character(debug) && nchar(debug) > 0)
}
