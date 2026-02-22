resolve_fixture_path <- function(path_value) {
  if (is.na(path_value) || !nzchar(path_value)) {
    return(path_value)
  }

  candidates <- c(
    path_value,
    test_path(path_value),
    file.path(test_path("data/input"), basename(path_value))
  )

  existing <- candidates[file.exists(candidates)]
  if (!length(existing)) {
    stop(
      "Unable to resolve fixture path '", path_value,
      "' from test root '", test_path(), "'."
    )
  }

  normalizePath(existing[[1]], winslash = "/", mustWork = TRUE)
}

prepare_samplefile_fixture <- function(samplefile_relpath) {
  sample_file <- test_path(samplefile_relpath)
  sample_dt <- data.table::fread(sample_file)

  if (!"family" %in% names(sample_dt) && "familyID" %in% names(sample_dt)) {
    sample_dt[, family := familyID]
  }
  if (!"gene" %in% names(sample_dt) && "genes" %in% names(sample_dt)) {
    sample_dt[, gene := genes]
  }

  if ("bamfile" %in% names(sample_dt)) {
    sample_dt[, bamfile := vapply(bamfile, resolve_fixture_path, character(1))]
  }

  samplefile_tmp <- tempfile(pattern = "samplefile_", fileext = ".tsv")
  data.table::fwrite(sample_dt, samplefile_tmp, sep = "\t")
  samplefile_tmp
}
