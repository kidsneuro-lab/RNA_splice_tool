library(testthat)
library(data.table)
library(cortar)

load_compare_splicing_fixtures <- function() {
  list(
    all_splicing_events = fread(test_path("mock_all_splicing_events.tsv")),
    sample_file = fread(test_path("mock_samplefile.tsv"))
  )
}

validate_fixture_alignment <- function(all_splicing_events, sample_file) {
  sample_ids <- sample_file$sampleID
  expect_true(all(paste0("count_", sample_ids) %in% names(all_splicing_events)))
  expect_true(all(paste0("pct_", sample_ids) %in% names(all_splicing_events)))
}

first_comparison_table <- function(comparisons) {
  non_null <- Filter(Negate(is.null), comparisons)
  expect_gt(length(non_null), 0)
  non_null[[1]]
}

# Integration test -------------------------------------------------------------
test_that("compareSplicing returns correct results for aligned mock data", {
  fixtures <- load_compare_splicing_fixtures()
  validate_fixture_alignment(fixtures$all_splicing_events, fixtures$sample_file)

  expect_no_error(
    result <- compareSplicing(
      fixtures$all_splicing_events,
      fixtures$sample_file,
      mode = "default",
      debug = FALSE
    )
  )

  expect_type(result, "list")
  expect_gt(length(Filter(Negate(is.null), result)), 0)
})

# Default mode ----------------------------------------------------------------
test_that("compareSplicing handles 'default' mode correctly", {
  fixtures <- load_compare_splicing_fixtures()
  result <- compareSplicing(
    fixtures$all_splicing_events,
    fixtures$sample_file,
    mode = "default",
    debug = FALSE
  )

  expect_type(result, "list")

  first_result <- first_comparison_table(result)
  expect_true("difference" %in% names(first_result))
  expect_true("controlavg" %in% names(first_result))
  expect_true("controlsd" %in% names(first_result))
})

# Panel mode ------------------------------------------------------------------
test_that("compareSplicing handles 'panel' mode correctly", {
  fixtures <- load_compare_splicing_fixtures()

  expect_no_error(
    result <- compareSplicing(
      fixtures$all_splicing_events,
      fixtures$sample_file,
      mode = "panel",
      debug = FALSE
    )
  )

  expect_type(result, "list")
  expect_gt(length(Filter(Negate(is.null), result)), 0)
})

# Research mode ---------------------------------------------------------------
test_that("compareSplicing handles 'research' mode correctly", {
  fixtures <- load_compare_splicing_fixtures()

  result <- compareSplicing(
    fixtures$all_splicing_events,
    fixtures$sample_file,
    mode = "research",
    debug = FALSE
  )

  expect_s3_class(result, "data.table")
  expect_true("controlavg" %in% names(result))
  expect_true("controlsd" %in% names(result))
})

# Unique events ---------------------------------------------------------------
test_that("compareSplicing reports uniqueness and sorted differences", {
  fixtures <- load_compare_splicing_fixtures()

  result <- compareSplicing(
    fixtures$all_splicing_events,
    fixtures$sample_file,
    mode = "default",
    debug = FALSE
  )

  first_result <- first_comparison_table(result)

  expect_true("unique" %in% names(first_result))

  diffs <- abs(first_result$difference)
  expect_true(all(diff(diffs) <= 0))
})

# Coverage filtering ----------------------------------------------------------
test_that("compareSplicing tracks control counts in default mode", {
  fixtures <- load_compare_splicing_fixtures()

  result <- compareSplicing(
    fixtures$all_splicing_events,
    fixtures$sample_file,
    mode = "default",
    debug = FALSE
  )

  first_result <- first_comparison_table(result)

  expect_true("controln" %in% names(first_result))
  expect_true(all(first_result$controln >= 0))
  expect_true(any(first_result$controln > 0))
})

# Error handling --------------------------------------------------------------
test_that("compareSplicing rejects unsupported modes", {
  fixtures <- load_compare_splicing_fixtures()

  expect_error(
    compareSplicing(
      fixtures$all_splicing_events,
      fixtures$sample_file,
      mode = "invalid_mode",
      debug = FALSE
    ),
    "Invalid mode"
  )
})
