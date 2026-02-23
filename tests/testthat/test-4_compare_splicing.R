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

test_that("calculate_unique_events identifies unique events correctly", {
  dt <- data.table(
    controlavg = c(0, 0.5, 0),
    pct_sample1 = c(0.3, 0.4, 0),
    pct_sample2 = c(0.2, 0.1, 0)
  )

  result <- calculate_unique_events(dt, c("pct_sample1", "pct_sample2"))
  expect_equal(result, c("2/2", "", ""))
})

test_that("calculate_unique_events handles empty data.table", {
  dt <- data.table(controlavg = numeric(0), pct_s1 = numeric(0))
  result <- calculate_unique_events(dt, "pct_s1")
  expect_equal(result, character(0))
})

test_that("calculate_unique_events handles partial family uniqueness", {
  dt <- data.table(
    controlavg = c(0),
    pct_sample1 = c(0.3),
    pct_sample2 = c(0)
  )

  result <- calculate_unique_events(dt, c("pct_sample1", "pct_sample2"))
  expect_equal(result, "1/2")
})

test_that("calculate_unique_events ignores NA family values", {
  dt <- data.table(
    controlavg = c(0, 0, NA_real_),
    pct_sample1 = c(NA_real_, 0.4, 0.2),
    pct_sample2 = c(0.1, NA_real_, 0.1)
  )

  result <- calculate_unique_events(dt, c("pct_sample1", "pct_sample2"))
  expect_equal(result, c("1/2", "1/2", ""))
})

test_that("filter_controls_by_coverage removes low-coverage controls", {
  dt <- data.table(
    gene = c("EMD", "EMD", "EMD", "DMD"),
    annotated = rep("canonical", 4),
    SJ_IR = rep("SJ", 4),
    count_ctrl1 = c(100, 120, 110, 5),
    count_ctrl2 = c(10, 5, 8, 500)
  )

  result <- filter_controls_by_coverage(
    events_dt = dt,
    gene_name = "EMD",
    ctrl_pct_cols = c("pct_ctrl1", "pct_ctrl2"),
    ctrl_read_cols = c("count_ctrl1", "count_ctrl2"),
    coverage_type = "het"
  )

  expect_equal(result$ctrl_pct_cols, "pct_ctrl1")
  expect_equal(result$ctrl_read_cols, "count_ctrl1")
})

test_that("filter_controls_by_coverage keeps all controls when threshold is 0", {
  dt <- data.table(
    gene = rep("EMD", 2),
    annotated = rep("canonical", 2),
    SJ_IR = rep("SJ", 2),
    count_ctrl1 = c(5, 3),
    count_ctrl2 = c(2, 1)
  )

  result <- filter_controls_by_coverage(
    events_dt = dt,
    gene_name = "EMD",
    ctrl_pct_cols = c("pct_ctrl1", "pct_ctrl2"),
    ctrl_read_cols = c("count_ctrl1", "count_ctrl2"),
    coverage_type = ""
  )

  expect_equal(length(result$ctrl_pct_cols), 2)
  expect_equal(length(result$ctrl_read_cols), 2)
})
