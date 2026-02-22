library(testthat)
library(data.table)
library(cortar)

# NOTE: This is an integration test for the full cortar() pipeline.

test_that("cortar produces correct output and runs successfully", {
  input_file <- prepare_samplefile_fixture("data/input/samples_default.tsv")
  output_dir <- withr::local_tempdir()
  expected_file <- test_path("data/expected/sample1_EMD_EMD_combined_full.tsv")
  output_file <- file.path(output_dir, "sample1_EMD_EMD_combined_full.tsv")

  expect_no_error(
    suppressWarnings(cortar(
      file = input_file,
      output_dir = output_dir,
      ria = FALSE
    ))
  )

  expect_true(file.exists(output_file))

  expected <- fread(expected_file, sep = "\t")
  actual <- fread(output_file, sep = "\t")

  expect_equal(actual, expected)
})
