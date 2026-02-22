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

test_that("cortar_batch does not depend on caller debug symbol", {
  clean_env <- new.env(parent = baseenv())
  clean_env$cortar_batch <- cortar_batch
  clean_env$batch_dir <- withr::local_tempdir()

  expect_identical(formals(cortar_batch)$debug, FALSE)
  expect_no_error(
    evalq(cortar_batch(folder = batch_dir, pattern = "\\.tsv$"), envir = clean_env)
  )
})
