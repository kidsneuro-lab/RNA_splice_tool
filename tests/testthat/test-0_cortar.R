library(testthat)
library(data.table)
library(cortar)

# NOTE: This is a process test for the `cortar` function.
# It checks that the function runs without errors and produces the expected output.
# It is NOT a unit test for the function's internal logic.

test_that("cortar produces correct output and runs successfully", {
  input_file <- test_path("data/input/samples.tsv")
  output_dir <- test_path("data/output")
  expected_file <- test_path("data/expected/sample1_EMD_EMD_combined_full.tsv")
  output_file <- file.path(output_dir, "sample1_EMD_EMD_combined_full.tsv")

  # Remove existing output file if present
  if (file.exists(output_file)) {
    file.remove(output_file)
  }

  # Run the function
  expect_error(
    cortar(
      file = input_file,
      output_dir = output_dir,
      ria = FALSE
    ),
    NA  # expect no error
  )

  # Check that file is created
  expect_true(file.exists(output_file))

  # Compare contents
  expected <- fread(expected_file, sep = "\t")
  actual <- fread(output_file, sep = "\t")

  expect_equal(actual, expected)
})