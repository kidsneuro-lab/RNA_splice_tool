#Load Mock Data ---------------------------------------------------------------
mock_splicing_events <- fread("tests/testthat/mock_all_splicing_events.tsv")

mock_Sample_File <- fread("tests/testthat/mock_samplefile.tsv")

#Integration Test -------------------------------------------------------------
test_that("compareSplicing returns correct results for mock data", {
  # Here you'll want to compare the result to a known expected output.
  # This might look like:
  # expect_equal(result, expected_output)
})

#Test default mode ------------------------------------------------------------
test_that("compareSplicing handles 'default' mode correctly", {
  # Add assertions specific to the 'default' mode
})

#Test panel mode --------------------------------------------------------------
test_that("compareSplicing handles 'panel' mode correctly", {
  # Add assertions specific to the 'panel' mode
})

#Test research mode -----------------------------------------------------------
test_that("compareSplicing handles 'research' mode correctly", {
  # Add assertions specific to the 'research' mode
})

#Test unique events -----------------------------------------------------------
test_that("compareSplicing identifies unique splicing events correctly", {
  # Here you'll want to create a specific mock data where you know which events should be marked as unique
  # Then you'll run the function and check that the 'unique' column has the expected values.
})

#Test order and sorting -------------------------------------------------------
test_that("compareSplicing orders and sorts results correctly", {
  # Check that the returned data is in the correct order.
})

#Test coverage filtering for controls -----------------------------------------
test_that("compareSplicing removes low coverage controls", {
  # Check that the returned data is in the correct order.
})

#Test error handling ----------------------------------------------------------
test_that("compareSplicing handles errors gracefully", {
  # Here you can use things like expect_error() to check that the function fails in a controlled manner for bad input
})

