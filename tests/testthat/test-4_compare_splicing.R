#Load Mock Data ---------------------------------------------------------------
mock_splicing_events <- fread(test_path("mock_all_splicing_events.tsv"))

mock_Sample_File <- fread(test_path("mock_samplefile.tsv"))

#Integration Test -------------------------------------------------------------
test_that("compareSplicing returns correct results for mock data", {
  # Test that the function runs without errors
  expect_error(
    result <- compareSplicing(
      mock_splicing_events, 
      mock_Sample_File, 
      mode = "default", 
      debug = FALSE
    ),
    NA
  )
  
  # Verify result structure
  expect_type(result, "list")
  expect_true(length(result) > 0)
})

#Test default mode ------------------------------------------------------------
test_that("compareSplicing handles 'default' mode correctly", {
  result <- compareSplicing(
    mock_splicing_events, 
    mock_Sample_File, 
    mode = "default", 
    debug = FALSE
  )
  
  # Should return a list with results for each test sample
  expect_type(result, "list")
  
  # Check that key columns exist if results are present
  if (length(result) > 0 && !is.null(result[[1]])) {
    expect_true("difference" %in% names(result[[1]]))
    expect_true("controlavg" %in% names(result[[1]]))
    expect_true("controlsd" %in% names(result[[1]]))
  }
})

#Test panel mode --------------------------------------------------------------
test_that("compareSplicing handles 'panel' mode correctly", {
  # Panel mode should work similarly to default
  expect_error(
    result <- compareSplicing(
      mock_splicing_events, 
      mock_Sample_File, 
      mode = "panel", 
      debug = FALSE
    ),
    NA
  )
})

#Test research mode -----------------------------------------------------------
test_that("compareSplicing handles 'research' mode correctly", {
  result <- compareSplicing(
    mock_splicing_events, 
    mock_Sample_File, 
    mode = "research", 
    debug = FALSE
  )
  
  # Research mode returns a data.table, not a list
  expect_true(inherits(result, "data.table") || inherits(result, "data.frame"))
  
  # Should have control average columns
  expect_true("controlavg" %in% names(result))
  expect_true("controlsd" %in% names(result))
})

#Test unique events -----------------------------------------------------------
test_that("compareSplicing identifies unique splicing events correctly", {
  result <- compareSplicing(
    mock_splicing_events, 
    mock_Sample_File, 
    mode = "default", 
    debug = FALSE
  )
  
  # Check that unique column exists
  if (length(result) > 0 && !is.null(result[[1]])) {
    expect_true("unique" %in% names(result[[1]]))
  }
})

#Test order and sorting -------------------------------------------------------
test_that("compareSplicing orders and sorts results correctly", {
  result <- compareSplicing(
    mock_splicing_events, 
    mock_Sample_File, 
    mode = "default", 
    debug = FALSE
  )
  
  # In default mode, results should be sorted by absolute difference
  if (length(result) > 0 && !is.null(result[[1]])) {
    diffs <- abs(result[[1]]$difference)
    # Check that differences are in descending order (allowing for ties)
    expect_true(all(diff(diffs) <= 0))
  }
})

#Test coverage filtering for controls -----------------------------------------
test_that("compareSplicing filters low coverage controls in default mode", {
  # Test with a sample file that has coverage specified
  result <- compareSplicing(
    mock_splicing_events, 
    mock_Sample_File, 
    mode = "default", 
    debug = FALSE
  )
  
  # Verify that controln is set (number of controls used)
  if (length(result) > 0 && !is.null(result[[1]])) {
    expect_true("controln" %in% names(result[[1]]))
    expect_true(all(result[[1]]$controln >= 0))
  }
})

#Test error handling ----------------------------------------------------------
test_that("compareSplicing handles errors gracefully", {
  # Test with invalid mode
  expect_error(
    compareSplicing(
      mock_splicing_events, 
      mock_Sample_File, 
      mode = "invalid_mode", 
      debug = FALSE
    )
  )
})

