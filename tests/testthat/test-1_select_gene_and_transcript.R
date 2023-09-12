library(testthat)

# TX_EXTRACTION ---------------------------------------------------------------

# Mock assembly data based on provided refseq_introns_exons_hg38 data
load(testthat::test_path("refseq_intron_exons_hg38_mock.rda"))
refseq_mock <- refseq_intron_exons_hg38_mock

# Mock gene data
mock_genes <- c("EMD", "COL2A1")

# 1. Basic functionality
test_that("tx_extraction correctly parses gene names", {
  result <- tx_extraction(mock_genes, refseq_mock)
  expected_output <- data.table::data.table(
    gene_name = c("EMD", "COL2A1"),
    tx = c("NM_000117.3", "NM_001844.5")
  )
  expect_equal(result, expected_output)
})

# 2. Mismatched or invalid input
test_that("tx_extraction throws an error for invalid gene names", {
  expect_error(tx_extraction(c("NON_EXISTENT"), refseq_mock),
               "Gene name `NON_EXISTENT` is invalid")
})

test_that("tx_extraction throws an error for invalid transcripts", {
  expect_error(tx_extraction(c("NM_001844.0"), refseq_mock),
               "Transcript identifier `NM_001844.0` is invalid")
})

# 3. Mix of gene names, transcripts, and blank
test_that("tx_extraction can handle mixed input", {
  mixed_input <- c("EMD", "NM_033150.3", "")
  result <- tx_extraction(mixed_input, refseq_mock)
  expected_output <- data.table::data.table(
    gene_name = c("EMD", "COL2A1"),
    tx = c("NM_000117.3", "NM_033150.3")
  )
  expect_equal(result, expected_output)
})

# 4. Edge cases
test_that("tx_extraction returns empty data.table for empty input", {
  result <- tx_extraction(character(0), refseq_mock)
  expect_equal(nrow(result), 0)
})

test_that("tx_extraction returns empty data.table for blank input", {
  result <- tx_extraction("", refseq_mock)
  expect_equal(nrow(result), 0)
})

# 5. Assembly Data (assuming refseq_introns_exons_hg38 has a similar structure to the mock data)
# Note: Modify this test to suit the actual data structure of refseq_introns_exons_hg38
test_that("tx_extraction handles different assemblies", {
  result <- tx_extraction(mock_genes, refseq_introns_exons_hg38)
  # Modify the expectation based on the expected behavior for this scenario
  expect_true(nrow(result) > 0)
})


# GENE_TO_GRANGE --------------------------------------------------------------
