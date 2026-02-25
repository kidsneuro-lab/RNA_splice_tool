library(testthat)

# TX_EXTRACTION ---------------------------------------------------------------

# Mock assembly data based on provided refseq_introns_exons_hg38 data
load(test_path("refseq_introns_exons_hg38_mock.rda"))
refseq_mock <- refseq_introns_exons_hg38_mock

# 1. Basic functionality
test_that("tx_extraction correctly parses gene names", {
  result <- tx_extraction(c("EMD", "COL2A1"), refseq_mock)
  expected_output <- data.table::data.table(
    gene_name = c("EMD", "COL2A1"),
    tx_version_id = c("NM_000117.3", "NM_001844.5")
  )
  expect_equal(result, expected_output)
})

test_that("tx_extraction correctly parses gene names", {
  result <- tx_extraction(c(2010, 1280), refseq_mock)
  expected_output <- data.table::data.table(
    gene_name = c("EMD", "COL2A1"),
    tx_version_id = c("NM_000117.3", "NM_001844.5")
  )
  expect_equal(result, expected_output)
})

test_that("tx_extraction correctly parses gene names", {
  result <- tx_extraction(c(2010, 'NM_033150.3'), refseq_mock)
  expected_output <- data.table::data.table(
    gene_name = c("EMD", "COL2A1"),
    tx_version_id = c("NM_000117.3", "NM_033150.3")
  )
  expect_equal(result, expected_output)
})

# Mix of gene names, transcripts, and blank
test_that("tx_extraction can handle mixed input", {
  result <- tx_extraction(c(2010, 'NM_033150.3',""), refseq_mock)
  expected_output <- data.table::data.table(
    gene_name = c("EMD", "COL2A1"),
    tx_version_id = c("NM_000117.3", "NM_033150.3")
  )
  expect_equal(result, expected_output)
})

# Mismatched or invalid input
test_that("tx_extraction throws an error for invalid gene names", {
  expect_error(tx_extraction(c("NON_EXISTENT"), refseq_mock),
               "Gene or Transcript identifier `NON_EXISTENT` is invalid")
})

test_that("tx_extraction throws an error for invalid transcripts", {
  expect_error(tx_extraction(c("NM_001844.0"), refseq_mock),
               "Gene or Transcript identifier `NM_001844.0` is invalid")
})

# 4. Edge cases
test_that("tx_extraction returns empty data.table for empty input", {
  result <- tx_extraction(character(0), refseq_mock)
  expect_equal(nrow(result), 0)
})

test_that("tx_extraction returns empty data.table for empty input", {
  result <- tx_extraction(c(""), refseq_mock)
  expect_equal(nrow(result), 0)
})

# 5. Assembly Data (assuming refseq_introns_exons_hg38 has a similar structure to the mock data)
# Note: Modify this test to suit the actual data structure of refseq_introns_exons_hg38
test_that("tx_extraction handles different assemblies", {
  result <- tx_extraction(c("EMD", "COL2A1"), refseq_introns_exons_hg38)
  # Modify the expectation based on the expected behavior for this scenario
  expect_true(nrow(result) > 0)
})


# GENE_TO_GRANGE --------------------------------------------------------------

test_that("Obtain granges for EMD and COL2A1 (UCSC annotation)", {
  gene_tx <- data.table::data.table(
    gene_name = c("EMD", "COL2A1"),
    tx_version_id = c("NM_000117.3", "NM_001844.5")
  )

  result = gene_to_GRange(
    gene_tx = gene_tx,
    annotation = 'UCSC',
    refseq_introns_exons = refseq_mock,
  )

  expected_result <- GenomicRanges::GRanges(
    seqnames = c('chrX','chr12'),
    IRanges::IRanges(
      start = c(154380882, 48004237),
      end = c(154381523, 48004476)
    ),
    strand = c('+','-')
  )

  expect_equal(result, expected_result)
})

test_that("Obtain granges for EMD and COL2A1 (Default annotation)", {
  gene_tx <- data.table::data.table(
    gene_name = c("EMD", "COL2A1"),
    tx_version_id = c("NM_000117.3", "NM_001844.5")
  )

  result = gene_to_GRange(
    gene_tx = gene_tx,
    annotation = '1000genomes',
    refseq_introns_exons = refseq_mock,
  )

  expected_result <- GenomicRanges::GRanges(
    seqnames = c('X','12'),
    IRanges::IRanges(
      start = c(154380882, 48004237),
      end = c(154381523, 48004476)
    ),
    strand = c('+','-')
  )

  expect_equal(result, expected_result)
})

test_that("Throw error if attempting to obtain GRanges for empty data.table", {
  gene_tx <- data.table::data.table(
    gene_name = c(),
    tx_version_id = c()
  )

  expect_error(
    gene_to_GRange(
      gene_tx = gene_tx,
      annotation = 'UCSC',
      refseq_introns_exons = refseq_mock,
    ),
    "Cannot obtain GRange for empty data.table"
  )
})

# INTRONS_TO_GRANGE --------------------------------------------------------------

test_that("Obtain intron granges for EMD (UCSC annotation)", {
  gene_tx <- data.table::data.table(
    gene_name = c("EMD"),
    tx_version_id = c("NM_000117.3")
  )

  result = introns_to_GRange(
    gene_tx = gene_tx,
    annotation = 'UCSC',
    refseq_intron_exons = refseq_mock,
  )

  expected_result <- GenomicRanges::GRanges(
    seqnames = c('chrX'),
    IRanges::IRanges(
      start = c(154380368),
      end = c(154380752)
    ),
    strand = c('+'),
    intron_no = 4L,
    gene = "EMD"
  )

  expect_equal(result$granges, expected_result)
})

# OTHER_INTRONS_TO_GRANGE --------------------------------------------------------------

test_that("Obtain other intron granges for EMD (UCSC annotation)", {
  gene_tx <- data.table::data.table(
    gene_name = c("COL2A1"),
    tx_version_id = c("NM_001844.5")
  )

  result = introns_other_tx_to_GRange(
    gene_tx = gene_tx,
    annotation = 'UCSC',
    refseq_intron_exons = refseq_mock,
  )

  expected_result <- GenomicRanges::GRanges(
    seqnames = c('chr12'),
    IRanges::IRanges(
      start = c(47980969),
      end = c(47981342)
    ),
    strand = c('-'),
    tx_id = 'NM_033150.3',
    intron_no = 36L,
    gene = "COL2A1"
  )

  expect_equal(result, expected_result)
})

# introns_jx_to_GRange --------------------------------------------------------------

test_that("Obtain other intron granges for EMD (UCSC annotation)", {
  gene_tx <- data.table::data.table(
    gene_name = c("COL2A1"),
    tx_version_id = c("NM_001844.5")
  )

  # Constants for junction boundaries
  INTRON_JUNCTION_UPSTREAM <- 4
  INTRON_JUNCTION_DOWNSTREAM <- 3

  result = introns_jx_to_GRange(
    gene_tx = gene_tx,
    annotation = 'UCSC',
    refseq_introns_exons = refseq_mock,
  )

  expected_starts_result <- GenomicRanges::GRanges(
    seqnames = c('chr12'),
    IRanges::IRanges(
      start = c(47976071-INTRON_JUNCTION_UPSTREAM),
      end = c(47976071+INTRON_JUNCTION_DOWNSTREAM)
    ),
    strand = c('-')
  )

  expected_ends_result <- GenomicRanges::GRanges(
    seqnames = c('chr12'),
    IRanges::IRanges(
      start = c(47976513-INTRON_JUNCTION_DOWNSTREAM),
      end = c(47976513+INTRON_JUNCTION_UPSTREAM)
    ),
    strand = c('-')
  )

  expect_equal(result$starts, expected_starts_result)
  expect_equal(result$ends, expected_ends_result)
})

# selectGenesTranscripts --------------------------------------------------------------

test_that("selectGenesTranscripts", {
  genes <- c('EMD')
  assembly <- 'hg38'
  annotation = 'UCSC'

  results <- selectGenesTranscripts(
    c('EMD'),
    'hg38',
    'UCSC'
  )

  expected_genes_result <- GenomicRanges::GRanges(
    seqnames = c('chrX'),
    IRanges::IRanges(
      start = c(154379295),
      end = c(154381523)
    ),
    strand = c('+')
  )

  expected_introns_result <- GenomicRanges::GRanges(
    seqnames = "chrX",
    ranges = IRanges::IRanges(
      start = c(154379567, 154379795, 154380020, 154380368, 154380803),
      end   = c(154379689, 154379941, 154380233, 154380752, 154380881)
    ),
    strand = "+",
    intron_no = 1:5,
    gene = "EMD"
  )

  expect_equal(results$genes, expected_genes_result)
  expect_equal(results$introns$granges, expected_introns_result)
  expect_equal(length(results$introns_other_tx), 0)
  expect_equal(length(results$junctions$starts), 5)
  expect_equal(length(results$junctions$ends), 5)
})
