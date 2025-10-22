library(testthat)
library(cortar)

# Test edge cases for tx_extraction ---------------------------------------------

test_that("tx_extraction handles invalid gene names", {
  expect_error(
    tx_extraction("INVALIDGENE999", refseq_introns_exons_hg38),
    "Gene name `INVALIDGENE999` is invalid"
  )
})

test_that("tx_extraction handles invalid transcript IDs", {
  expect_error(
    tx_extraction("NM_999999.99", refseq_introns_exons_hg38),
    "Transcript identifier `NM_999999.99` is invalid"
  )
})

test_that("tx_extraction handles empty input", {
  expect_error(
    tx_extraction(character(0), refseq_introns_exons_hg38)
  )
})

test_that("tx_extraction handles mixed gene and transcript input", {
  # This should work with valid entries
  result <- tx_extraction(
    c("DMD", "NM_004006.3"),
    refseq_introns_exons_hg38,
    debug = FALSE
  )
  
  expect_s3_class(result, "data.table")
  expect_equal(nrow(result), 2)
  expect_true(all(c("gene_name", "tx") %in% names(result)))
})

# Test edge cases for GRanges operations ---------------------------------------

test_that("create_granges handles single-element input", {
  gr <- create_granges(
    seqnames = "chr1",
    starts = 100,
    ends = 200,
    strands = "+",
    annotation = "1000genomes"
  )
  
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 1)
})

test_that("create_granges handles empty metadata list", {
  gr <- create_granges(
    seqnames = "chr1",
    starts = 100,
    ends = 200,
    strands = "+",
    annotation = "1000genomes",
    metadata = list()
  )
  
  expect_s4_class(gr, "GRanges")
  expect_equal(length(GenomicRanges::mcols(gr)), 0)
})

test_that("create_granges handles multiple metadata columns", {
  gr <- create_granges(
    seqnames = c("chr1", "chr2"),
    starts = c(100, 200),
    ends = c(200, 300),
    strands = c("+", "-"),
    annotation = "1000genomes",
    metadata = list(
      gene = c("GENE1", "GENE2"),
      intron_no = c(1, 2),
      tx_id = c("NM_001", "NM_002")
    )
  )
  
  expect_s4_class(gr, "GRanges")
  expect_equal(length(GenomicRanges::mcols(gr)), 3)
  expect_true(all(c("gene", "intron_no", "tx_id") %in% names(GenomicRanges::mcols(gr))))
})

# Test edge cases for coverage threshold ---------------------------------------

test_that("get_coverage_threshold handles boundary values", {
  expect_equal(get_coverage_threshold("0"), 0)
  expect_equal(get_coverage_threshold("1"), 1)
  expect_equal(get_coverage_threshold("1000"), 1000)
})

test_that("get_coverage_threshold handles decimal values", {
  expect_equal(get_coverage_threshold("50.5"), 50.5)
  expect_equal(get_coverage_threshold("99.99"), 99.99)
})

# Test edge cases for parameter validation -------------------------------------

test_that("validate_parameter handles multiple allowed values", {
  expect_silent(
    validate_parameter("mode", "default", 
                      allowed_values = c("default", "panel", "research"))
  )
  expect_silent(
    validate_parameter("mode", "panel", 
                      allowed_values = c("default", "panel", "research"))
  )
  expect_silent(
    validate_parameter("mode", "research", 
                      allowed_values = c("default", "panel", "research"))
  )
})

test_that("validate_parameter provides helpful error for invalid value", {
  expect_error(
    validate_parameter("assembly", "hg37", allowed_values = c("hg19", "hg38")),
    "Parameter 'assembly' has invalid value 'hg37'. Allowed values are: hg19, hg38"
  )
})

test_that("validate_parameter handles whitespace-only strings", {
  expect_error(
    validate_parameter("file", "   "),
    "Parameter 'file' is required"
  )
})

# Test edge cases for BAM file validation --------------------------------------

test_that("validate_bam_files handles empty vector", {
  expect_silent(validate_bam_files(character(0)))
})

test_that("validate_bam_files provides detailed error for multiple missing files", {
  error_msg <- tryCatch(
    validate_bam_files(c("/missing1.bam", "/missing2.bam", "/missing3.bam")),
    error = function(e) e$message
  )
  
  expect_match(error_msg, "missing1.bam")
  expect_match(error_msg, "missing2.bam")
  expect_match(error_msg, "missing3.bam")
})

# Test edge cases for debug mode -----------------------------------------------

test_that("is_debug_enabled handles various input types", {
  expect_false(is_debug_enabled(0))
  expect_false(is_debug_enabled(1))
  expect_false(is_debug_enabled(NA))
  expect_false(is_debug_enabled(list()))
})

test_that("is_debug_enabled handles path edge cases", {
  expect_true(is_debug_enabled("/"))
  expect_true(is_debug_enabled("."))
  expect_true(is_debug_enabled(".."))
  expect_true(is_debug_enabled("a"))
})

# Test %nin% operator edge cases -----------------------------------------------

test_that("%nin% handles mixed types", {
  expect_true(1 %nin% c("1", "2", "3"))
  expect_false("1" %nin% c("1", "2", "3"))
})

test_that("%nin% handles NA values correctly", {
  expect_true(NA %nin% c(1, 2, 3))
  expect_true(1 %nin% c(NA, 2, 3))  # NA doesn't match 1
})

test_that("%nin% handles special numeric values", {
  expect_true(Inf %nin% c(1, 2, 3))
  expect_true(-Inf %nin% c(1, 2, 3))
  expect_false(Inf %nin% c(1, 2, Inf))
})

# Test edge cases for constants ------------------------------------------------

test_that("coverage constants are positive", {
  expect_true(DEFAULT_COVERAGE_HET > 0)
  expect_true(DEFAULT_COVERAGE_HOM_HEMI > 0)
  expect_true(DEFAULT_COVERAGE_NONE >= 0)
})

test_that("junction constants are positive", {
  expect_true(INTRON_JUNCTION_UPSTREAM > 0)
  expect_true(INTRON_JUNCTION_DOWNSTREAM > 0)
})

test_that("coverage threshold hierarchy is logical", {
  expect_true(DEFAULT_COVERAGE_HET > DEFAULT_COVERAGE_HOM_HEMI)
  expect_true(DEFAULT_COVERAGE_HOM_HEMI > DEFAULT_COVERAGE_NONE)
})
