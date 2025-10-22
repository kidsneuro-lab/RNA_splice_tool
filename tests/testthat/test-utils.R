library(testthat)
library(cortar)

# Test is_debug_enabled() -----------------------------------------------------

test_that("is_debug_enabled returns TRUE for valid debug paths", {
  expect_true(is_debug_enabled("/path/to/debug"))
  expect_true(is_debug_enabled("debug"))
  expect_true(is_debug_enabled("./output/debug"))
})

test_that("is_debug_enabled returns FALSE for invalid debug values", {
  expect_false(is_debug_enabled(FALSE))
  expect_false(is_debug_enabled(""))
  expect_false(is_debug_enabled(TRUE))
  expect_false(is_debug_enabled(NULL))
})

# Test get_coverage_threshold() -----------------------------------------------

test_that("get_coverage_threshold returns correct values for standard types", {
  expect_equal(get_coverage_threshold("het"), 60)
  expect_equal(get_coverage_threshold("hom"), 30)
  expect_equal(get_coverage_threshold("hemi"), 30)
  expect_equal(get_coverage_threshold(""), 0)
})

test_that("get_coverage_threshold handles custom numeric values", {
  expect_equal(get_coverage_threshold("50"), 50)
  expect_equal(get_coverage_threshold("100"), 100)
  expect_equal(get_coverage_threshold("15"), 15)
})

# Test %nin% operator ---------------------------------------------------------

test_that("%nin% operator works correctly", {
  expect_true("a" %nin% c("b", "c", "d"))
  expect_false("a" %nin% c("a", "b", "c"))
  expect_true(1 %nin% c(2, 3, 4))
  expect_false(1 %nin% c(1, 2, 3))
})

test_that("%nin% handles edge cases", {
  expect_true("x" %nin% character(0))
  expect_false(NA %nin% c(NA, "a", "b"))
  expect_true("test" %nin% NULL)
})

# Test validate_parameter() ---------------------------------------------------

test_that("validate_parameter works for allowed values", {
  expect_error(
    validate_parameter("mode", "invalid", allowed_values = c("default", "panel", "research")),
    "Parameter 'mode' has invalid value 'invalid'"
  )
  expect_silent(validate_parameter("mode", "default", allowed_values = c("default", "panel", "research")))
  expect_silent(validate_parameter("assembly", "hg38", allowed_values = c("hg19", "hg38")))
})

test_that("validate_parameter checks for NULL and empty values", {
  expect_error(
    validate_parameter("file", NULL),
    "Parameter 'file' is required but was not provided"
  )
  expect_error(
    validate_parameter("file", ""),
    "Parameter 'file' is required but was not provided"
  )
})

test_that("validate_parameter checks file existence", {
  # Test with non-existent file
  expect_error(
    validate_parameter("file", "/nonexistent/file.txt", check_file_exists = TRUE),
    "File '/nonexistent/file.txt' does not exist"
  )
})

# Test validate_bam_files() ---------------------------------------------------

test_that("validate_bam_files detects missing files", {
  expect_error(
    validate_bam_files(c("/nonexistent/file1.bam", "/nonexistent/file2.bam")),
    "The following BAM files do not exist"
  )
})

# Test create_granges() -------------------------------------------------------

test_that("create_granges creates basic GRanges correctly", {
  gr <- create_granges(
    seqnames = c("chr1", "chr2"),
    starts = c(100, 200),
    ends = c(200, 300),
    strands = c("+", "-"),
    annotation = "1000genomes"
  )
  
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 2)
  expect_equal(as.character(GenomicRanges::seqnames(gr)), c("chr1", "chr2"))
  expect_equal(BiocGenerics::start(gr), c(100, 200))
  expect_equal(BiocGenerics::end(gr), c(200, 300))
})

test_that("create_granges handles metadata correctly", {
  gr <- create_granges(
    seqnames = c("chr1"),
    starts = c(100),
    ends = c(200),
    strands = c("+"),
    annotation = "1000genomes",
    metadata = list(gene = "DMD", intron_no = 1)
  )
  
  expect_s4_class(gr, "GRanges")
  expect_true("gene" %in% names(GenomicRanges::mcols(gr)))
  expect_true("intron_no" %in% names(GenomicRanges::mcols(gr)))
  expect_equal(GenomicRanges::mcols(gr)$gene, "DMD")
  expect_equal(GenomicRanges::mcols(gr)$intron_no, 1)
})

test_that("create_granges sets UCSC seqlevel style correctly", {
  gr <- create_granges(
    seqnames = c("1"),
    starts = c(100),
    ends = c(200),
    strands = c("+"),
    annotation = "UCSC"
  )
  
  expect_s4_class(gr, "GRanges")
  # After setting UCSC style, chromosome names should have "chr" prefix
  expect_match(as.character(GenomicRanges::seqnames(gr))[1], "chr", ignore.case = TRUE)
})
