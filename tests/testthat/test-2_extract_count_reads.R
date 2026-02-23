library(testthat)
library(cortar)

test_that("read_sj_sample remaps strand encodings and adjusts BED coordinates", {
  sj_file <- test_path("data", "input", "mock_sj.out.tab")
  ir_file <- test_path("data", "input", "mock_ir.bed")

  genes.GRanges <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1, end = 1000),
    strand = "*"
  )
  introns.GRanges <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 301, end = 400),
    strand = "+"
  )

  sample_reads <- read_sj_sample(
    sjfile = sj_file,
    irfile = ir_file,
    genes.GRanges = genes.GRanges,
    introns.GRanges = introns.GRanges
  )

  expect_equal(length(sample_reads$sj), 3L)
  expect_equal(as.character(BiocGenerics::strand(sample_reads$sj)), c("*", "+", "-"))
  expect_equal(GenomicRanges::start(sample_reads$ir), 301L)
  expect_equal(GenomicRanges::end(sample_reads$ir), 400L)
  expect_equal(GenomicRanges::mcols(sample_reads$ir)$ir_score, 7)
})

test_that("read_sj_sample handles empty SJ and IR files", {
  sj_file <- test_path("data", "input", "mock_sj_empty.out.tab")
  ir_file <- test_path("data", "input", "mock_ir_empty.bed")

  genes.GRanges <- GenomicRanges::GRanges()
  introns.GRanges <- GenomicRanges::GRanges()

  expect_no_error(
    sample_reads <- read_sj_sample(
      sjfile = sj_file,
      irfile = ir_file,
      genes.GRanges = genes.GRanges,
      introns.GRanges = introns.GRanges
    )
  )

  expect_equal(length(sample_reads$sj), 0L)
  expect_equal(length(sample_reads$ir), 0L)
})
