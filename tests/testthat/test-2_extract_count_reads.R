library(testthat)
library(cortar)

test_that("read_sj_sample remaps strand encodings and adjusts IR starts", {
  sj_file <- tempfile(fileext = ".tab")
  ir_file <- tempfile(fileext = ".bed")

  writeLines(
    c(
      "chr1\t100\t200\t0\t0\t0\t11\t0\t0",
      "chr1\t210\t300\t1\t0\t0\t12\t0\t0",
      "chr1\t320\t410\t2\t0\t0\t13\t0\t0"
    ),
    sj_file
  )
  writeLines("chr1\t300\t400\tir_1\t7\t+", ir_file)

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
  expect_equal(GenomicRanges::mcols(sample_reads$ir)$ir_score, 7)
})

test_that("read_sj_sample handles empty SJ and IR files", {
  sj_file <- tempfile(fileext = ".tab")
  ir_file <- tempfile(fileext = ".bed")
  file.create(sj_file)
  file.create(ir_file)

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
