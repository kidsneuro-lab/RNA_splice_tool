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

test_that("extractCountReads does not multiply counts when a gene is listed more than once", {
  # In default mode genes come from file$genes, one entry per sample row, so the
  # same gene appears once for every sample that lists it
  bamfile <- system.file("extdata", "test_sample.bam", package = "cortar")

  genes_tx <- selectGenesTranscripts("LEMD3", "hg38", "UCSC")

  count_reads <- function(genes.GRanges) {
    reads <- suppressMessages(extractCountReads(
      genes.GRanges = genes.GRanges,
      introns.GRanges = genes_tx$introns$granges,
      intron_starts.GRanges = genes_tx$junctions$starts,
      intron_ends.GRanges = genes_tx$junctions$ends,
      bamfiles = bamfile,
      sjfiles = "",
      irfiles = "",
      sample_names = "test_sample",
      assembly = "hg38",
      annotation = "UCSC",
      paired = TRUE,
      stranded = 2,
      input = "bamfile",
      debug = FALSE
    ))
    as.data.table(reads)
  }

  single <- count_reads(genes_tx$genes)
  duplicated <- count_reads(c(genes_tx$genes, genes_tx$genes))

  expect_gt(sum(single$count_test_sample), 0)
  expect_equal(duplicated, single)
})

test_that("selectGenesTranscripts returns one range per gene listed more than once", {
  single <- selectGenesTranscripts("LEMD3", "hg38", "UCSC")
  duplicated <- selectGenesTranscripts(c("LEMD3", "LEMD3"), "hg38", "UCSC")

  expect_equal(duplicated, single)
})
