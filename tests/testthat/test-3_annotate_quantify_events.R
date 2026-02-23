library(testthat)
library(data.table)
library(cortar)

test_that("eventAnnotation labels canonical splicing correctly", {
  dt <- data.table(
    intron_jxn_start = 5,
    intron_jxn_end = 5,
    strand = "+",
    SJ_IR = "SJ"
  )

  expect_equal(eventAnnotation(dt), "canonical exon 5-6 splicing")
})

test_that("eventAnnotation labels intron retention correctly", {
  dt <- data.table(
    intron_jxn_start = 3,
    intron_jxn_end = 3,
    strand = "+",
    SJ_IR = "IR"
  )

  expect_equal(eventAnnotation(dt), "intron 3 retention")
})

test_that("eventAnnotation labels exon skipping correctly", {
  dt <- data.table(
    intron_jxn_start = 2,
    intron_jxn_end = 5,
    strand = "+",
    SJ_IR = "SJ"
  )

  expect_equal(eventAnnotation(dt), "exon 3-4-5 skipping")
})

test_that("eventAnnotation labels cryptic splice sites correctly", {
  dt <- data.table(
    intron_jxn_start = c(3, 3, NA, NA),
    intron_jxn_end = c(NA, NA, 4, 4),
    strand = c("+", "-", "+", "-"),
    SJ_IR = rep("SJ", 4)
  )

  expected <- c(
    "exon 3 ~ cryptic acceptor",
    "cryptic donor ~ exon 4",
    "cryptic donor ~ exon 5",
    "exon 4 ~ cryptic acceptor"
  )

  expect_equal(eventAnnotation(dt), expected)
})

test_that("eventAnnotation handles empty data.table", {
  dt <- data.table(
    intron_jxn_start = integer(0),
    intron_jxn_end = integer(0),
    strand = character(0),
    SJ_IR = character(0)
  )

  expect_equal(eventAnnotation(dt), character(0))
})

test_that("framed returns TRUE for in-frame canonical junctions", {
  rfsq <- refseq_introns_exons_hg38[
    !is.na(region_start) & !is.na(region_end),
  ]
  expect_gt(nrow(rfsq), 0)

  dt <- data.table(
    SJ_IR = "SJ",
    start = rfsq$region_start[1],
    end = rfsq$region_end[1]
  )

  expect_equal(framed(dt, "hg38"), "TRUE")
})

test_that("framed returns empty string for intron retention events", {
  dt <- data.table(SJ_IR = "IR", start = 100, end = 200)

  expect_equal(framed(dt, "hg38"), "")
})

test_that("framed returns NA for unannotated junctions", {
  dt <- data.table(SJ_IR = "SJ", start = -999, end = -998)
  result <- framed(dt, "hg38")

  expect_true(is.na(result))
})

test_that("framed handles empty data.table", {
  dt <- data.table(SJ_IR = character(0), start = integer(0), end = integer(0))

  expect_equal(framed(dt, "hg38"), character(0))
})
