library(testthat)
library(data.table)
library(cortar)

make_mock_comparison <- function(sample_id = "sample1", gene_name = "GENE1") {
  comparison <- data.table(
    chr = c("chr1", "chr1"),
    start = c(100L, 220L),
    end = c(150L, 280L),
    width = c(51L, 61L),
    strand = c("+", "+"),
    intron_jxn_start = c(1L, 2L),
    intron_jxn_end = c(1L, 2L),
    intron_no = c(1L, 2L),
    annotated = c("canonical", "canonical"),
    SJ_IR = c("SJ", "SJ"),
    assembly = c("hg38", "hg38"),
    gene = c(gene_name, gene_name),
    count_placeholder = c(12L, 8L),
    pct_placeholder = c(0.60, 0.40),
    controlavg = c(0.45, 0.50),
    controlsd = c(0.10, 0.08),
    difference = c(0.15, -0.10),
    controln = c(2L, 2L),
    unique = c("1/1", ""),
    zscore = c(1.5, -1.2),
    event = c("canonical exon 1-2 splicing", "canonical exon 2-3 splicing"),
    frame_conserved = c("TRUE", "TRUE"),
    coverage_flag = c("pass", "pass"),
    note = c("", "")
  )

  setnames(
    comparison,
    c("count_placeholder", "pct_placeholder"),
    c(paste0("count_", sample_id), paste0("pct_", sample_id))
  )

  comparison
}

test_that("generateReport creates Excel and TSV outputs with expected names", {
  withr::with_tempdir({
    export_dir <- getwd()
    sample_file <- data.table(
      sampletype = "test",
      sampleID = "sample1",
      family = "family1",
      genes = "GENE1"
    )
    comparisons <- list(make_mock_comparison("sample1", "GENE1"))

    expect_no_error(
      generateReport(
        comparisons = comparisons,
        Sample_File = sample_file,
        Export = export_dir,
        mode = "default",
        prefix = "tier6_",
        debug = FALSE
      )
    )

    expect_true(file.exists(file.path(export_dir, "tier6_sample1_GENE1_combined_dt_.xlsx")))
    expect_true(file.exists(file.path(export_dir, "tier6_sample1_GENE1_combined_full.tsv")))
  })
})

test_that("generateReport preserves report column names", {
  withr::with_tempdir({
    export_dir <- getwd()
    sample_file <- data.table(
      sampletype = "test",
      sampleID = "sample1",
      family = "family1",
      genes = "GENE1"
    )
    comparison <- make_mock_comparison("sample1", "GENE1")

    generateReport(
      comparisons = list(comparison),
      Sample_File = sample_file,
      Export = export_dir,
      mode = "default",
      prefix = "tier6_",
      debug = FALSE
    )

    tsv_output <- data.table::fread(file.path(export_dir, "tier6_sample1_GENE1_combined_full.tsv"))
    xlsx_output <- data.table::as.data.table(
      openxlsx::read.xlsx(file.path(export_dir, "tier6_sample1_GENE1_combined_dt_.xlsx"))
    )

    expect_identical(names(tsv_output), names(comparison))
    expect_identical(names(xlsx_output), names(comparison))
  })
})
