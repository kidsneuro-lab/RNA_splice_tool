#Load Mock Data
metadata <- data.table(sampleID = "AC" ,
                       familyID = 538,
                       gene = "COL4A5",
                       transcript = "NM_001001",
                       alignment = "var1_6bp",
                       bamfile_directory = "test/directory/bamfiles/" ,
                       sjir_directory = "test/directory/sjir/" ,
                       provider = "AGRF",
                       RNAseq = "rRNA depletion",
                       assembly = "hg38",
                       age_of_biopsy = 0.2,
                       sex = "male",
                       tissue = "pbmcs",
                       chx_dmso = "chx")

db_path <- "/Volumes/PRJ-GDT/NGS_data/AGRF/GRCh38/RNASeq/KNC_RNAseq_sample_database.txt"

db <- fread("/Volumes/PRJ-GDT/NGS_data/AGRF/GRCh38/RNASeq/KNC_RNAseq_sample_database.txt")


#Integration Test -------------------------------------------------------------

test_that("integration test produces correct TSV", {
  mock_samplefile <- fread(test_path("create_new_samplefile_mock.tsv"))
  result <- create_new_samplefile(db_path, metadata, age_of_biopsy_range = 100, sex_matched = F, output = "")
  expect_equal(mock_samplefile, result) # valid_tsv_format is a function you'd need to write to check the validity of your final TSV
})

#Test generate_sample_metadata ------------------------------------------------

test_that("manually input sample metadata is correctly processed", {
  result <- generate_sample_metadata(sampleID = "AC" ,
                         familyID = 538,
                         gene = "COL4A5",
                         transcript = "NM_001001",
                         alignment = "var1_6bp",
                         bamfile_directory = "test/directory/bamfiles/" ,
                         sjir_directory = "test/directory/sjir/" ,
                         provider = "AGRF",
                         RNAseq = "rRNA depletion",
                         assembly = "hg38",
                         age_of_biopsy = 0.2,
                         sex = "male",
                         tissue = "pbmcs",
                         chx_dmso = "chx")

  expect_equal(result, metadata)
})

test_that("generate_sample_metadata flags if sampleID is a non-blank character", {})
test_that("generate_sample_metadata flags if familyID is an integer", {})
test_that("generate_sample_metadata flags if gene is a known or unknown gene", {})
test_that("generate_sample_metadata flags if alignment is in list of alignments or not", {})
test_that("generate_sample_metadata flags if bamfile_directory exists or not", {})
test_that("generate_sample_metadata flags if sjir_directory exists or not", {})
test_that("generate_sample_metadata flags if provider is in list of providers or not", {})
test_that("generate_sample_metadata flags if RNAseq is in list of RNAseq types or not", {})
test_that("generate_sample_metadata flags if assembly is in list of assemblies or not", {})
test_that("generate_sample_metadata flags if age_of_biopsy is between 0 and 140", {})
test_that("generate_sample_metadata flags if age_of_biopsy is fetal", {})
test_that("generate_sample_metadata flags if sex is in list of sexes or not", {})
test_that("generate_sample_metadata flags if tissue is in list of tissues or not", {})
test_that("generate_sample_metadata flags if chx is 'chx' or 'dmso' or not", {})


#Test generate_bam_path -------------------------------------------------------

test_that("generate_paths returns correct paths", {
  result <- generate_paths(metadata)
  expect_equal(result, c(bam = "test/directory/bamfiles/AC/AC_var1_6bp.sorted.subset.bam",
                         sj = "test/directory/sjir/AC/AC_var1_6bp/SJ.out.tab",
                         ir = "test/directory/sjir/AC/AC_var1_6bp/AC_var1_6bp.intron.average.bed.gz"))
})

test_that("generate_bam_path returns correct bamfile path", {
  result <- generate_bam_path(metadata)
  expect_equal(result, "test/directory/bamfiles/AC/AC_var1_6bp.sorted.subset.bam")
})

test_that("generate_sj_path returns correct sjfile path", {
  result <- generate_sj_path(metadata)
  expect_equal(result, "test/directory/sjir/AC/AC_var1_6bp/SJ.out.tab")
})

test_that("generate_ir_path returns correct irfile path", {
  result <- generate_ir_path(metadata)
  expect_equal(result, "test/directory/sjir/AC/AC_var1_6bp/AC_var1_6bp.intron.average.bed.gz")
})


#Test select_controls ---------------------------------------------------------

test_that("select_controls excludes rows with given gene", {
  result <- select_controls(db, metadata, age_of_biopsy_range = 100, sex_matched = F)
  expect_false(metadata$gene %in% result$gene)
})

test_that("select_controls excludes rows with gene == 'unsolved'", {
  result <- select_controls(db, metadata, age_of_biopsy_range = 100, sex_matched = F)
  expect_false("unsolved" %in% result$gene)
})

test_that("select_controls includes all samples with given familyID", {})

test_that("select_controls includes samples within age_of_biopsy_range", {
  age_of_biopsy_range <- 100
  result <- select_controls(db, metadata, age_of_biopsy_range, sex_matched = F)
  expect_true(unique(result$age_of_biopsy > metadata$age_of_biopsy - age_of_biopsy_range &
              result$age_of_biopsy < metadata$age_of_biopsy + age_of_biopsy_range))
})

test_that("select_controls includes sex matched samples when TRUE", {
  result <- select_controls(db, metadata, age_of_biopsy_range = 100, sex_matched = T)
  expect_true(unique(result$sex == metadata$sex))
})

test_that("select_controls includes only chx_dmso matched samples", {
  result <- select_controls(db, metadata, age_of_biopsy_range = 100, sex_matched = F)
  expect_true(unique(result$chx_dmso == metadata$chx_dmso))
})

test_that("select_controls includes only tissue matched samples", {
  result <- select_controls(db, metadata, age_of_biopsy_range = 100, sex_matched = F)
  expect_true(unique(result$tissue == metadata$tissue))
})

test_that("select_controls includes only RNAseq matched samples", {
  result <- select_controls(db, metadata, age_of_biopsy_range = 100, sex_matched = F)
  expect_true(unique(result$RNAseq == metadata$RNAseq))
})

test_that("select_controls includes only assembly matched samples", {
  result <- select_controls(db, metadata, age_of_biopsy_range = 100, sex_matched = F)
  expect_true(unique(result$assembly == metadata$assembly))
})


#Test generate_samplefile -----------------------------------------------------

test_that("generate_samplefile includes correct columns", {
  db_subset <- select_controls(db, metadata, age_of_biopsy_range = 100, sex_matched = F) #This is bad, don't make tests reliant on other tests
  paths <- generate_paths(metadata) #This is bad, don't make tests reliant on other tests
  result <- generate_samplefile(db_subset, metadata, paths, use.db.paths = T)
  expected_columns <- c("sampleID", "familyID", "sampletype", "genes", "bamfile", "sjfile", "irfile")
  expect_equal(colnames(result), expected_columns)
})

test_that("generate_samplefile uses db_subset file paths if selected", {})

test_that("generate_samplefile creates file paths if selected", {})

#Test export_samplefile -------------------------------------------------------
test_that("export_samplefile creates a file at the given output path", {
  mock_samplefile <- data.table(
    col1 = c(1, 2, 3),
    col2 = c("A", "B", "C")
  )
  tmp_output <- tempfile(fileext = ".tsv")
  export_samplefile(mock_samplefile, tmp_output)
  expect_true(file.exists(tmp_output))
  file.remove(tmp_output)
})

test_that("export_samplefile writes the correct content to the file", {
  mock_samplefile <- data.table(
    col1 = c(1, 2, 3),
    col2 = c("A", "B", "C")
  )
  tmp_output <- tempfile(fileext = ".tsv")
  export_samplefile(mock_samplefile, tmp_output)
  written_content <- fread(tmp_output)
  expect_equal(mock_samplefile, written_content)
  file.remove(tmp_output)
})

