#Integration Test

generate_sample_metadata <- function(sampleID,
                                     familyID,
                                     gene,
                                     transcript,
                                     alignment,
                                     bamfile_directory,
                                     sjir_directory,
                                     provider,
                                     RNAseq,
                                     assembly,
                                     age_of_biopsy,
                                     sex,
                                     tissue,
                                     chx_dmso){

  metadata <- data.table(sampleID = sampleID ,
                         familyID = familyID,
                         gene = gene,
                         transcript = transcript,
                         alignment = alignment,
                         bamfile_directory = bamfile_directory ,
                         sjir_directory = sjir_directory ,
                         provider = provider,
                         RNAseq = RNAseq,
                         assembly = assembly,
                         age_of_biopsy = age_of_biopsy,
                         sex = sex,
                         tissue = tissue,
                         chx_dmso = chx_dmso)

}

create_new_samplefile <- function(db_path, metadata, age_of_biopsy_range, sex_matched = F, output = ""){

  db <- fread(db_path)

  paths <- generate_paths(metadata)

  db_subset <- select_controls(db, metadata, age_of_biopsy_range = 100, sex_matched = F)

  samplefile <- generate_samplefile(db_subset, metadata, paths)

  if(stringr::str_length(output) > 0){
    export_samplefile(samplefile, output)
    message(paste0("Exported samplefile to: '", output,"'"))
  }

  samplefile

}

generate_paths <- function(metadata){
  bam_path <- generate_bam_path(metadata)
  sj_path <- generate_sj_path(metadata)
  ir_path <- generate_ir_path(metadata)

  c(bam = bam_path, sj = sj_path, ir = ir_path)
}

select_controls <- function(db, metadata, age_of_biopsy_range, sex_matched){
  db_subset_incl_excl <- db[(age_of_biopsy > as.numeric(metadata$age_of_biopsy) - age_of_biopsy_range) &
                            (age_of_biopsy < as.numeric(metadata$age_of_biopsy) + age_of_biopsy_range) &
                             assembly == metadata$assembly &
                             RNAseq == metadata$RNAseq &
                             tissue == metadata$tissue &
                             chx_dmso == metadata$chx_dmso &
                             genes != metadata$gene &
                             genes != "unsolved"]

  if(sex_matched == T){
    db_subset_incl_excl <- db_subset_incl_excl[sex == metadata$sex]
  }

  db_family <- db[familyID == metadata$familyID & chx_dmso == metadata$chx_dmso]

  db_subset <- unique(rbindlist(list(db_subset_incl_excl, db_family)))

  db_subset
}

generate_samplefile <- function(db_subset, metadata, paths, use.db.paths = T){

  test_sample <- data.table(sampleID = metadata$sampleID,
                            familyID = metadata$familyID,
                            sampletype = "test",
                            genes = metadata$gene,
                            bamfile = paths[["bam"]],
                            sjfile = paths[["sj"]],
                            irfile = paths[["ir"]])

  control_subset <- db_subset[,sampletype := ""]

  if(use.db.paths == T){
    control_subset <- db_subset[,.(sampleID,
                                   familyID,
                                   sampletype,
                                   genes,
                                   bamfile = bamfile_6bp_csj,
                                   sjfile = sjfile,
                                   irfile = irfile)]
  }

  samplefile <- unique(rbindlist(list(test_sample, control_subset)))

  samplefile

}

export_samplefile <- function(samplefile, output){

  fwrite(samplefile, output, sep = "\t")

}


generate_bam_path <- function(metadata){
  paste0(metadata$bamfile_directory,metadata$sampleID,"/",metadata$sampleID,"_",metadata$alignment,".sorted.subset.bam")
}

generate_sj_path <- function(metadata){
  paste0(metadata$sjir_directory,metadata$sampleID,"/",metadata$sampleID,"_",metadata$alignment,"/SJ.out.tab")
}

generate_ir_path <- function(metadata){
  paste0(metadata$sjir_directory,metadata$sampleID,"/",metadata$sampleID,"_",metadata$alignment,"/",metadata$sampleID,"_",metadata$alignment,".intron.average.bed.gz")
}
