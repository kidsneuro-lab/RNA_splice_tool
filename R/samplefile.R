# Load data on startup
pool <- pool::dbPool(
  drv = RSQLite::SQLite(),
  dbname = "/Volumes/research-data/PRJ-GDT/NGS_data/KNC_RNAseq_samples.sqlite3"
)

# example: get_all_sample_ids(pool, "AGRF", "pbmcs","dmso")
get_all_sample_ids <- function(pool, provider = "", tissue = "", treatment = ""){
  basecall <- "SELECT sampleID FROM samples"
  provider_call <- ifelse(provider != "", paste0('provider="',provider,'"'), "")
  tissue_call <- ifelse(tissue != "", paste0('tissue="',tissue,'"'), "")
  treatment_call <- ifelse(treatment != "", paste0('chx_dmso="',treatment,'"'), "")
  option_calls <- c(provider_call, tissue_call, treatment_call)
  if(length(option_calls[option_calls != ""])){
    where_call <- paste0(c(option_calls[option_calls != ""]), collapse = " and ")
  }else{
    where_call <- ""
  }
  if(where_call != ""){
    where_call <- paste0(" WHERE ", where_call)
  }
  call <- paste0(basecall,where_call)
  ids <- DBI::dbGetQuery(pool, call)
  ids
}

# example: get_all_tissue_treatments(pool, "AGRF")
get_all_tissue_treatments <- function(pool, provider){
  call <- paste0('SELECT tissue, chx_dmso FROM samples WHERE provider="',provider,'"')
  tissues <- DBI::dbGetQuery(pool, call)
  unique(tissues)
}


get_sample_information <- function(pool, sample){
  DBI::dbGetQuery(pool, paste0("SELECT sampleID, familyID, genes, transcript, csj, sampleID, familyID, sampletype, age_at_biopsy, sex, tissue, chx_dmso, genes, transcript, provider FROM samples WHERE sampleID = '",sample,"'"))
}

get_control_samples <- function(pool, tissue, treatment, provider, sex, age, genes, transcript, ageMatch = F, sexMatch = F, includeUnsolved = F){
  basecall <- "SELECT sampleID, familyID, genes, transcript, csj, sampleID, familyID, sampletype, age_at_biopsy, sex, tissue, chx_dmso, genes, transcript, provider FROM samples WHERE "
  provider_call <- paste0('provider="',provider,'" and ')
  tissue_call <- paste0('tissue="',tissue,'" and ')
  treatment_call <- paste0('chx_dmso="',treatment,'" and ')
  sex_call <- paste0('sex="',sex,'" and ')
  genes_call <- paste0('genes!="',genes,'" and ')
  transcript_call <- paste0('transcript!="',transcript,'"')
  unsolved_call <- paste0('genes!="unsolved" and ')
  additional_call <- ""
  if(sexMatch == T){
    additional_call <- paste0(sex_call)
  }
  if(includeUnsolved == F){
    additional_call <- paste0(additional_call, unsolved_call)
  }
  where_call <- paste0(additional_call, provider_call, tissue_call, treatment_call, genes_call, transcript_call)
  call <- paste0(basecall, where_call)
  DBI::dbGetQuery(pool, call)
}

# make a samplefile
# samplefile <- make_samplefile(pool, "LeMa_DMD")
make_samplefile <- function(pool, sample, sexMatch = F, ageMatch = F, includeUnsolved = F){
  sample_data <- get_sample_information(pool, sample)
  sample_data$test <- "test"
  tissue <- sample_data$tissue
  treatment <- sample_data$chx_dmso
  provider <- sample_data$provider
  sex <- sample_data$sex
  age <- sample_data$age_at_biopsy
  genes <- sample_data$genes
  transcript <- sample_data$transcript
  control_data <- get_control_samples(pool, tissue, treatment, provider, sex, age, genes, transcript, ageMatch, sexMatch, includeUnsolved)
  control_data$test <- "control"
  samplefile <- rbind(sample_data[,c(1,16,2:5)], control_data[,c(1,16,2:5)])
  rundata <- rbind(sample_data[,6:15], control_data[,6:15])
  list('samplefile' = samplefile, 'rundata' = rundata)
}

# '/Volumes/research-data/PRJ-GDT/NGS_data/Rhett/alignments/results'
# example:
get_bamfiles <- function(samplefile, folder, append = T, run = "_var1_6bp", suffix = ".sorted.subset.bam", csj = T, force = F){
  if(csj == F){
    bamfiles <- apply(samplefile$samplefile, 1, function(x) paste0(folder,"/",x['sampleID'],"/",x['sampleID'],run,"/",x['sampleID'],run,suffix))
    #bamfiles <- apply(samplefile$samplefile, 1, function(x) paste0(folder,"/",x['sampleID'],run,suffix))
  } else {
    samplefile$samplefile$csj[samplefile$samplefile$csj == "NA"] <- ""
    bamfiles <- apply(samplefile$samplefile, 1, function(x) paste0(folder,"/",x['sampleID'],"/",x['sampleID'],run,x['csj'],"/",x['sampleID'],run,x['csj'],suffix))
  }
  missing_files <- sapply(bamfiles, file.exists) == FALSE

  if(force == T){
    bamfiles <- bamfiles[missing_files == FALSE]
    samplefile$samplefile <- samplefile$samplefile[missing_files == FALSE,]
    samplefile$rundata <- samplefile$rundata[missing_files == FALSE,]
  }

  if(append == T){
    samplefile$samplefile <- samplefile$samplefile[,1:5]
    samplefile$samplefile$bamfiles <- bamfiles
    bamfiles <- samplefile
  }

  if(sum(missing_files) > 0 & force == F){
    message("To ignore this error, review the missing_files listed, then resubmit with the flag 'drop_missing_bamfiles = TRUE'")
    stop(paste0("Missing file: ", c(names(missing_files[missing_files == TRUE])), collapse = "\n  "))
  }

  list('bamfiles' = bamfiles, 'missing' = missing_files)
}



finish_samplefile <- function(bamfiles, tests, coverage = 1){
  test_index <- bamfiles$bamfiles$samplefile$sampleID %in% tests
  bamfiles$bamfiles$samplefile$coverage <- coverage
  bamfiles$bamfiles$samplefile$test[test_index] <- 'test'
  bamfiles$bamfiles$samplefile$genes[test_index == FALSE] <- ''
  bamfiles$bamfiles$samplefile$transcript[test_index == FALSE] <- ''
  names(bamfiles$bamfiles$samplefile) <- c('sampleID','sampletype','familyID','genes','transcript','bamfiles', 'coverage')
  bamfiles
}

create_samplefile <- function(pool, sample, folder, output_dir, csj = F, sexMatch = F, ageMatch = F, includeUnsolved = F, drop_missing_bamfiles = F){
  samplefile <- make_samplefile(pool, sample, sexMatch = F, ageMatch = F, includeUnsolved = F)
  bamfiles <- get_bamfiles(samplefile, folder, csj = csj, force = drop_missing_bamfiles)
  final <- finish_samplefile(bamfiles, sample)

  data.table::fwrite(final$bamfiles$samplefile,paste0(output_dir,sample,"_samplefile.tsv"), sep = "\t")
  data.table::fwrite(final$bamfiles$rundata,paste0(output_dir,sample,"_rundata.tsv"), sep = "\t")

  final
}

# get_all_sample_ids(pool)

# create_samplefile(pool, "609_NF1_P", "/Volumes/research-data/PRJ-GDT/NGS_data/Rhett/alignments/results", "../test/")

# create_samplefile(
#   pool,
#   "422_POMT2_P_D",
#   "/Volumes/research-data/PRJ-GDT/NGS_data/AGRF/GRCh38/Alignments",
#   "/Volumes/research-data/PRJ-GDT/NGS_data/Rhett/NMD_ASO_Targets"
# )
