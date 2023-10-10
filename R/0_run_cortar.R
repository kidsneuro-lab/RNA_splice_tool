#' Run cortar
#'
#' `cortar()` runs the entire cortar pipeline with specified parameters,
#' and returns excel & pdf reports.
#'
#' @param file A file path, pointing to the cortar samplefile - see readme for
#'     more information
#' @param mode One of `"default"` (default),`"panel"`, or `"research"` - see
#'     readme for more information
#' @param assembly Assembly used for alignment: either `"hg38"` (default) or
#'     `"hg19"`
#' @param annotation  Annotation used for alignment: either `"UCSC"` or
#'     `"1000genomes"`
#' @param paired Is the RNA-seq paired-end?: `TRUE`/`FALSE`
#' @param stranded Strandedness of the RNA-seq: `0` for unstranded, `1`
#'     for forward stranded or `2` for reverse stranded
#' @param subset Does the RNA-Seq need to be subsetted to the genes
#'     of interest? `TRUE`/`FALSE` (Optional but significantly improves speed
#'     of subsequent analyses - not currently available)
#' @param output_dir A directory path, pointing to the desired location for
#'     export of cortar results (e.g. `"output/"`)
#' @param genelist A character vector with genes/RefSeq transcripts of interest
#'     (Only for panel or research modes; default = NULL)
#' @param prefix A character vector to be appended to the beginning of the
#'     output file names
#'
#' @export
#'
#' @examples
#' #### == COMING SOON == ####
#'
cortar <- function(file,
                   mode = "default",
                   assembly = "hg38",
                   annotation = "UCSC",
                   input_type = "bamfile",
                   paired = TRUE,
                   stranded = 2,
                   subset = NULL,
                   output_dir = "~",
                   genelist = NULL,
                   prefix = "",
                   debug = F,
                   ria = F) {

  # Error catching
  # file
  if (!file.exists(file)) {
    stop(
      "File '", file, "' does not exist or is non-readable. getwd()=='",
      getwd(), "'"
    )
  }

  # mode
  if (mode != "default") {
    if (mode == "panel" | mode == "research") {
    } else {
      stop("'", mode, "' is not an available cortar mode ('default','panel','research')")
    }
  }

  # assembly
  if (assembly %nin% c("hg19", "hg38")) {
    stop("Assembly '", assembly, "' is not an available assembly ('hg38','hg19')")
  }

  # annotation
  if (annotation %nin% c("1000genomes", "UCSC")) {
    stop("Annotation '", annotation, "' is not an available annotation ('1000genomes','UCSC')")
  }

  # paired
  if (class(paired) != "logical") {
    stop("Paired must be a logical: TRUE, FALSE.
         Supplied:", paired, "
         Note: You may need to remove quotation marks.")
  }

  # output_dir
  if (!dir.exists(output_dir)) {
    message("Directory '", output_dir, "' does not exist.
    Would you like to create this directory?
    \t 1. Yes
    \t 2. No")
    selection <- readline(prompt = "Selection: ")
    if (selection == "1") {
      dir.create(output_dir)
      message("Directory '", output_dir, "' created.")
      message("")
    } else {
      stop(
        "Directory '", output_dir, "' does not exist. getwd()=='",
        getwd(), "'"
      )
    }
  }

  # Initialisation messages
  message(paste0("Running cortar "))
  message(paste0("        file: ", file))
  message(paste0("        mode: ", mode))
  message(paste0("    assembly: ", assembly))
  message(paste0("  annotation: ", annotation))
  message(paste0("      paired: ", paired))
  message(paste0("    stranded: ", stranded))
  message(paste0("      output: ", output_dir))
  message("")

  # Read in cortar samplefile
  file <- data.table::fread(file)
  for(bamfile in file$bamfile){
    if(!file.exists(bamfile)){
      stop(
        "File does not exist or is non-readable. path = '",
        bamfile, "'"
      )
    }
  }

  # Initialise debug directory
  if(debug == T){
    debug <- paste0(output_dir,"debug")
    message(paste0("RUNNING IN DEBUG MODE! All output will be saved to: '", debug,"'"))
    message("")
    if(!dir.exists(debug)){
      dir.create(debug)
    }
    fwrite(as.data.table(file),paste0(debug,"/","0_samplefile.tsv"), sep = "\t")
  }else{
    debug <- ""
  }

  # Select data of interest
  if (input_type == "bamfile"){
    file$sjfile <- ""
    file$irfile <- ""
  } else if (input_type == "sj"){
    file$bamfile <- ""
  }

  # Select genes and transcripts of interest
  # A genelist must be provided for panel or research mode
  if (mode == "panel" | mode == "research") {
    genes_tx <- selectGenesTranscripts(
      genes = genelist,
      assembly = assembly,
      annotation = annotation,
      debug = debug
    )
  }

  # Genes specified in the cortar samplefile are used in default
  if (mode == "default") {
    genes_tx <- selectGenesTranscripts(
      genes = file$genes,
      assembly = assembly,
      annotation = annotation,
      debug = debug
    )
  }

  # Reads for the genes and transcripts of interest are extracted and counted
  reads <- extractCountReads(
    genes.GRanges = genes_tx[[1]],
    introns.GRanges = genes_tx[[2]][[1]],
    intron_starts.GRanges = genes_tx[[4]][[1]],
    intron_ends.GRanges = genes_tx[[4]][[2]],
    bamfiles = file$bamfile,
    sjfiles = file$sjfile,
    irfiles = file$irfile,
    sample_names = file$sampleID,
    assembly = assembly,
    annotation = annotation,
    paired = paired,
    stranded = stranded,
    input = input_type,
    debug = debug
  )

  # Events supported by extracted reads are annotated and quantified
  events <- annotateQuantifyEvents(
    ids = file$sampleID,
    combined_sj = reads,
    introns.GRanges = genes_tx[[2]][[1]],
    introns_other_tx.GRanges = genes_tx[[3]],
    introns = genes_tx[[2]][[2]],
    assembly = assembly,
    debug = debug,
    ria = ria
  )

  # Comparisons of events between test samples and controls are performed
  comparisons <- compareSplicing(
    all_splicing_events = events,
    Sample_File = file,
    mode = mode,
    debug = debug
  )

  # Reports in excel and pdf formats are generated for comparisons
  generateReport(
    comparisons = comparisons,
    Sample_File = file,
    Export = output_dir,
    mode = mode,
    prefix = prefix,
    debug = debug
  )

  # Termination message
  message("")
  message(paste("Done! Reports saved in:", output_dir))
  message("")
}

#' Run multiple instances of cortar
#'
#' `cortar_batch()` runs the entire cortar pipeline with specified parameters,
#' and returns excel & pdf reports for multiple samplefiles.
#'
#' @param folder A character vector of the full file path, pointing to a folder
#'     containing one or more cortar samplefiles
#' @param pattern An optional regular expression. Only file names which match
#'     the regular expression will be returned (default: `"*.tsv"`).
#' @inheritParams cortar
#'
#' @export
#'
#' @examples
#' #### == COMING SOON == ####
#'
cortar_batch <- function(folder,
                           pattern = "*.tsv",
                           mode = "default",
                           assembly = "hg38",
                           annotation = "UCSC",
                           input_type = "sj",
                           paired = T,
                           stranded = 2,
                           subset = F,
                           output_dir = "~",
                           genelist = NULL,
                           prefix = "",
                           debug = debug,
                           ria = T){
  batches_in <- sapply(list.files(folder, pattern = pattern),function(x){paste0(folder,"/",x)})
  batches_out <- sapply(list.files(folder, pattern = pattern),function(x){paste0(output_dir,"/",strsplit(x,"\\.")[[1]][1])})
  for(batch in seq(1,length(batches_in))){
    if(!dir.exists(batches_out[batch])){
      message("'",batches_out[batch],"' created.")
      dir.create(batches_out[batch])
    }
    cortar(file = batches_in[batch],
           mode = mode,
           assembly = assembly,
           annotation = annotation,
           input_type = input_type,
           paired = paired,
           stranded = stranded,
           subset = subset,
           output_dir = batches_out[batch],
           genelist = genelist,
           prefix = prefix,
           debug = debug,
           ria = ria)
  }
}

#' Extract gene coordinates for RNA-seq subsetting
#'
#' Returns the gene coordinates for a given entrez_gene_symbol or
#' RefSeq/Ensembl transcript to be used for RNA-seq subsetting
#'
#' @param genes A character vector with the entrez_gene_symbols,
#' RefSeq transcript ids, or Ensembl gene or transcript ids for the genes for
#' analysis.
#' @param hg Either 19 or 38
#' @param overhang Number of nucleotides flanking the gene to be included
#' (default: 1000nt)
#'
#' @return A character vector.
#' @export
#'
#' @examples
#' == COMING SOON==
#'

subsetBamfiles <- function(genes, hg, overhang = 1000){

  # Read in cortar samplefile
  # file <- data.table::fread(file)

  # Select correct gene annotation for chosen assembly
  if (hg == 38) {
    Refseq_Genes <- refseq_introns_exons_hg38
    Ensembl_Genes <- ensembl_allgenes_chr1_Y_hg38
  } else if (hg == 19) {
    Refseq_Genes <- refseq_introns_exons_hg19
    Ensembl_Genes <- ensembl_allgenes_chr1_Y_hg19
  }

  genes_tx <- data.table("gene_name" = character(), "tx" = character())

  for (gene in seq(1, length(genes))) {
    if (length(grep("NM_[0-9]+\\.[0-9]+", genes[gene])) > 0) {
      tx <- stringr::str_extract_all(genes[gene], "NM_[0-9]+\\.[0-9]+")[[1]][1]
      gene_name <- unique(Refseq_Genes[tx_version_id == tx & region_type == c("intron"), "gene_name"])[[1]]
    } else {
      gene_name <- genes[gene]
      tx <- unique(Refseq_Genes[gene_name == genes[gene] & canonical == 1, "tx_version_id"])[[1]]
    }
    genes_tx <- rbind(genes_tx, data.table("gene_name" = gene_name, "tx" = tx))
  }

  subsetgenes <- Ensembl_Genes[`Gene name` %in% genes_tx$gene_name]

  # Format gene input
  forsubset <- paste0("  - \"'","\'","chr",subsetgenes$`Chromosome/scaffold name`,":",subsetgenes$`Gene start (bp)`-overhang,"-",
                      subsetgenes$`Gene end (bp)`+overhang,"\'","'\"")
  forsubset <- paste(forsubset,collapse="\n")

  cat(forsubset)
}


#' Run the cortar test
#'
#' This function runs a test for the `cortar` package by copying the necessary files from the
#' package's extdata directory to a specified test directory. After copying, it processes a sample
#' TSV file and then calls the `cortar` function.
#'
#' @param test_path Character string indicating the path to the directory where the test will run.
#'   By default, it uses the current working directory (`getwd()`).
#'
#' @return None. The function will perform the tests and modify files, but does not return any values.
#' @export
#'
#' @examples
#' \dontrun{
#' run_cortar_test("/path/to/your/test/directory")
#' }

run_cortar_test <- function(test_path = getwd()){
  extdata_path <- system.file("extdata", package="cortar")

  # Check if the folder exists
  if (!dir.exists(extdata_path)) {
    stop("The specified folder does not exist in the package.")
  }

  dest_path <- paste0(test_path,"/cortar_test")

  dir.create(dest_path, showWarnings = T, recursive = T)
  dir.create(paste0(dest_path,"/output"), showWarnings = T, recursive = T)
  file.copy(from = extdata_path, to = dest_path, recursive = T)

  test_samplefile.tsv <- fread(paste0(dest_path,"/extdata/test_samplefile.tsv"))

  test_samplefile.tsv$bamfile <- paste0(dest_path,"/extdata/",test_samplefile.tsv$bamfile)

  fwrite(test_samplefile.tsv, paste0(dest_path,"/extdata/test_samplefile.tsv"), sep = "\t")

  cortar::cortar(paste0(dest_path,"/extdata/test_samplefile.tsv"),
         input_type = "bamfile",
         output_dir = paste0(dest_path,"/output"))

}
