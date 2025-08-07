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
                   ria = T) {
                   # reads in absentia - count multi-exon skipping as an event for a skipped intron

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

#' Subset BAM Files Based on Gene Annotations
#'
#' The `subsetBamfiles` function retrieves and formats genomic ranges for specified genes
#' based on the chosen genome assembly (hg19 or hg38). It supports gene names, Ensembl gene IDs,
#' and transcript IDs (RefSeq). The function adjusts the genomic coordinates by a specified
#' overhang and outputs the formatted ranges as strings.
#'
#' @param genes A character vector of gene identifiers. This can include:
#'   \itemize{
#'     \item Gene names (e.g., `"EMD"`, `"DMD"`)
#'     \item Ensembl gene IDs (e.g., `"ENSG00000231514"`)
#'     \item Transcript IDs (e.g., `"NM_XXXXXX"`)
#'   }
#' @param hg An integer specifying the genome assembly version. Supported values are `19` (hg19) and `38` (hg38).
#' @param overhang A numeric value indicating the number of base pairs to subtract from the start position
#'   and add to the end position of each genomic range. Default is `1000`.
#'
#' @return The function outputs formatted genomic range strings to the console using `cat()`.
#'   Each range is formatted as:
#'   \code{  - "''chr<chrom>:<start>-<end>''"}
#'
#' @examples
#' \dontrun{
#' # Example 1: Obtaining genomic ranges for multiple genes based on gene names
#' subsetBamfiles(genes = c('EMD', 'DMD'), hg = 38, overhang = 0)
#' # Output:
#' #   - "''chrX:154379567-154380881''"
#' #   - "''chrX:31121931-33211281''"
#'
#' # Example 2: Obtaining genomic ranges for a single gene based on gene name
#' subsetBamfiles(genes = c('MPP5'), hg = 38, overhang = 0)
#' # Output:
#' #   - "''chr14:67240713-67336061''"
#'
#' # Example 3: Obtaining genomic ranges based on Ensembl gene ID
#' subsetBamfiles(genes = c('ENSG00000231514'), hg = 38, overhang = 0)
#' # Output:
#' #   - "''chrY:26626520-26627159''"
#'
#' # Example 4: Throwing an error if gene is not found
#' subsetBamfiles(genes = c('NOT_FOUND'), hg = 38, overhang = 0)
#' # Error:
#' # [NOT_FOUND] not found in Ensembl or Refseq. Please check input values.
#' }
#'
#' @details
#' The function processes each gene identifier to determine its corresponding genomic coordinates.
#' It handles different types of identifiers by checking against RefSeq and Ensembl gene annotations.
#' If a gene cannot be located in either annotation set, the function will terminate and throw an error.
#'
#' The `overhang` parameter allows users to extend or reduce the genomic range boundaries,
#' which can be useful for including additional upstream or downstream regions.
#'
#' @seealso
#' \code{\link[testthat]{test_that}}, \code{\link[data.table]{data.table}}
#'
#' @import data.table
#' @export
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

  # First check if the gene can be located in

  gene_coordinates <- list()

  for (gene_counter in seq(1, length(genes))) {
    gene <- genes[gene_counter]

    if (grepl("^ENSG", gene)) {
      result <- Ensembl_Genes[`Gene stable ID` == gene,
                              .(chrom = `Chromosome/scaffold name`,
                                start = `Gene start (bp)`,
                                end = `Gene end (bp)`)]
      gene_coordinates[[gene]] <- as.list(result)

    } else if (grepl("^NM_", gene) || grepl("^ENST", gene)) {
      tx <- gsub(pattern = "\\.\\d+$", "", gene)
      result <- Refseq_Genes[tx_id == tx & region_type == 'intron' & (region_no==1 | last_region==1)][order(ifelse(strand == '+', last_region, -last_region))]
      result_list <- as.list(result[,.(start = min(region_start),
                                       end = max(region_end)), by = .(chrom)])
      gene_coordinates[[gene]] <- result_list

    } else if (gene %in% Refseq_Genes$gene_name) {
      result <- Refseq_Genes[gene_name == gene & canonical == 1 & region_type == 'intron' & (region_no==1 | last_region==1)][order(ifelse(strand == '+', last_region, -last_region))]
      result_list <- as.list(result[,.(start = min(region_start),
                                       end = max(region_end)), by = .(chrom)])

      gene_coordinates[[gene]] <- result_list

    } else if (gene %in% Ensembl_Genes$`Gene name`) {
      result <- Ensembl_Genes[`Gene name` == gene,
                              .(chrom = `Chromosome/scaffold name`,
                                start = `Gene start (bp)`,
                                end = `Gene end (bp)`)]
      gene_coordinates[[gene]] <- as.list(result)

    } else (
      stop(paste0("[", gene, "] not found in Ensembl or Refseq. Please check input values."))
    )
  }

  updated_gene_coordinates <- lapply(gene_coordinates, function(sublist) {
    # Append 'chr' to the chrom value. Ensure MT chromosome formatting is correct
    sublist$chrom <- ifelse(paste0("chr", sublist$chrom) == 'chrMT', 'chrM', paste0("chr", sublist$chrom))

    # Subtract overhang from the start position
    sublist$start <- sublist$start - overhang

    # Add overhang to the end position
    sublist$end <- sublist$end + overhang

    # Return the modified sublist
    return(sublist)
  })

  formatted_gene_coords <- sapply(updated_gene_coordinates, function(sublist) {
    # Create the formatted string with double single quotes
    formatted_string <- paste0("  - ", '"', "''", sublist$chrom, ":", sublist$start, "-", sublist$end, "''", '"')
    return(formatted_string)
  })

  formatted_gene_coords <- unname(formatted_gene_coords)

  cat(paste(formatted_gene_coords, collapse = "\n"))
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
