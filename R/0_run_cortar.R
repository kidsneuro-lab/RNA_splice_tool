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
                   debug = FALSE,
                   ria = TRUE) {
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
  if (mode %nin% c("default", "panel", "research")) {
    stop("'", mode, "' is not an available cortar mode ('default','panel','research')")
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
  if (!is.logical(paired)) {
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
  if (isTRUE(debug)) {
    debug <- paste0(output_dir,"debug")
    message(paste0("RUNNING IN DEBUG MODE! All output will be saved to: '", debug,"'"))
    message("")
    if(!dir.exists(debug)){
      dir.create(debug)
    }
    fwrite(as.data.table(file),paste0(debug,"/","0_samplefile.tsv"), sep = "\t")
  }else{
    debug <- FALSE
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
  if (mode == "panel" || mode == "research") {
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
    genes.GRanges = genes_tx$genes,
    introns.GRanges = genes_tx$introns$granges,
    intron_starts.GRanges = genes_tx$junctions$starts,
    intron_ends.GRanges = genes_tx$junctions$ends,
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
    introns.GRanges = genes_tx$introns$granges,
    introns_other_tx.GRanges = genes_tx$introns_other_tx,
    introns = genes_tx$introns$metadata,
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
                           paired = TRUE,
                           stranded = 2,
                           subset = FALSE,
                           output_dir = "~",
                           genelist = NULL,
                           prefix = "",
                           debug = FALSE,
                           ria = TRUE){
  batches_in <- sapply(list.files(folder, pattern = pattern),function(x){paste0(folder,"/",x)})
  batches_out <- sapply(list.files(folder, pattern = pattern),function(x){paste0(output_dir,"/",strsplit(x,"\\.")[[1]][1])})
  for(batch in seq_along(batches_in)){
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

#' Create BED Coordinates for Genes from RefSeq Exon Annotations
#'
#' `subsetBamfiles()` looks up each requested gene in the bundled RefSeq exon annotation
#' (`refseq_introns_exons_hg19` or `refseq_introns_exons_hg38`), computes the gene span
#' (min exon start to max exon end), applies an optional flanking `overhang`, and prints
#' BED-style coordinates to the console.
#'
#' The output is intended for tools that accept BED intervals (tab-delimited:
#' `chrom`, `start`, `end`). Note that BED uses a **0-based, half-open** start coordinate,
#' so this function subtracts 1 from the computed start after applying `overhang`.
#'
#' @param genes Character vector of gene identifiers to look up in RefSeq exon annotations.
#'   By default these are RefSeq `gene_name` values; if `use_ncbi_gene_id = TRUE`, these
#'   are matched to RefSeq `gene_id` values instead.
#' @param hg Genome assembly to use. Must be `"19"` (hg19) or `"38"` (hg38).
#' @param overhang Non-negative numeric scalar. Number of base pairs to extend each gene
#'   span upstream (subtract from start) and downstream (add to end). Default is `1000`.
#' @param use_ncbi_gene_id Logical; if `TRUE`, match `genes` against RefSeq `gene_id`
#'   (NCBI Gene ID field in the annotation). If `FALSE` (default), match against `gene_name`.
#'
#' @return Invisibly returns `NULL`. The function prints one BED interval per gene to the
#'   console via `cat()`, with fields separated by tabs:
#'   \code{chr<chrom>\t<start0>\t<end>}
#'
#' @details
#' The gene span is obtained for the specified gene(s) with specified overhang subtracted
#' and added at each gene's start and end coordinates respectively. The function can either
#' determine the gene spans using NCBI gene id or HGNC names
#'
#' Chromosomes are formatted with a `chr` prefix. `MT` is converted to `chrM` (i.e.,
#' `chrMT` becomes `chrM`) for compatibility with common BAM/BED conventions.
#'
#' The function stops if:
#' \itemize{
#'   \item none of the requested genes are found in the selected RefSeq annotation, or
#'   \item only a subset are found (it reports which identifiers were missing).
#' }
#'
#' @examples
#' \dontrun{
#' # Using gene symbols (gene_name) on hg38 with no overhang
#' subsetBamfiles(genes = c("EMD", "DMD"), hg = "38", overhang = 0)
#'
#' # Using NCBI gene IDs (gene_id) on hg38
#' subsetBamfiles(genes = c("1234", "5678"), hg = "38", overhang = 500,
#'                use_ncbi_gene_id = TRUE)
#' }
#'
#' @import data.table
#' @export
subsetBamfiles <- function(genes, hg = c('19','38'), overhang = 1000, use_ncbi_gene_id = FALSE){
  stopifnot(length(genes) > 0L)

  hg <- match.arg(hg, c("19", "38"))

  if (!is.numeric(overhang) || length(overhang) != 1L || is.na(overhang) || overhang < 0) {
    stop("`overhang` must be a single non-negative numeric value.")
  }
  overhang <- as.integer(overhang)

  if (any(duplicated(genes))) {
    stop("Duplicate input genes. Please check.")
  }

  # Select correct gene annotation for chosen assembly. Avoid modifying the global table
  Refseq_Genes <- data.table::copy(
    if (hg == "38") refseq_introns_exons_hg38 else refseq_introns_exons_hg19
  )

  Refseq_Genes[, requested_identifier := if (isTRUE(use_ncbi_gene_id)) gene_id else gene_name]

  refseq_genes <- Refseq_Genes[region_type == 'exon',
                              .(start = min(region_start),
                                end = max(region_end)), by = .(gene_id,
                                                               gene_name,
                                                               requested_identifier,
                                                               chrom)]

  # Create a data.table based on input genes
  requested_identifiers <- data.table::data.table(
    requested_identifier = genes,
    input_order_index = seq_along(genes)
  )

  # Obtain the refseq genes entry for each input gene
  matched_gene_spans <- refseq_genes[requested_identifiers, on = .(requested_identifier)]

  if (nrow(matched_gene_spans) == 0) {
    stop("No NCBI genes found corresponding to genes of interest")
  }

  # No. of genes where NCBI gene entry not found
  genes_with_no_NCBI_entry <- matched_gene_spans[is.na(gene_id)]

  if (nrow(genes_with_no_NCBI_entry) > 0) {
    stop(sprintf("Number of NCBI genes does not match genes of interest. NCBI genes not found for: %s", paste(genes_with_no_NCBI_entry[,requested_identifier], collapse = ',')))
  }

  # Ensure output is in the same order as input
  data.table::setorder(matched_gene_spans, input_order_index)

  matched_gene_spans[,chrom := ifelse(startsWith(chrom, "chr"), chrom, paste0("chr", chrom))]
  matched_gene_spans[, chrom := ifelse(chrom == "chrMT", "chrM", chrom)]

  # Apply overhang, clamping start to 0 to avoid negative coordinates
  matched_gene_spans[, start := pmax(0L, start - 1 - as.integer(overhang))]
  matched_gene_spans[, end   := end + as.integer(overhang)]

  # BED is 0-based start, so subtract 1 from start
  formatted_gene_coords <- sprintf(
    "%s\t%d\t%d",
    matched_gene_spans$chrom,
    matched_gene_spans$start,
    matched_gene_spans$end
  )

  cat(paste(formatted_gene_coords, collapse = "\n"))
  invisible(formatted_gene_coords)
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

  dir.create(dest_path, showWarnings = TRUE, recursive = TRUE)
  dir.create(paste0(dest_path,"/output"), showWarnings = TRUE, recursive = TRUE)
  file.copy(from = extdata_path, to = dest_path, recursive = TRUE)

  test_samplefile.tsv <- fread(paste0(dest_path,"/extdata/test_samplefile.tsv"))

  test_samplefile.tsv$bamfile <- paste0(dest_path,"/extdata/",test_samplefile.tsv$bamfile)

  fwrite(test_samplefile.tsv, paste0(dest_path,"/extdata/test_samplefile.tsv"), sep = "\t")

  cortar::cortar(paste0(dest_path,"/extdata/test_samplefile.tsv"),
         input_type = "bamfile",
         output_dir = paste0(dest_path,"/output"))

}
