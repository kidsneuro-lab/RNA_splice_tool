#' Run cortar
#'
#' `cortar()` runs the entire cortar pipeline with specified parameters,
#' and returns excel & pdf reports.
#'
#' @param file A file path, pointing to the cortar samplefile.
#' @param assembly Assembly used for alignment: either `"hg38"` or `"hg19"`
#' @param annotation  Annotation used for alignment: either `"UCSC"` or
#'     `"1000genomes"`
#' @param paired Is the RNA-seq paired-end?: `TRUE`/`FALSE`
#' @param stranded Strandedness of the RNA-seq: `0` for unstranded, `1`
#'     for forward stranded or `2` for reverse stranded
#' @param subset Does the RNA-Seq need to be subsetted to the genes
#'     of interest? (Optional, but improves speed of subsequent runs.)
#'     `TRUE`/`FALSE`
#' @param output_dir A directory path, pointing to the desired location for
#'     export of cortar results (e.g. `"output/"`)
#'
#' @export
#'
#' @examples
#' #### == COMING SOON == ####
#'

cortar <- function(file,
                   assembly,
                   annotation,
                   paired,
                   stranded,
                   subset = FALSE,
                   output_dir = "~") {
  message(paste0("Running cortar "))
  message(paste0("        file: ", file))
  message(paste0("    assembly: ", assembly))
  message(paste0("  annotation: ", annotation))
  message(paste0("      paired: ", paired))
  message(paste0("    stranded: ", stranded))
  message(paste0("      output: ", output_dir))
  message("")


  file <- data.table::fread(file)

  genes_tx <- selectGenesTranscripts(
    genes = file$genes,
    assembly = assembly,
    annotation = annotation
  )

  if (subset == T) {
    subsetBamfiles
  }

  reads <- extractCountReads(
    genes.GRanges = genes_tx[[1]],
    introns.GRanges = genes_tx[[2]][[1]],
    intron_starts.GRanges = genes_tx[[4]][[1]],
    intron_ends.GRanges = genes_tx[[4]][[2]],
    bamfiles = file$bamfile,
    sample_names = file$sampleID,
    assembly = assembly,
    paired = T,
    stranded = 2
  )

  events <- annotateQuantifyEvents(
    ids = file$sampleID,
    combined_sj = reads,
    introns.GRanges = genes_tx[[2]][[1]],
    introns_other_tx.GRanges = genes_tx[[3]],
    introns = genes_tx[[2]][[2]],
    assembly = assembly
  )

  comparisons <- (compareSplicing(
    all_splicing_events = events,
    Sample_File = file
  ))

  generateReport(
    comparisons = comparisons,
    Sample_File = file,
    Export = output_dir
  )

  message("")
  message(paste("Done! Reports saved in:", output_dir))

}
