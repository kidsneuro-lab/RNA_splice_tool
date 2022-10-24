#' Split a string
#'
#' @param string A character vector with, at most, one element.
#' @inheritParams stringr::str_split
#'
#' @return A character vector.
#' @export
#'
#' @examples
#' x <- "alfa,bravo,charlie,delta"
#' str_split_one(x, pattern = ",")
#' str_split_one(x, pattern = ",", n = 2)
#'
#' y <- "192.168.0.1"
#' str_split_one(y, pattern = stringr::fixed("."))
#'
# INPUT = Vector of gene names and/or transcripts
# OUTPUT = Subsetted Refseq table and GRanges object
#'
cortar <- function(file,
                   assembly,
                   annotation,
                   paired,
                   stranded,
                   subset = F,
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
