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
#'
#'
#--Extract all junctions and intron retention (non-split >4bp intron/>4bp exon) reads----
#'
extractCountReads <- function(genes.GRanges,
                              introns.GRanges,
                              intron_starts.GRanges,
                              intron_ends.GRanges,
                              bamfiles,
                              sample_names,
                              assembly,
                              annotation,
                              paired,
                              stranded) {

  # inputs:
  #   - genes.GRanges
  #   - sample_file (ID)
  #   - intron_starts.GRanges
  #   - intron_ends.GRanges
  #   - introns.GRanges
  #   - genome

  message("Extracting and counting reads...")

  if (assembly == "hg19") {
    if (annotation == "UCSC") {
      Genome_Assembly <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    } else if (annotation == "1000genomes") {
      Genome_Assembly <- BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
    }
  } else if (assembly == "hg38") {
    Genome_Assembly <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  }

  sj <- list()
  ir <- list()

  param <- Rsamtools::ScanBamParam(
    which = genes.GRanges, flag = Rsamtools::scanBamFlag(
      isDuplicate = FALSE,
      isSecondaryAlignment = FALSE,
      isPaired = paired
    )
  )

  for (sample_number in 1:length(sample_names)) {
    sample_name <- sample_names[sample_number]
    message("\t", sample_name)

    if (paired == F) {
      alignment <- GenomicAlignments::readGAlignments(
        file = bamfiles[sample_number],
        param = param
      )
    } else if (paired == T) {
      alignment <- GenomicAlignments::readGAlignmentPairs(
        file = bamfiles[sample_number],
        param = param,
        strandMode = stranded
      )
    }

    overlaps_intron_starts <- GenomicRanges::countOverlaps(
      intron_starts.GRanges,
      alignment,
      minoverlap = 8,
      type = "any"
    )

    overlaps_intron_ends <- GenomicRanges::countOverlaps(
      intron_ends.GRanges,
      alignment,
      minoverlap = 8,
      type = "any"
    )

    GenomicRanges::mcols(introns.GRanges)["ir_s"] <- overlaps_intron_starts
    GenomicRanges::mcols(introns.GRanges)["ir_e"] <- overlaps_intron_ends
    GenomicRanges::mcols(introns.GRanges)["ir"] <- (overlaps_intron_starts + overlaps_intron_ends) / 2
    ir[[sample_name]] <- introns.GRanges
    sj[[sample_name]] <- GenomicAlignments::summarizeJunctions(alignment, genome = Genome_Assembly)
    BiocGenerics::strand(sj[[sample_name]]) <- GenomicRanges::mcols(sj[[sample_name]])[, "intron_strand"]
  }
  #--Aggregate and count junctional reads for each sample------------------------
  combined_sj <- unique(unlist(GenomicRanges::GRangesList(unlist(sj))))
  combined_ir <- unique(unlist(GenomicRanges::GRangesList(unlist(ir))))

  GenomicRanges::mcols(combined_sj)[c(
    "score", "plus_score", "minus_score", "intron_motif",
    "intron_strand"
  )] <- NULL
  GenomicRanges::mcols(combined_ir)[c("ir", "intron_no", "genes")] <- NULL

  for (sample_number in 1:length(sample_names)) {
    sample_name <- sample_names[sample_number]
    # message("\t", sample_name)

    GenomicRanges::mcols(combined_sj)$SJ_IR <- "SJ"
    GenomicRanges::mcols(combined_sj)[paste0("count_", sample_name)] <- 0
    qryhits <- GenomicRanges::findOverlaps(sj[[sample_name]], combined_sj, type = "equal")
    GenomicRanges::mcols(combined_sj[S4Vectors::subjectHits(qryhits)])[paste0("count_", sample_name)] <-
      GenomicRanges::mcols(sj[[sample_name]][S4Vectors::queryHits(qryhits)])[, "score"]

    GenomicRanges::mcols(combined_ir)$SJ_IR <- "IR"
    GenomicRanges::mcols(combined_ir)[paste0("count_", sample_name)] <- 0
    qryhits <- GenomicRanges::findOverlaps(ir[[sample_name]], combined_ir, type = "equal")
    GenomicRanges::mcols(combined_ir[S4Vectors::queryHits(qryhits)])[paste0("count_", sample_name)] <-
      GenomicRanges::mcols(ir[[sample_name]][S4Vectors::queryHits(qryhits)])[, "ir"]
  }
  combined_sj <- c(combined_sj, combined_ir)
  message("")
  return(combined_sj)
}
