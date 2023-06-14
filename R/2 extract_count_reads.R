#' Extract and count split-reads from .bam files
#'
#' `selectGenesTranscripts` extracts split-reads from control and test .bam
#' files, counts them and aggregates the data into a single data table.
#'
#' @param genes.GRanges A GRanges object containing the coordinates of the
#'     genes of interest. Created by `selectGenesTranscripts()`
#' @param introns.GRanges A GRanges object containing the coordinates of the
#'     introns of the genes of interest. Created by `selectGenesTranscripts()`
#' @param intron_starts.GRanges A GRanges object containing the coordinates
#'     of the upstream exon-intron junctions of the genes of interest.
#'     Created by `selectGenesTranscripts()`
#' @param intron_ends.GRanges A GRanges object containing the coordinates
#'     of the downstream exon-intron junctions of the genes of interest.
#'     Created by `selectGenesTranscripts()`
#' @param bamfiles A character vector of file paths, pointing to the .bam files
#'     for analysis. Included in a properly formatted cortar samplefile.
#' @param sample_names A character vector of sample names. Included in a
#'     properly formatted cortar samplefile.
#' @param assembly Assembly used for alignment: either `"hg38"` or `"hg19"`
#' @param annotation  Annotation used for alignment: either `"UCSC"` or
#'     `"1000genomes"`
#' @param paired Is the RNA-seq paired-end?: `TRUE`/`FALSE`
#' @param stranded Strandedness of the RNA-seq: `0` for unstranded, `1`
#'     for forward stranded or `2` for reverse stranded
#'
#'
#' @returns A GRanges object consisting of aggregated junctional read count data
#' from all control and test samples. The GRanges object has at least three
#' metadata columns:
#' * SJ_IR: whether the entry is counting split (SJ) or non-split (IR) reads at
#' the junction.
#' * gene: the name of the gene which the junction occurs within or `NA` if the
#' junction is not with a gene
#' *count_sample_name: the count for each junction for a given sample. There
#' will be as many count_sample_name metadata columns as there are
#' sample_names specified
#' @export
#'
#' @examples
#' #### == COMING SOON == ####
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
  message("Extracting and counting reads...")

  if (assembly == "hg19") {
    if (annotation == "UCSC") {
      Genome_Assembly <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    } else if (annotation == "1000genomes") {
      Genome_Assembly <- BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
    }
  } else if (assembly == "hg38") {
    if (annotation == "UCSC"){
      Genome_Assembly <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    } else if (annotation == "NCBI"){
      Genome_Assembly <- BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
    }
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


  combined_sj <- unique(unlist(GenomicRanges::GRangesList(unlist(sj))))
  combined_ir <- unique(unlist(GenomicRanges::GRangesList(unlist(ir))))

  GenomicRanges::mcols(combined_sj)[c(
    "score", "plus_score", "minus_score", "intron_motif",
    "intron_strand"
  )] <- NULL
  GenomicRanges::mcols(combined_ir)[c("ir", "intron_no", "genes")] <- NULL

  for (sample_number in 1:length(sample_names)) {
    sample_name <- sample_names[sample_number]

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
  GenomicRanges::mcols(combined_sj)[c("ir_s", "ir_e")] <- NULL

  message("")
  return(combined_sj)
}
