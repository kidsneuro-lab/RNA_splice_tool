#' Select genes and transcripts and build genomic ranges
#'
#' `selectGenesTranscripts()` resolves the requested genes/transcripts and
#' returns the genomic ranges used by downstream cortar steps.
#'
#' @param genes Character vector of identifiers. Supported formats are gene
#'   symbols, Entrez Gene IDs, RefSeq transcript IDs with version
#'   (for example `"NM_000546.6"`), and RefSeq transcript IDs without version
#'   (for example `"NM_000546"`).
#' @param assembly Genome assembly used for alignment. Supported values are
#'   `"hg38"` and `"hg19"`.
#' @param annotation Sequence naming style used by the aligned BAM files.
#'   Supported values are `"UCSC"` and `"1000genomes"`.
#' @param debug Optional output directory for writing intermediate TSV files.
#'   Use `""` (default) to disable debug output.
#'
#' @returns A named list with:
#'   * `genes`: `GRanges` covering selected canonical transcripts.
#'   * `introns`: list with `granges` (`GRanges`) and `metadata`
#'     (`data.table`) for introns in selected transcripts.
#'   * `introns_other_tx`: `GRanges` of introns from other transcripts of the
#'     same genes.
#'   * `junctions`: combined `GRanges` of intron start/end junction windows.
#'
#' @examples
#' \dontrun{
#' # Example 1: single-gene selection (mirrors test fixture usage)
#' res <- selectGenesTranscripts(c("EMD"), "hg38", "UCSC")
#' res$genes
#' length(res$introns$granges)
#'
#' # Example 2: alternate sequence style for a transcript-rich gene
#' res2 <- selectGenesTranscripts(c("COL2A1"), "hg38", "1000genomes")
#' length(res2$introns_other_tx)
#' length(res2$junctions)
#' }
#'
#' @export
selectGenesTranscripts <- function(genes,
                                   assembly,
                                   annotation,
                                   debug = "") {
  # Initialisation message
  futile.logger::flog.info("Selecting genes and transcripts...")
  futile.logger::flog.debug(sprintf("assembly: %s", assembly))
  futile.logger::flog.debug(sprintf("annotation: %s", annotation))

  if (assembly == "hg38") {
    refseq_introns_exons <- refseq_introns_exons_hg38
  } else if (assembly == "hg19") {
    refseq_introns_exons <- refseq_introns_exons_hg19
  }

  gene_tx <- tx_extraction(genes, refseq_introns_exons, debug)
  genes.GRanges <- gene_to_GRange(gene_tx, annotation, refseq_introns_exons, debug)
  introns.GRanges <- introns_to_GRange(gene_tx, annotation, refseq_introns_exons, debug)
  introns_other_tx.GRanges <- introns_other_tx_to_GRange(gene_tx, annotation, refseq_introns_exons, debug)
  introns_jx.GRanges <- introns_jx_to_GRange(gene_tx, annotation, refseq_introns_exons, debug)

  return(list(
    genes = genes.GRanges,
    introns = introns.GRanges,
    introns_other_tx = introns_other_tx.GRanges,
    junctions = unlist(introns_jx.GRanges)
  ))
}

#' Resolve input identifiers to RefSeq transcripts
#'
#' `tx_extraction()` maps each input identifier to a RefSeq transcript record.
#' Gene-level identifiers return canonical transcripts (`canonical == 1L`),
#' while transcript identifiers return the matching transcript entry.
#'
#' @param genes Character vector of gene symbols, Entrez Gene IDs, RefSeq
#'   transcript IDs with version, or RefSeq transcript IDs without version.
#' @param refseq_introns_exons `data.table` of RefSeq exon/intron annotations
#'   for a genome assembly (for example `refseq_introns_exons_hg38`).
#' @param debug Optional output directory for writing intermediate TSV files.
#'   Use `""` (default) to disable debug output.
#'
#' @returns A `data.table` with one row per input identifier and columns:
#'   * `gene_name`: resolved HGNC-style gene symbol.
#'   * `tx_version_id`: resolved RefSeq transcript ID with version.
#'
#' @details
#' Identifier matching uses indexed/keyed `data.table` joins for efficient
#' lookup across large annotation tables.
#'
#' @examples
#' \dontrun{
#' # Example 1: gene symbols
#' tx_extraction(c("EMD", "COL2A1"), refseq_introns_exons_hg38)
#'
#' # Example 2: Entrez Gene ID + transcript ID (with version)
#' tx_extraction(c(2010, "NM_033150.3"), refseq_introns_exons_hg38)
#' }
#'
#' @export
tx_extraction <- function(genes,
                          refseq_introns_exons,
                          debug = "") {
  refseq_genes <- refseq_introns_exons[region_type == 'exon',
                                      .(start = min(region_start),
                                        end = max(region_end)), by = .(gene_id,
                                                                       gene_name,
                                                                       tx_id,
                                                                       tx_version_id,
                                                                       canonical,
                                                                       chrom)]

  # PERF: add keys/indices so lookups are not full scans each iteration
  data.table::setkey(refseq_genes, gene_id)
  data.table::setindex(refseq_genes, tx_version_id)
  data.table::setindex(refseq_genes, tx_id)
  data.table::setindex(refseq_genes, gene_name)

  genes_info <- list()

  for (counter in seq_along(genes)) {
    gene_or_transcript <- genes[[counter]]

    # Skip empty lines in input samplesheet where Gene is ""
    # This can occur in default mode where only a single sample has a gene associated with it
    if ((is.character(gene_or_transcript) & gene_or_transcript == "") | is.na(gene_or_transcript)) {
      next
    }

    futile.logger::flog.debug(sprintf("Gene or Transcript: %s", as.character(gene_or_transcript)))

    if (length(gene_or_transcript) > 1) {
      stop(sprintf("Element %d has length > 1", counter))
    }

    # Coerce to character if numeric (e.g. Entrez Gene IDs passed as integers)
    isGeneID    <- grepl("^\\d+$", gene_or_transcript)
    isTranscriptWithVersion <- grepl("^NM_[0-9]+\\.[0-9]+$", gene_or_transcript)
    isTranscriptNoVersion <- grepl("^NM_[0-9]+$", gene_or_transcript)
    isGeneName  <- grepl("^[A-Za-z0-9]{3,10}$",   gene_or_transcript)

    gene_or_transcript <- ifelse(isGeneID, as.numeric(gene_or_transcript), gene_or_transcript)

    futile.logger::flog.debug(sprintf("isGeneID: %s | isTranscriptWithVersion: %s | isTranscriptNoVersion: %s | isGeneName: %s",
                                      isGeneID, isTranscriptWithVersion, isTranscriptNoVersion, isGeneName))

    # PERF: else-if chain prevents multiple scans/joins for the same element
    # data.table join syntax:
    #   refseq_genes[.(gene_or_transcript), on = "gene_id"]
    #
    # .() is shorthand for list() and creates a one-row lookup table.
    # This performs a join of refseq_genes to that value using gene_id.
    #
    # Because gene_id is keyed / indexed, data.table uses fast binary search
    # instead of scanning the entire table.
    #
    # This is significantly faster than:
    #   refseq_genes[gene_id == gene_or_transcript]
    # which may trigger a full table scan for each lookup.
    #
    # In large tables or inside loops, the join syntax reduces
    # repeated O(n) scans to O(log n) lookups.

    if (isGeneID) {
      genes_info[[counter]] <-
        refseq_genes[.(gene_or_transcript), on = "gene_id"][canonical == 1L,
                                                            .(gene_name, tx_version_id)
        ]

    } else if (isTranscriptWithVersion) {
      genes_info[[counter]] <-
        refseq_genes[.(gene_or_transcript), on = "tx_version_id",
                     .(gene_name, tx_version_id)
        ]

    } else if (isTranscriptNoVersion) {
      genes_info[[counter]] <-
        refseq_genes[.(gene_or_transcript), on = "tx_id",
                     .(gene_name, tx_version_id)
        ]

    } else if (isGeneName) {
      genes_info[[counter]] <-
        refseq_genes[.(gene_or_transcript), on = "gene_name"][canonical == 1L,
                                                              .(gene_name, tx_version_id)
        ]

    }

    if (is.na(genes_info[[counter]]$gene_name) | is.na(genes_info[[counter]]$tx)) {
      stop(sprintf("Gene or Transcript identifier `%s` is invalid", gene_or_transcript))
    }
  }

  if(is_debug_enabled(debug)){
    fwrite(as.data.table(genes_tx),paste0(debug,"/","1_tx_extraction.tsv"), sep = "\t")
  }

  data.table::rbindlist(genes_info, use.names = FALSE, fill = FALSE)
}

#' Convert selected transcripts to gene-level genomic ranges
#'
#' Internal helper that converts selected transcript rows into gene-level
#' `GRanges` spanning the minimum and maximum exon coordinates.
#'
#' @param gene_tx `data.table` with at least `gene_name` and `tx_version_id`
#'   columns, typically returned by [tx_extraction()].
#' @param annotation Sequence naming style used by aligned BAM files. If set to
#'   `"UCSC"`, seqlevels are converted to UCSC style.
#' @param refseq_introns_exons `data.table` containing RefSeq exon/intron
#'   annotations.
#' @param debug Optional output directory for writing intermediate TSV files.
#'   Use `FALSE` (default) to disable debug output.
#'
#' @returns A `GRanges` object with one range per selected transcript/gene row.
#'
#' @examples
#' \dontrun{
#' gene_tx <- data.table::data.table(
#'   gene_name = c("EMD", "COL2A1"),
#'   tx_version_id = c("NM_000117.3", "NM_001844.5")
#' )
#'
#' # Example 1: UCSC chromosome naming
#' cortar:::gene_to_GRange(gene_tx, "UCSC", refseq_introns_exons_hg38)
#'
#' # Example 2: 1000 Genomes-style chromosome naming
#' cortar:::gene_to_GRange(gene_tx, "1000genomes", refseq_introns_exons_hg38)
#' }
#' @keywords internal
gene_to_GRange <- function(gene_tx, annotation, refseq_introns_exons, debug = FALSE) {
  if (nrow(gene_tx) == 0) {
    stop("Cannot obtain GRange for empty data.table")
  }

  refseq_genes <- refseq_introns_exons[region_type == 'exon',
                                      .(start = min(region_start),
                                        end = max(region_end)), by = .(gene_id,
                                                                       gene_name,
                                                                       tx_id,
                                                                       tx_version_id,
                                                                       canonical,
                                                                       strand,
                                                                       chrom)]

  genomic_ranges_for_selected_genes <- refseq_genes[gene_tx, on = .(gene_name, tx_version_id)]

  genes.GRanges <- GenomicRanges::GRanges(
    seqnames = genomic_ranges_for_selected_genes$chrom,
    IRanges::IRanges(
      start = genomic_ranges_for_selected_genes$start,
      end = genomic_ranges_for_selected_genes$end,
    ),
    strand = genomic_ranges_for_selected_genes$strand
  )

  if (annotation == "UCSC") {
    GenomeInfoDb::seqlevelsStyle(genes.GRanges) <- "UCSC"
  }

  if(is_debug_enabled(debug)){
    fwrite(as.data.table(genes.GRanges),paste0(debug,"/","2_gene_to_GRange.tsv"), sep = "\t")
  }

  return(genes.GRanges)
}

#' Convert selected-transcript introns to genomic ranges
#'
#' Internal helper that extracts introns from selected transcripts and returns
#' both range and tabular representations.
#'
#' @param gene_tx `data.table` with transcript selection, including
#'   `tx_version_id`.
#' @param annotation Sequence naming style used by aligned BAM files. If set to
#'   `"UCSC"`, seqlevels are converted to UCSC style.
#' @param refseq_intron_exons `data.table` containing RefSeq exon/intron
#'   annotations.
#' @param debug Optional output directory for writing intermediate TSV files.
#'   Use `FALSE` (default) to disable debug output.
#'
#' @returns A list with:
#'   * `granges`: intron `GRanges` for selected transcripts.
#'   * `metadata`: intron `data.table` used to construct `granges`.
#'
#' @examples
#' \dontrun{
#' # Example 1: introns for EMD (UCSC style)
#' gene_tx <- data.table::data.table(
#'   gene_name = "EMD",
#'   tx_version_id = "NM_000117.3"
#' )
#' introns <- cortar:::introns_to_GRange(gene_tx, "UCSC", refseq_introns_exons_hg38)
#' introns$granges
#'
#' # Example 2: same transcript with 1000 Genomes sequence style
#' introns2 <- cortar:::introns_to_GRange(gene_tx, "1000genomes", refseq_introns_exons_hg38)
#' introns2$metadata[, .N]
#' }
#' @keywords internal
introns_to_GRange <- function(gene_tx, annotation, refseq_intron_exons, debug = FALSE) {
  introns <- refseq_intron_exons[tx_version_id %in% gene_tx$tx & region_type == c("intron")]

  introns.GRanges <- GenomicRanges::GRanges(
    seqnames = introns$chrom,
    IRanges::IRanges(
      start = introns$region_start,
      end = introns$region_end
    ),
    strand = introns$strand
  )

  GenomicRanges::mcols(introns.GRanges)["intron_no"] <- introns$region_no
  GenomicRanges::mcols(introns.GRanges)["gene"] <- introns$gene_name

  if (annotation == "UCSC") {
    GenomeInfoDb::seqlevelsStyle(introns.GRanges) <- "UCSC"
  }

  if(is_debug_enabled(debug)){
    fwrite(as.data.table(introns.GRanges),paste0(debug,"/","3_introns_to_GRange_GRanges.tsv"), sep = "\t")
    fwrite(as.data.table(introns),paste0(debug,"/","3_introns_to_GRange_introns.tsv"), sep = "\t")
  }

  return(list(granges = introns.GRanges, metadata = introns))
}

#' Convert alternative-transcript introns to genomic ranges
#'
#' Internal helper that retrieves introns from non-selected transcripts within
#' the same genes as `gene_tx`.
#'
#' @param gene_tx `data.table` with selected gene/transcript rows, including
#'   `gene_name` and `tx_version_id`.
#' @param annotation Sequence naming style used by aligned BAM files. If set to
#'   `"UCSC"`, seqlevels are converted to UCSC style.
#' @param refseq_intron_exons `data.table` containing RefSeq exon/intron
#'   annotations.
#' @param debug Optional output directory for writing intermediate TSV files.
#'   Use `FALSE` (default) to disable debug output.
#'
#' @returns A `GRanges` object of introns from other transcripts, with metadata
#' columns `tx_id`, `intron_no`, and `gene`.
#'
#' @examples
#' \dontrun{
#' gene_tx <- data.table::data.table(
#'   gene_name = "COL2A1",
#'   tx_version_id = "NM_001844.5"
#' )
#'
#' # Example 1: alternative-transcript introns in UCSC style
#' other_introns <- cortar:::introns_other_tx_to_GRange(
#'   gene_tx, "UCSC", refseq_introns_exons_hg38
#' )
#' other_introns
#'
#' # Example 2: same query with 1000 Genomes sequence style
#' other_introns2 <- cortar:::introns_other_tx_to_GRange(
#'   gene_tx, "1000genomes", refseq_introns_exons_hg38
#' )
#' length(other_introns2)
#' }
#' @keywords internal
introns_other_tx_to_GRange <- function(gene_tx, annotation, refseq_intron_exons, debug = FALSE) {
  introns <- refseq_intron_exons[gene_name %in% gene_tx$gene_name & tx_version_id %nin% gene_tx$tx_version_id & region_type == c("intron")]

  introns.GRanges <- GenomicRanges::GRanges(
    seqnames = introns$chrom,
    IRanges::IRanges(
      start = introns$region_start,
      end = introns$region_end
    ),
    strand = introns$strand
  )

  GenomicRanges::mcols(introns.GRanges)["tx_id"] <- introns$tx_version_id
  GenomicRanges::mcols(introns.GRanges)["intron_no"] <- introns$region_no
  GenomicRanges::mcols(introns.GRanges)["gene"] <- introns$gene_name

  if (annotation == "UCSC") {
    GenomeInfoDb::seqlevelsStyle(introns.GRanges) <- "UCSC"
  }

  if(is_debug_enabled(debug)){
    fwrite(as.data.table(introns.GRanges),paste0(debug,"/","4_introns_other_tx_to_GRange.tsv"), sep = "\t")
  }

  return(introns.GRanges)
}

#' Build intron-junction windows for selected transcripts
#'
#' Internal helper that creates `GRanges` windows around intron starts and ends
#' (exon-intron boundaries) using `INTRON_JUNCTION_UPSTREAM` and
#' `INTRON_JUNCTION_DOWNSTREAM`.
#'
#' @param gene_tx `data.table` with transcript selection, including
#'   `tx_version_id`.
#' @param annotation Sequence naming style used by aligned BAM files. If set to
#'   `"UCSC"`, seqlevels are converted to UCSC style.
#' @param refseq_introns_exons `data.table` containing RefSeq exon/intron
#'   annotations.
#' @param debug Optional output directory for writing intermediate TSV files.
#'   Use `FALSE` (default) to disable debug output.
#'
#' @returns A list with:
#'   * `starts`: `GRanges` centred on intron starts.
#'   * `ends`: `GRanges` centred on intron ends.
#'
#' @examples
#' \dontrun{
#' gene_tx <- data.table::data.table(
#'   gene_name = "COL2A1",
#'   tx_version_id = "NM_001844.5"
#' )
#'
#' # Example 1: junction windows for COL2A1 (UCSC style)
#' jx <- cortar:::introns_jx_to_GRange(gene_tx, "UCSC", refseq_introns_exons_hg38)
#' c(length(jx$starts), length(jx$ends))
#'
#' # Example 2: same query with 1000 Genomes sequence style
#' jx2 <- cortar:::introns_jx_to_GRange(gene_tx, "1000genomes", refseq_introns_exons_hg38)
#' jx2$starts
#' }
#' @keywords internal
introns_jx_to_GRange <- function(gene_tx, annotation, refseq_introns_exons, debug = FALSE) {
  introns <- refseq_introns_exons[tx_version_id %in% gene_tx$tx_version_id & region_type == c("intron")]

  intron_starts.GRanges <- GenomicRanges::GRanges(
    seqnames = introns$chrom,
    IRanges::IRanges(
      start = introns$region_start - INTRON_JUNCTION_UPSTREAM,
      end = introns$region_start + INTRON_JUNCTION_DOWNSTREAM
    ),
    strand = introns$strand
  )

  intron_ends.GRanges <- GenomicRanges::GRanges(
    seqnames = introns$chrom,
    IRanges::IRanges(
      start = introns$region_end - INTRON_JUNCTION_DOWNSTREAM,
      end = introns$region_end + INTRON_JUNCTION_UPSTREAM
    ),
    strand = introns$strand
  )

  if (annotation == "UCSC") {
    GenomeInfoDb::seqlevelsStyle(intron_starts.GRanges) <- "UCSC"
    GenomeInfoDb::seqlevelsStyle(intron_ends.GRanges) <- "UCSC"
  }

  if(is_debug_enabled(debug)){
    fwrite(as.data.table(intron_starts.GRanges),paste0(debug,"/","5_introns_jx_to_GRange_starts.tsv"), sep = "\t")
    fwrite(as.data.table(intron_ends.GRanges),paste0(debug,"/","5_introns_jx_to_GRange_ends.tsv"), sep = "\t")
  }

  return(list(starts = intron_starts.GRanges, ends = intron_ends.GRanges))
}
