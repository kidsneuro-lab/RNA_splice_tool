#' Annotate and Quantify Splicing Events
#'
#' Internal function that annotates detected splicing junctions and quantifies
#' their relative proportions. Identifies canonical splicing, alternative splicing,
#' exon skipping, and cryptic splice site usage.
#'
#' @param ids Character vector of sample identifiers
#' @param combined_sj GRanges object with combined splice junction data
#' @param introns.GRanges GRanges object with canonical intron coordinates
#' @param introns_other_tx.GRanges GRanges object with alternative transcript introns
#' @param introns Data.table with intron metadata
#' @param assembly Genome assembly version ("hg38" or "hg19")
#' @param debug Debug parameter (path or FALSE)
#' @param ria Logical indicating whether to include "reads in absentia"
#'
#' @return A data.table with annotated and quantified splicing events
#' @keywords internal
annotateQuantifyEvents <- function(ids, combined_sj, introns.GRanges, introns_other_tx.GRanges, introns, assembly, debug, ria) {

  message("Annotating and quantifying events...")

  combined_sj_sorted <- GenomeInfoDb::sortSeqlevels(combined_sj)
  combined_sj_sorted <- sort(combined_sj_sorted)

  events_by_intron <- list()

  # Extract and Annotate All Events at Canonical Junctions
  for(query_intron in seq_len(nrow(introns))) {
    intron_name <- paste0(introns$gene_name[query_intron], " intron ", introns$region_no[query_intron])

    # For all samples extract all the reads which overlap the query intron jxns
    qryhits_start <- GenomicRanges::findOverlaps(introns.GRanges[query_intron], combined_sj_sorted, type = "start")
    qryhits_end <- GenomicRanges::findOverlaps(introns.GRanges[query_intron], combined_sj_sorted, type = "end")

    # Combine the overlapping reads for the query intron
    query_intron.GRanges <- c(
      combined_sj_sorted[S4Vectors::subjectHits(qryhits_start)],
      combined_sj_sorted[S4Vectors::subjectHits(qryhits_end)]
    )

    # Include reads in absentia
    if (isTRUE(ria)) {
      qryhits_within <- GenomicRanges::findOverlaps(introns.GRanges[query_intron], combined_sj_sorted, type = "within")
      query_intron.GRanges <- c(
        query_intron.GRanges,
        combined_sj_sorted[S4Vectors::subjectHits(qryhits_within)]
      )
    }

    # Identify the non-canonical splicing junctions
    GenomicRanges::mcols(query_intron.GRanges)$"intron_jxn_start" <- NA
    GenomicRanges::mcols(query_intron.GRanges)$"intron_jxn_end" <- NA

    qryhits_start_allintrons <- GenomicRanges::findOverlaps(introns.GRanges, query_intron.GRanges, type = "start")
    qryhits_end_allintrons <- GenomicRanges::findOverlaps(introns.GRanges, query_intron.GRanges, type = "end")

    GenomicRanges::mcols(query_intron.GRanges[S4Vectors::subjectHits(qryhits_start_allintrons)])$"intron_jxn_start" <-
      GenomicRanges::mcols(introns.GRanges[S4Vectors::queryHits(qryhits_start_allintrons)])[, "intron_no"]
    GenomicRanges::mcols(query_intron.GRanges[S4Vectors::subjectHits(qryhits_end_allintrons)])$"intron_jxn_end" <-
      GenomicRanges::mcols(introns.GRanges[S4Vectors::queryHits(qryhits_end_allintrons)])[, "intron_no"]

    # Annotate cryptics that are actually other isoforms
    GenomicRanges::mcols(query_intron.GRanges)$"annotated" <- ""

    qryhits_equal_othertx <- GenomicRanges::findOverlaps(introns_other_tx.GRanges, query_intron.GRanges, type = "equal")

    if (length(qryhits_equal_othertx) > 0) {
      GenomicRanges::mcols(query_intron.GRanges[S4Vectors::subjectHits(qryhits_equal_othertx)])$"annotated" <- "alternative"
    }

    # Annotate canonical splice junctions
    qryhits_equal_allintrons <- GenomicRanges::findOverlaps(introns.GRanges, query_intron.GRanges, type = "equal")

    GenomicRanges::mcols(query_intron.GRanges[S4Vectors::subjectHits(qryhits_equal_allintrons)])$"annotated" <- "canonical"

    # Convert intron GRange into data.table
    query_intron.dt <- data.table::as.data.table(query_intron.GRanges)
    query_intron.dt <- unique(query_intron.dt)

    # Calculate the proportion of splicing at the intron each event represents
    for (sample_name in ids) {
      query_intron.dt[, paste0("pct_", sample_name)] <-
        data.table::nafill((query_intron.dt[[paste0("count_", sample_name)]] /
          sum(query_intron.dt[[paste0("count_", sample_name)]])),
        type = "const",
        fill = 0,
        nan = NA
        )
    }

    # Annotate with assembly, gene, intron_no
    query_intron.dt$assembly <- assembly
    query_intron.dt$gene <- introns$gene_name[query_intron]
    query_intron.dt$intron_no <- introns$region_no[query_intron]

    # Annotate with event
    query_intron.dt$event <- eventAnnotation(query_intron.dt)

    # Annotate with frame
    query_intron.dt$frame_conserved <- framed(query_intron.dt, assembly)

    # Append intron analysis to all events by intron
    events_by_intron[[intron_name]] <- query_intron.dt
  }


  # Merge all introns into a single data.table
  all_splicing_events <- data.table::rbindlist(events_by_intron)

  message("")

  if(is_debug_enabled(debug)){
    fwrite(as.data.table(all_splicing_events),paste0(debug,"/","7_all_splicing_events.tsv"), sep = "\t")
  }

  return(all_splicing_events)
}

#' Annotate Splicing Events
#'
#' Internal helper function that converts event information into human-readable
#' format (e.g., "canonical splicing", "exon skipping", "cryptic splice-site use").
#'
#' @param query_intron.dt A data.table containing splicing event information
#'
#' @return A character vector with annotated event descriptions
#' @keywords internal
#Takes a data-table with events and returns them in a human readable format
eventAnnotation <- function(query_intron.dt){
    n_events <- nrow(query_intron.dt)
    if (n_events == 0) {
      return(character(0))
    }

    start_jxn <- query_intron.dt$intron_jxn_start
    end_jxn <- query_intron.dt$intron_jxn_end
    strand <- as.character(query_intron.dt$strand)
    sj_ir <- query_intron.dt$SJ_IR

    min_range <- pmin(start_jxn, end_jxn)
    max_range <- pmax(start_jxn, end_jxn)
    range_diff <- max_range - min_range

    skip_labels <- vapply(seq_len(n_events), function(i) {
      if (!is.na(min_range[i]) && !is.na(max_range[i]) && range_diff[i] > 0) {
        paste0(
          "exon ",
          paste0(seq.int(min_range[i] + 1, max_range[i]), collapse = "-"),
          " skipping"
        )
      } else {
        ""
      }
    }, character(1))

    events <- data.table::fcase(
      !is.na(start_jxn) & !is.na(end_jxn) & range_diff == 0 & sj_ir == "SJ",
      paste0("canonical exon ", min_range, "-", max_range + 1, " splicing"),
      !is.na(start_jxn) & !is.na(end_jxn) & range_diff == 0 & sj_ir != "SJ",
      paste0("intron ", min_range, " retention"),
      !is.na(start_jxn) & !is.na(end_jxn) & range_diff > 0,
      skip_labels,
      !is.na(start_jxn) & is.na(end_jxn) & strand == "-",
      paste0("cryptic donor ~ exon ", start_jxn + 1),
      !is.na(start_jxn) & is.na(end_jxn) & strand == "+",
      paste0("exon ", start_jxn, " ~ cryptic acceptor"),
      !is.na(start_jxn) & is.na(end_jxn) & strand == "*",
      "cryptic (strand unknown)",
      is.na(start_jxn) & !is.na(end_jxn) & strand == "-",
      paste0("exon ", end_jxn, " ~ cryptic acceptor"),
      is.na(start_jxn) & !is.na(end_jxn) & strand == "+",
      paste0("cryptic donor ~ exon ", end_jxn + 1),
      is.na(start_jxn) & !is.na(end_jxn) & strand == "*",
      "cryptic (strand unknown)",
      is.na(start_jxn) & is.na(end_jxn),
      "unannotated junctions",
      default = NA_character_
    )

    return(events)
}

#' Calculate Frame Conservation for Splicing Events
#'
#' Internal helper function that determines whether a splicing event maintains
#' the reading frame.
#'
#' @param query_intron.dt A data.table containing splicing event information
#' @param assembly Genome assembly version ("hg38" or "hg19")
#'
#' @return A vector indicating frame conservation status (TRUE/FALSE/NA/"")
#' @keywords internal
#Calculate event frame - this could be pre-computed for annotated
framed <- function(query_intron.dt, assembly){

    if (assembly == "hg38") {
        rfsq <- refseq_introns_exons_hg38
    } else if (assembly == "hg19") {
        rfsq <- refseq_introns_exons_hg19
    } else {
        stop("Unsupported assembly: ", assembly)
    }

    n_events <- nrow(query_intron.dt)
    if (n_events == 0) {
      return(character(0))
    }

    frame <- rep("", n_events)
    is_sj <- query_intron.dt$SJ_IR == "SJ"
    starts <- query_intron.dt$start
    ends <- query_intron.dt$end

    start_annotated <- starts %in% rfsq$region_start
    end_annotated <- ends %in% rfsq$region_end

    for(event in which(is_sj)){
      if(start_annotated[event] && end_annotated[event]){
        pairend <- unique(rfsq$region_end[rfsq$region_start == starts[event]])
        if(ends[event] %in% pairend){
          frame[event] <- TRUE
        }else{
          dist2authentic <- abs(ends[event] - pairend)
          frame[event] <- any(dist2authentic %% 3 == 0)
        }
      } else if(start_annotated[event]){
        pairend <- unique(rfsq$region_end[rfsq$region_start == starts[event]])
        dist2authentic <- abs(ends[event] - pairend)
        frame[event] <- any(dist2authentic %% 3 == 0)
      } else if(end_annotated[event]){
        pairstart <- unique(rfsq$region_start[rfsq$region_end == ends[event]])
        dist2authentic <- abs(starts[event] - pairstart)
        frame[event] <- any(dist2authentic %% 3 == 0)
      } else{
        frame[event] <- NA
      }
    }

    return(frame)
}
