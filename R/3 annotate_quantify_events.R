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

annotateQuantifyEvents <- function(ids, combined_sj, introns.GRanges, introns_other_tx.GRanges, introns, assembly) {

    message("Annotating and quantifying events...")
      #--Extract, annotate, and quantify all events at canonical junctions-----------
  combined_sj_sorted <- GenomeInfoDb::sortSeqlevels(combined_sj)
  combined_sj_sorted <- sort(combined_sj)

  events_by_intron <- list()

  # Extract and Annotate All Events at Canonical Junctions
  for(query_intron in seq(1, nrow(introns))) {
    intron_name <- paste0(introns$gene_name[query_intron], " intron ", introns$region_no[query_intron])
    # message(paste0("\t", intron_name))

    # For all samples extract all the reads which overlap the query intron jxns
    qryhits_start <- GenomicRanges::findOverlaps(introns.GRanges[query_intron], combined_sj_sorted, type = "start")
    qryhits_end <- GenomicRanges::findOverlaps(introns.GRanges[query_intron], combined_sj_sorted, type = "end")

    # Combine the overlapping reads for the query intron
    query_intron.GRanges <- c(
      combined_sj_sorted[S4Vectors::subjectHits(qryhits_start)],
      combined_sj_sorted[S4Vectors::subjectHits(qryhits_end)]
    )

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
    query_intron.dt$frame_conserved <- ""
    # query_intron.dt$frame_conserved <- framed(query_intron.dt, assembly)

    # Append intron analysis to all events by intron
    events_by_intron[[intron_name]] <- query_intron.dt
  }


  # Merge all introns into a single data.table
  all_splicing_events <- data.table::rbindlist(events_by_intron)

  message("")
  return(all_splicing_events)
}

#Takes a data-table with events and returns them in a human readable format
eventAnnotation <- function(query_intron.dt){
    event <- 1
    events <- c()
    for(event in seq(1,nrow(query_intron.dt))){
        #Extract intron start and end ranges and strand
        exon_range_start <- query_intron.dt[event,intron_jxn_start]
        exon_range_end <- query_intron.dt[event,intron_jxn_end]
        exon_ranges <- c(exon_range_start, exon_range_end)
        strand <- as.character(query_intron.dt$strand[event])

        #normal splicing
        if (!is.na(exon_range_start) & !is.na(exon_range_end)) {
            if(max(exon_ranges)-min(exon_ranges) == 0){
                if(query_intron.dt$SJ_IR[event] == "SJ"){
                    events[event] <- (paste0("canonical exon ",min(exon_ranges),"-",
                                             max(exon_ranges)+1," splicing"))

                    #intron retention
                }else{
                    events[event] <- (paste0("intron ",min(exon_ranges)," retention"))
                }

                #exon skipping
            }else{
                events[event] <- paste0("exon ", paste0(seq(min(exon_ranges)+1,max(exon_ranges)),collapse="-")," skipping")
            }

            #cryptic splice-site use
        }else if (!is.na(exon_range_start)){
            if (strand == "-"){
                events[event] <- (paste("CD", " ~ ", "exon ",
                                        exon_range_start+1, sep="", collapse=""))
            }
            if (strand == "+"){
                events[event] <- (paste("exon ",exon_range_start, " ~ ",
                                        "CA", sep="", collapse=""))
            }
        }else if (!is.na(exon_range_end)){
            if (strand == "-"){
                events[event] <- (paste("exon ", exon_range_end, " ~ ",
                                        "CA", sep="", collapse=""))
            }
            if (strand == "+"){
                events[event] <- (paste("CD", " ~ ", "exon ",
                                        exon_range_end+1, sep="", collapse=""))
            }

            #catching errors
        }else{
            events[event] <- ("unannotated junctions")
        }
    }
    return(events)
}

#Calculate event frame
framed <- function(query_intron.dt, assembly){

    if (assembly == "hg38") {
        rfsq <- refseq_introns_exons_hg38
    } else if (assembly == "hg19") {
        rfsq <- refseq_introns_exons_hg19.tsv
    }
    Mendelian_Intron_PTCs <- mendelian_intron_ptcs_hg38

    frame <- c()



    for(event in seq(1,nrow(query_intron.dt))){
        if(query_intron.dt$SJ_IR[event] == "SJ"){
            #is the start/end annotated
            #both annotated
            if(query_intron.dt$start[event] %in% rfsq$region_start & query_intron.dt$end[event] %in% rfsq$region_end){
                pairstart <- unique(rfsq$region_start[which(rfsq$region_end == query_intron.dt$end[event])])
                pairend <- unique(rfsq$region_end[which(rfsq$region_start == query_intron.dt$start[event])])
                if(query_intron.dt$end[event] == pairend){
                    frame[event] <- TRUE
                }else{
                    dist2authentic <- abs(query_intron.dt$end[event]-pairend)
                    frame[event] <- dist2authentic%%3 == 0
                }
            }
            #start annotated
            else if(query_intron.dt$start[event] %in% rfsq$region_start){
                pairend <- unique(rfsq$region_end[which(rfsq$region_start == query_intron.dt$start[event])])
                dist2authentic <- abs(query_intron.dt$end[event]-pairend)
                frame[event] <- dist2authentic%%3 == 0
            }
            #end annotated
            else if(query_intron.dt$end[event] %in% rfsq$region_end){
                pairstart <- unique(rfsq$region_start[which(rfsq$region_end == query_intron.dt$end[event])])
                dist2authentic <- abs(query_intron.dt$start[event]-pairstart)
                frame[event] <- dist2authentic%%3 == 0
            }
            #unannotated junctions
            else{
                frame[event] <- NA
            }

        # }else if(query_intron.dt$SJ_IR[event] == "IR"){
        #     frame[event] <- ""
        #     boundaries <<- c(query_intron.dt$start[event],query_intron.dt$end[event])
        #     for(gpos in seq(1,2)){
        #         #Find gene, canonical transcript, and overlapping regions
        #         mendelian_introns_gene <<- Mendelian_Intron_PTCs[gene == query_intron.dt$gene[event]]
        #         #Find gene, canonical transcript, overlapping region and strand
        #         for(regions in 1:nrow(Mendelian_Intron_PTCs)){
        #             if(between(boundaries[gpos],
        #                        mendelian_introns_gene$region_start[regions],
        #                        mendelian_introns_gene$region_end[regions])){
        #                 print(regions)
        #                 frame[event] <- mendelian_introns_gene$frame[regions]
        #                 if(is.na(frame[event])){
        #                     frame[event] <- ""
        #                     break
        #                 }else if(frame[event] == TRUE){
        #                     if(sum(mendelian_introns_gene[regions,c(window1ptc,window2ptc,window3ptc)]) == Inf){
        #                         frame[event] <- FALSE
        #                         break
        #                     }else{
        #                         frame[event] <- ""
        #                         break
        #                     }
        #                 }
        #             }
        #         }
        #     }
        }else{
            frame[event] <- ""
        }
    }
    return(frame)
}
