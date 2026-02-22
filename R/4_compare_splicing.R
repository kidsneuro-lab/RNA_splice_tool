#' Compare Splicing Between Samples
#'
#' Internal function that performs statistical comparisons of splicing events
#' between test samples and controls. Calculates differences, standard deviations,
#' and identifies unique events.
#'
#' @param all_splicing_events Data.table containing all annotated splicing events
#' @param Sample_File Data.table with sample metadata (from cortar samplefile)
#' @param mode Analysis mode: "default", "panel", or "research"
#' @param debug Debug parameter (path or FALSE)
#'
#' @return A list of comparison results (one per test sample) or a single
#'   data.table for research mode
#' @keywords internal
compareSplicing <- function(all_splicing_events, Sample_File, mode, debug) {
    message("Comparing samples...")
    comparisons <- list()
    valid_modes <- c("default", "panel", "research")
    if (mode %nin% valid_modes) {
      stop(
        "Invalid mode '", mode, "'. Expected one of: ",
        paste(valid_modes, collapse = ", ")
      )
    }

    family_col <- if ("family" %in% names(Sample_File)) {
      "family"
    } else if ("familyID" %in% names(Sample_File)) {
      "familyID"
    } else {
      stop("Sample_File must contain either a 'family' or 'familyID' column.")
    }

    gene_col <- if ("gene" %in% names(Sample_File)) {
      "gene"
    } else if ("genes" %in% names(Sample_File)) {
      "genes"
    } else {
      stop("Sample_File must contain either a 'gene' or 'genes' column.")
    }

    family_values <- Sample_File[[family_col]]
    gene_values <- Sample_File[[gene_col]]

if(mode == "default" | mode == "panel"){
  #--Compare splicing between test and controls and Generate Report--------------
  for (sample_number in seq(1, nrow(Sample_File))) {

    # Initialise new copy of the all_splicing_events dataset
    all_splicing_events_sample <- all_splicing_events
    if (Sample_File$sampletype[sample_number] == "test"){

      # Identify columns for the proband and family members
      proband <- Sample_File$sampleID[sample_number]
      family <- Sample_File$sampleID[which(
        family_values == family_values[sample_number]
      )]
      message("\t", proband)
      familycols <- paste0("pct_", family)
      familyreadcols <- paste0("count_", family)

      # Identify columns for the controls
      ctrls <- Sample_File$sampleID[which(
        family_values != family_values[sample_number]
      )]
      if(mode == "default"){
        ctrls <- Sample_File$sampleID[which(
          family_values != family_values[sample_number] &
          gene_values != gene_values[sample_number]
        )]
      }
      ctrlscols <- paste0("pct_", ctrls)
      ctrlsreadcols <- paste0("count_", ctrls)

      #In default mode, add a filter to remove controls with a median coverage less than threshold over all canonical exons
      if(mode == "default"){
        canon_splicing_counts <- all_splicing_events_sample[gene == gene_values[sample_number] &
                                                            annotated == "canonical" &
                                                            SJ_IR == "SJ", ..ctrlsreadcols]

        coverage_threshold <- get_coverage_threshold(Sample_File$coverage[sample_number])
        control_medians <- sapply(canon_splicing_counts, median, na.rm = TRUE)
        keep_controls <- !is.na(control_medians) & control_medians > coverage_threshold
        ctrlscols <- ctrlscols[keep_controls]
        ctrlsreadcols <- ctrlsreadcols[keep_controls]

      }
      # Initialise various columns
      # Proband
      all_splicing_events_sample$proband <- Sample_File$sampleID[sample_number]

      # Control average pct, read count, sd, and n
      if (length(ctrlscols) == 0 || length(ctrlsreadcols) == 0) {
        all_splicing_events_sample$controlavg <- NA_real_
        all_splicing_events_sample$controlavgreads <- NA_real_
        all_splicing_events_sample$controlsd <- NA_real_
        all_splicing_events_sample$controln <- 0
      } else {
        all_splicing_events_sample$controlavg <- rowMeans(
          all_splicing_events_sample[, ..ctrlscols]
        )
        all_splicing_events_sample$controlavgreads <- rowMeans(
          all_splicing_events_sample[, ..ctrlsreadcols]
        )
        all_splicing_events_sample <- cbind(all_splicing_events_sample,
          controlsd = apply(all_splicing_events_sample[, ..ctrlscols], 1, sd)
        )
        all_splicing_events_sample$controln <- length(ctrlscols)
      }

      # Difference between proband and control average
      all_splicing_events_sample$difference <- all_splicing_events_sample[, paste0("pct_", proband), with = F] -
        all_splicing_events_sample$controlavg

      # Compute the standard deviation thresholds
      all_splicing_events_sample$two_sd <- abs(
        all_splicing_events_sample$difference
      ) > all_splicing_events_sample$controlsd * 2
      all_splicing_events_sample$three_sd <- abs(
        all_splicing_events_sample$difference
      ) > all_splicing_events_sample$controlsd * 3
      all_splicing_events_sample$four_sd <- abs(
        all_splicing_events_sample$difference
      ) >
        all_splicing_events_sample$controlsd * 4

      # Identify unique events
      for (i in seq(1, nrow(all_splicing_events_sample))) {
        event_unique_count <- 0
        for (member in familycols) {
          control_is_zero <- !is.na(all_splicing_events_sample$controlavg[i]) &&
            all_splicing_events_sample$controlavg[i] == 0
          member_value <- all_splicing_events_sample[, ..member][i]
          member_has_signal <- !is.na(member_value) && member_value != 0
          if (control_is_zero && member_has_signal) {
            event_unique_count <- event_unique_count + 1
          }
        }
        if (event_unique_count >= 1) {
          all_splicing_events_sample$unique[i] <- paste(
            event_unique_count, "/", length(familycols),
            sep = ""
          )
        } else {
          all_splicing_events_sample$unique[i] <- ""
        }
      }

      # Order columns and sort by the greatest difference
      all_splicing_events_sample <- all_splicing_events_sample[order(
        abs(all_splicing_events_sample$difference),
        decreasing = T
      ),
      c(
        "assembly",
        "proband",
        "seqnames",
        "start",
        "end",
        "width",
        "strand",
        "gene",
        "event",
        "annotated",
        "frame_conserved",
        "unique",
        "difference",
        familycols,
        "controlavg",
        "controlsd",
        "controln",
        familyreadcols,
        "controlavgreads",
        "two_sd",
        "three_sd",
        "four_sd",
        "intron_no",
        "SJ_IR"
      ),
      with = F
      ]
      comparisons[[sample_number]] <- all_splicing_events_sample
    }}
  }else if(mode == "research"){
      all_splicing_events_sample <- all_splicing_events
      ctrls <- Sample_File$sampleID
      ctrlscols <- paste0("pct_", ctrls)
      ctrlsreadcols <- paste0("count_", ctrls)

      # Control average pct, read count, sd, and n
      all_splicing_events_sample$controlavg <- rowMeans(
        all_splicing_events_sample[, ..ctrlscols]
      )
      all_splicing_events_sample$controlavgreads <- rowMeans(
        all_splicing_events_sample[, ..ctrlsreadcols]
      )
      all_splicing_events_sample <- cbind(all_splicing_events_sample,
                                          controlsd = apply(all_splicing_events_sample[, ..ctrlscols], 1, sd)
      )
      all_splicing_events_sample$controln <- length(ctrlscols)

    all_splicing_events_sample <- all_splicing_events_sample[order(
      gene,controlavg,
      decreasing = F
    ),
    c(
      "assembly",
      "seqnames",
      "start",
      "end",
      "width",
      "strand",
      "gene",
      "event",
      "annotated",
      "frame_conserved",
      "controlavg",
      "controlsd",
      "controln",
      "controlavgreads",
      "intron_no",
      "SJ_IR"
    ),
    with = F
    ]
    comparisons <- all_splicing_events_sample
  }
  message("")

  if(is_debug_enabled(debug)){
    fwrite(as.data.table(comparisons),paste0(debug,"/","8_comparisons.tsv"), sep = "\t")
  }

  return(comparisons)
}
