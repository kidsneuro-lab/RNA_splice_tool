#' Compare Splicing Between Samples
#'
#' Internal helper that filters control samples by median canonical-splice
#' coverage for the proband gene.
#'
#' @param events_dt Data.table containing annotated splicing events
#' @param gene_name Character string of the proband gene
#' @param ctrl_pct_cols Character vector of control pct columns
#' @param ctrl_read_cols Character vector of control read-count columns
#' @param coverage_type Character string of coverage type
#'
#' @return A named list containing filtered control column vectors
#' @keywords internal
filter_controls_by_coverage <- function(events_dt,
                                        gene_name,
                                        ctrl_pct_cols,
                                        ctrl_read_cols,
                                        coverage_type) {
  if (length(ctrl_pct_cols) == 0 || length(ctrl_read_cols) == 0) {
    return(list(
      ctrl_pct_cols = ctrl_pct_cols,
      ctrl_read_cols = ctrl_read_cols
    ))
  }

  canon_splicing_counts <- events_dt[
    gene == gene_name & annotated == "canonical" & SJ_IR == "SJ",
    ..ctrl_read_cols
  ]

  coverage_threshold <- get_coverage_threshold(coverage_type)
  control_medians <- sapply(canon_splicing_counts, median, na.rm = TRUE)
  keep_controls <- !is.na(control_medians) & control_medians > coverage_threshold

  list(
    ctrl_pct_cols = ctrl_pct_cols[keep_controls],
    ctrl_read_cols = ctrl_read_cols[keep_controls]
  )
}

#' Identify Unique Splicing Events
#'
#' Determines which splicing events are present in family members but absent in
#' controls (control average is 0 and family member value is non-zero).
#'
#' @param events_dt A data.table with splicing event data
#' @param family_pct_cols Character vector of family pct column names
#'
#' @return Character vector indicating family uniqueness per event
#' @keywords internal
calculate_unique_events <- function(events_dt, family_pct_cols) {
  n_events <- nrow(events_dt)
  n_family <- length(family_pct_cols)
  if (n_events == 0 || n_family == 0) {
    return(character(n_events))
  }

  family_matrix <- as.matrix(events_dt[, ..family_pct_cols])
  member_has_signal <- !is.na(family_matrix) & family_matrix != 0
  control_is_zero <- !is.na(events_dt$controlavg) & events_dt$controlavg == 0
  control_is_zero_matrix <- matrix(control_is_zero, nrow = n_events, ncol = n_family)

  unique_counts <- rowSums(control_is_zero_matrix & member_has_signal)

  return(ifelse(
    unique_counts >= 1,
    paste0(unique_counts, "/", n_family),
    ""
  ))
}

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

if (mode == "default" || mode == "panel") {
  #--Compare splicing between test and controls and Generate Report--------------
  for (sample_number in seq_len(nrow(Sample_File))) {

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
        filtered_controls <- filter_controls_by_coverage(
          events_dt = all_splicing_events_sample,
          gene_name = gene_values[sample_number],
          ctrl_pct_cols = ctrlscols,
          ctrl_read_cols = ctrlsreadcols,
          coverage_type = Sample_File$coverage[sample_number]
        )
        ctrlscols <- filtered_controls$ctrl_pct_cols
        ctrlsreadcols <- filtered_controls$ctrl_read_cols

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
      all_splicing_events_sample$difference <- all_splicing_events_sample[, paste0("pct_", proband), with = FALSE] -
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
      all_splicing_events_sample$unique <- calculate_unique_events(
        all_splicing_events_sample,
        familycols
      )

      # Order columns and sort by the greatest difference
      all_splicing_events_sample <- all_splicing_events_sample[order(
        abs(all_splicing_events_sample$difference),
        decreasing = TRUE
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
      with = FALSE
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
      decreasing = FALSE
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
    with = FALSE
    ]
    comparisons <- all_splicing_events_sample
  }
  message("")

  if(is_debug_enabled(debug)){
    fwrite(as.data.table(comparisons),paste0(debug,"/","8_comparisons.tsv"), sep = "\t")
  }

  return(comparisons)
}
