compareSplicing <- function(all_splicing_events, Sample_File, mode) {
    message("Comparing samples...")
    comparisons <- list()

if(mode == "default" | mode == "panel"){
  #--Compare splicing between test and controls and Generate Report--------------
  for (sample_number in seq(1, nrow(Sample_File))) {
    # Initialise new copy of the all_splicing_events dataset
    all_splicing_events_sample <- all_splicing_events
    if (Sample_File$sampletype[sample_number] == "test"){
      # Identify columns for the proband and family members
      proband <- Sample_File$sampleID[sample_number]
      family <- Sample_File$sampleID[which(
        Sample_File$family == Sample_File$family[sample_number]
      )]
      message("\t", proband)
      familycols <- paste0("pct_", family)
      familyreadcols <- paste0("count_", family)

      # Identify columns for the controls
      ctrls <- Sample_File$sampleID[which(
        Sample_File$family != Sample_File$family[sample_number]
      )]
      ctrlscols <- paste0("pct_", ctrls)
      ctrlsreadcols <- paste0("count_", ctrls)

      # Initialise various columns
      # Proband
      all_splicing_events_sample$proband <- Sample_File$sampleID[sample_number]

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
          if (all_splicing_events_sample$controlavg[i] == 0 &
            all_splicing_events_sample[, ..member][i] != 0) {
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
  return(comparisons)
}

