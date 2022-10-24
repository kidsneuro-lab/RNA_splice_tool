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
generateReport <- function(comparisons, Sample_File, Export) {
    message("Generating reports...")
#--Compare splicing between test and controls and Generate Report--------------
  for (sample_number in seq(1, nrow(Sample_File))) {
    all_splicing_events_sample <- comparisons[[sample_number]]
    # Initialise new copy of the all_splicing_events dataset
    if (Sample_File$sampletype[sample_number] == "test") {
      proband <- Sample_File$sampleID[sample_number]
      family <- Sample_File$sampleID[which(
        Sample_File$family == Sample_File$family[sample_number]
      )]
      message("\t", family)
      familycols <- paste0("pct_", family)
      familyreadcols <- paste0("count_", family)

      # Set report outputs and parameters
      report <- T
      splicing_diagnostics_report <- F
      full_all_genes_report <- F
      testgenes <- unique(Sample_File$genes[sample_number])
      # normalSpliceMap(all_splicing_events_sample, familycols[1], proband, testgenes, export = Export)
      # normalSpliceMap(all_splicing_events_sample, "norm_proband", testgenes)

      # Generate filtered excel spreadsheet +/- summary html
      # Excel spreadsheet
      if (report == TRUE) {
        generate.excel(all_splicing_events_sample[gene %in% testgenes],
          length(familycols),
          gene = testgenes,
          Export,
          Sample_File$sampleID[sample_number]
        )
      if (report == FALSE) {
        data.table::fwrite(all_splicing_events_sample, # [which(combined_dt_intron_test$genes == testgenes)],
          file = paste0(
            Export, "/", Sample_File$sampleID[sample_number], "_", testgenes, "_combined_full",
            ".tsv"
          ), sep = "\t"
        )
      }

        #--HTML summary report---- You are up to here! Nearly done! -----------
        if (splicing_diagnostics_report == TRUE) {
          report_table <- all_splicing_events_sample

          names(report_table)[which(names(report_table) == paste0("pct_", proband))] <- "probandpct"

          report_table$seqnames <- paste0(report_table$seqnames, ":", report_table$start, "-", report_table$end)
          report_table$difference <- formattable::percent(report_table$difference, digits = 1)
          report_table$probandpct <- formattable::percent(report_table$probandpct, digits = 1)
          report_table$controlavg <- formattable::percent(report_table$controlavg, digits = 1)


          sig_introns <- unique(report_table[two_sd == TRUE & abs(difference) > 0.05, intron_no])

          no_introns <- length(sig_introns)
          strand <- Refseq_Genes[tx_id %in% unique(Sample_File$transcript), strand][1]
          full_gene_figure(testgenes, unique(Sample_File$transcript), sig_introns, strand)
          # Close the graphics device
          dev.off()



          normaltable <- report_table[two_sd == TRUE & intron_no %in% sig_introns & SJ_IR == "SJ" & annotated == "canonical"]

          normalevents <- normaltable[, .(seqnames, event, difference, probandpct, controlavg, frame_conserved)][order(abs(difference), decreasing = T)]

          names(normalevents) <- c("splice-junction", "event", "difference", "proband", "controls", "frame")

          normalChangeBarPlot(normaltable)

          splicingFrameConsequences(normaltable, sig_introns[[1]])

          sig_introns_list <- list()

          for (intron in sig_introns) {
            message("## Intron ", intron)

            introntable <- report_table[two_sd == TRUE & introns == intron]
            insigrow <- colSums(table[two_sd == FALSE & introns == intron, .(difference, probandpct, controlavg)])
            insigrow <- data.table::as.data.table(t(c("", "NS events < 2sd", insigrow, "")))
            insigrow$difference <- formattable::percent(insigrow$difference, digits = 1)
            insigrow$probandpct <- formattable::percent(insigrow$probandpct, digits = 1)
            insigrow$controlavg <- formattable::percent(insigrow$controlavg, digits = 1)


            events <- introntable[, .(seqnames, event, difference, probandpct, controlavg, frame_conserved)][order(abs(difference), decreasing = T)]

            names(events) <- c("splice-junction", "event", "difference", "proband", "controls", "frame")

            sig_introns_list[[intron]] <- rbind(events, insigrow, use.names = FALSE)


            sds <- table[, lapply(.(two_sd, three_sd, four_sd), sum)]


            rmarkdown::render(
              input = "R/splicing_diagnostics_report.Rmd",
              output_file = paste(Export, "/", Sample_File$sampleID[sample_number], "_", testgenes, "_combined_full",
                ".html",
                sep = ""
              ),
              params = list(
                "normal" = normalevents,
                "aberrant" = sig_introns_list,
                "testgenes" = testgenes,
                "control_no" = length(ctrls),
                "sample" = Sample_File$sampleID[[sample_number]],
                "transcript" = Sample_File$transcript[[sample_number]],
                "refseq" = Refseq_Genes,
                "strand" = strand,
                "sds" = sds,
                "details" = "../../../Reports/v0.4 Cortar 05.22/AGRF Batch 3/Blood/AGRF_blood_batch_3_details.txt"
              )
            )
          }

          # Full dataset export - no report option
        } else if (full_all_genes_report == T) {
          # Exporting the combineddt dataframe to an excel spreadsheet
          openxlsx::write.xlsx(all_splicing_events_sample,
            paste(Export, Sample_File$sampleID[sample_number],
              "_combined_dt_", Sample_File$assembly[1],
              ".xlsx",
              sep = ""
            ),
            asTable = T,
            overwrite = T
          )
        }
      }
    }
  }
}

normalSpliceMap <- function(table, familycols, proband, genes, export){

    # table <- all_splicing_events_sample
    # probands <- familycols
    # genes <- testgenes

    #png(paste0("output/figs/",proband,"_normalSpliceMap.png"), width = 10, height = 6, units="cm", res = 150)
    pdf(paste0(export, proband,"_normalSpliceMap.pdf"), width = 5, height = 3)

    filtered_table <- as.data.table(table[SJ_IR == "SJ" & annotated == 'canonical' & gene == genes])
    print(filtered_table)

    probpct <- as.vector(filtered_table[, ..familycols])

    myplot1 <- ggplot2::ggplot() +
        geom_ribbon(aes(ymin = filtered_table$controlavg - filtered_table$controlsd*2,
                        ymax = filtered_table$controlavg + filtered_table$controlsd*2,
                        x = filtered_table$intron_no), fill = "grey70", color = "grey70") +
        #geom_col(aes(x=filtered_table$intron_no, y=abs(filtered_table$difference)/filtered_table$controlavg), fill = "maroon") +
        #geom_col(aes(x=filtered_table$intron_no, y=abs(filtered_table$difference)), fill = "black") +
        #geom_line(aes(x=filtered_table$intron_no, y=filtered_table$controlsd*2), color = "black") +
        geom_point(aes(x=filtered_table$intron_no, y=probpct), color = "red") +
        geom_point(aes(x=filtered_table$intron_no, y=filtered_table$controlavg), color = "blue") +
        geom_line(aes(x=filtered_table$intron_no, y=probpct), color = "red") +
        geom_line(aes(x=filtered_table$intron_no, y=filtered_table$controlavg), color = "blue") +
        scale_y_continuous(breaks=seq(0,1.0,0.1), limits = c(-0.5,1.5)) +
        #scale_x_continuous(breaks=seq(1,23,1), limits = c(1,23)) +
        ggtitle(paste0(proband)) + xlab("intron") + ylab("proportion of all splicing") +
        theme_minimal()

    print(myplot1)

    dev.off()

    #png(paste0("output/figs/",proband,"_normalSpliceMap_bar.png"), width = 10, height = 6, units="cm", res = 150)
    pdf(paste0(export, proband,"_normalSpliceMap_bar.pdf"), width = 5, height = 3)

    filtered_table <- as.data.table(table[SJ_IR == "SJ" & annotated == 'canonical' & gene == genes])

    probpct <- as.vector(filtered_table[, ..familycols])

    myplot1 <- ggplot() +
        geom_ribbon(aes(ymin = filtered_table$controlsd*-2,
                        ymax = filtered_table$controlsd*2,
                        x = filtered_table$intron_no), fill = "grey70", color = "grey70") +
        geom_col(aes(x=filtered_table$intron_no, y=filtered_table$difference/filtered_table$controlavg), fill = "maroon", width = 0.5) +
        geom_col(aes(x=filtered_table$intron_no, y=filtered_table$difference), fill = "black", width = 0.5) +
        #geom_line(aes(x=filtered_table$intron_no, y=filtered_table$controlsd*2), color = "black") +
        #geom_line(aes(x=filtered_table$intron_no, y=filtered_table$controlsd*-2), color = "black") +
        geom_hline(aes(yintercept=0.5),linetype=2, color = "grey70") +
        geom_hline(aes(yintercept=0)) +
        scale_y_continuous(breaks=seq(-1.0,1.0,0.1), limits = c(-1.5,round(max(filtered_table$difference,filtered_table$difference/filtered_table$controlavg,filtered_table$controlsd*2),digits=1)+0.05)) +
        #scale_x_continuous(breaks=seq(1,23,1), limits = c(1,23)) +
        ggtitle(paste0(proband)) + xlab("intron") + ylab("change in normal splicing prop.") +
        theme_minimal()

    print(myplot1)

    dev.off()
}

#Generate Report
generate.excel <- function(data, familymembers, gene, export, sample){

    # Create an Excel workbook object and add a worksheet
    wb <- openxlsx::createWorkbook()
    sht <- openxlsx::addWorksheet(wb, gene)

    # Create a percent style
    pct <- openxlsx::createStyle(numFmt='PERCENTAGE')
    twodp <- openxlsx::createStyle(numFmt='0.00')
    centre <- openxlsx::createStyle(halign = "center")
    headerStyle <- openxlsx::createStyle(textRotation = 45, fgFill = "#4F81BD",
                               textDecoration = "bold", fontColour = "white")

    # Add data to the worksheet
    openxlsx::writeDataTable(wb, sht, data)

    # Add the percent, centre, and header styles to the desired cells
    openxlsx::addStyle(wb, sht, style=pct, cols = c(13:(16+familymembers-1)),
             rows = 2:(nrow(data)+1), gridExpand=TRUE)
    #note this may change with additional family members
    openxlsx::addStyle(wb, sht, style=centre, cols = c(1,6:8,10:12),
             rows = 2:(nrow(data)+1), gridExpand=TRUE)
    openxlsx::addStyle(wb, sht, headerStyle, cols = 1:(ncol(data)+1),
             rows = 1, gridExpand = TRUE)
    openxlsx::addStyle(wb, sht, twodp, cols = 22, rows = 2:(nrow(data)+1),
             gridExpand = TRUE)

    # Set column widths for event and proband - Currently not working 20211013
    #width_vec <- apply(data, 2, function(x) max(nchar(as.character(x)) + 2,
    #na.rm = TRUE))
    #setColWidths(wb, sht, cols = c(9,12), widths = width_vec[c(9,12)])

    # Add conditionalFormatting to difference and percentage columns
    openxlsx::conditionalFormatting(wb, sht, cols=c(14:(16+familymembers-1)),
                          rows = 2:(nrow(data)+1),
                          rule = NULL, style = c("#FCFCFF","#63BE7B"),
                          type = "colourScale")

    openxlsx::conditionalFormatting(wb, sht, cols=c(13), rows = 2:(nrow(data)+1),
                          rule = NULL, style=c("#F8696B","#FCFCFF","#F8696B"),
                          type = "colourScale")

    # Export report
    openxlsx::saveWorkbook(wb, paste(export,"/",sample,"_",gene,"_combined_dt_",
                           ".xlsx", sep=""),
                 overwrite = T)

}



