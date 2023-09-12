generateReport <- function(comparisons, Sample_File, Export, mode, prefix) {
    message("Generating reports...")

  if(mode == "default" | mode == "panel"){
    for (sample_number in seq(1, nrow(Sample_File))) {
      # Initialise new copy of the all_splicing_events dataset
      if (Sample_File$sampletype[sample_number] == "test") {
        all_splicing_events_sample <- comparisons[[sample_number]]
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
        full_all_genes_report <- T
        figure <- F

        if(mode == "panel"){
          testgenes <- unique(all_splicing_events_sample$gene)
          testgenename <- "panel"
        }else{
          testgenes <- unique(Sample_File$genes[sample_number])
          testgenename <- unique(Sample_File$genes[sample_number])
        }
        # will need a for loop
        if(figure == T){
          for(gene in testgenes){
            normalSpliceMap(all_splicing_events_sample, familycols[1], proband, gene, export = Export, mode = mode, prefix = prefix)
          }
        }
        # Generate filtered excel spreadsheet +/- summary html
        # Excel spreadsheet
        if (report == TRUE) {
          generate.excel(all_splicing_events_sample[gene %in% testgenes],
            length(familycols),
            gene = testgenename,
            Export,
            Sample_File$sampleID[sample_number],
            prefix = prefix
          )
        }
        if (full_all_genes_report == TRUE) {
          data.table::fwrite(all_splicing_events_sample,
            file = paste0(
              Export, "/",prefix,Sample_File$sampleID[sample_number], "_", testgenename, "_combined_full",
              ".tsv"
            ), sep = "\t"
          )
        }
      }
    }
   }else if(mode == "research"){
    all_splicing_events_sample <- comparisons
    testgenes <- unique(all_splicing_events_sample$gene)
    proband <- "splicing_analysis"
    for(gene in testgenes){
      normalSpliceMap(all_splicing_events_sample,
                      familycols[1],
                      proband,
                      gene,
                      export = Export,
                      mode = mode,
                      prefix = prefix)
    }
    data.table::fwrite(all_splicing_events_sample,
                       file = paste0(
                         Export, "/",prefix,"splicing_analysis","_combined_full",
                         ".tsv"
                       ), sep = "\t"
    )
  }
}

normalSpliceMap <- function(table, familycols, proband, genes, export, mode, prefix){

  pdf(paste0(export,"/",prefix,proband,"_",genes,"_normalSpliceMap.pdf"), width = 5, height = 3)

    filtered_table <- as.data.table(table[SJ_IR == "SJ" & annotated == 'canonical' & gene == genes])

    if(mode == "default" | mode == "panel"){
      probpct <- unlist(as.vector(filtered_table[, ..familycols]))
      probcolour <- "red"
    }else if (mode == "research"){
      probpct <- c(0)
      filtered_table$difference <- c(0)
      probcolour <- "transparent"
    }

    myplot1 <- ggplot2::ggplot() +
        geom_ribbon(aes(ymin = filtered_table$controlavg - filtered_table$controlsd*2,
                        ymax = filtered_table$controlavg + filtered_table$controlsd*2,
                        x = filtered_table$intron_no), fill = "grey70", color = "grey70") +
        geom_point(aes(x=filtered_table$intron_no, y=probpct), color = probcolour) +
        geom_point(aes(x=filtered_table$intron_no, y=filtered_table$controlavg), color = "blue") +
        geom_line(aes(x=filtered_table$intron_no, y=probpct), color = probcolour) +
        geom_line(aes(x=filtered_table$intron_no, y=filtered_table$controlavg), color = "blue") +
        scale_y_continuous(breaks=seq(0,1.0,0.1), limits = c(-0.5,1.5)) +
        ggtitle(paste0(proband,"_",genes)) + xlab("intron") + ylab("proportion of all splicing") +
        theme_minimal()

    print(myplot1)

    dev.off()

    pdf(paste0(export,"/",prefix,proband,"_",genes,"_normalSpliceMap_bar.pdf"), width = 5, height = 3)

    myplot1 <- ggplot() +
        geom_ribbon(aes(ymin = filtered_table$controlsd*-2,
                        ymax = filtered_table$controlsd*2,
                        x = filtered_table$intron_no), fill = "grey70", color = "grey70") +
        geom_col(aes(x=filtered_table$intron_no, y=filtered_table$difference/filtered_table$controlavg), fill = "maroon", width = 0.5) +
        geom_col(aes(x=filtered_table$intron_no, y=filtered_table$difference), fill = "black", width = 0.5) +
        geom_hline(aes(yintercept=0.5),linetype=2, color = "grey70") +
        geom_hline(aes(yintercept=0)) +
        scale_y_continuous(breaks=seq(-1.0,1.0,0.1), limits = c(-1.5,round(max(filtered_table$difference,filtered_table$difference/filtered_table$controlavg,filtered_table$controlsd*2),digits=1)+0.05)) +
        ggtitle(paste0(proband,"_",genes)) + xlab("intron") + ylab("change in normal splicing prop.") +
        theme_minimal()

    print(myplot1)

    dev.off()
}

#Generate Report
generate.excel <- function(data, familymembers, gene, export, sample, prefix){

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

    # Add conditionalFormatting to difference and percentage columns
    openxlsx::conditionalFormatting(wb, sht, cols=c(14:(16+familymembers-1)),
                          rows = 2:(nrow(data)+1),
                          rule = NULL, style = c("#FCFCFF","#63BE7B"),
                          type = "colourScale")

    openxlsx::conditionalFormatting(wb, sht, cols=c(13), rows = 2:(nrow(data)+1),
                          rule = NULL, style=c("#F8696B","#FCFCFF","#F8696B"),
                          type = "colourScale")

    # Export report
    openxlsx::saveWorkbook(wb, paste(export,"/",prefix,sample,"_",gene,"_combined_dt_",
                           ".xlsx", sep=""),
                 overwrite = T)

}
