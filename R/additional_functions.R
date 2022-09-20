#Takes a data-table with events and returns them in a human readable format
eventAnnotation <- function(intron.dt){
    event <- 1
    events <- c()
    for(event in seq(1,nrow(intron.dt))){
        #Extract intron start and end ranges and strand
        exon_range_start <- query_intron.dt[event,intron_jxn_start]
        exon_range_end <- query_intron.dt[event,intron_jxn_end]
        exon_ranges <- c(exon_range_start, exon_range_end)
        strand <- as.character(intron.dt$strand[event])

        #normal splicing
        if (!is.na(exon_range_start) & !is.na(exon_range_end)) {
            if(max(exon_ranges)-min(exon_ranges) == 0){
                if(intron.dt$SJ_IR[event] == "SJ"){
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
framed <- function(intron.dt){
    rfsq <- Refseq_Genes
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

        }else{
            frame[event] <- ""
        }
    }
    return(frame)
}

#Generate Report
generate.report <- function(data, familymembers, gene, export, sample){

    # Create an Excel workbook object and add a worksheet
    wb <- createWorkbook()
    sht <- addWorksheet(wb, gene)

    # Create a percent style
    pct <- createStyle(numFmt='PERCENTAGE')
    twodp <- createStyle(numFmt='0.00')
    centre <- createStyle(halign = "center")
    headerStyle <- createStyle(textRotation = 45, fgFill = "#4F81BD",
                               textDecoration = "bold", fontColour = "white")

    # Add data to the worksheet
    writeDataTable(wb, sht, data)

    # Add the percent, centre, and header styles to the desired cells
    addStyle(wb, sht, style=pct, cols = c(13:(16+familymembers-1)),
             rows = 2:(nrow(data)+1), gridExpand=TRUE)
    #note this may change with additional family members
    addStyle(wb, sht, style=centre, cols = c(1,6:8,10:12),
             rows = 2:(nrow(data)+1), gridExpand=TRUE)
    addStyle(wb, sht, headerStyle, cols = 1:(ncol(data)+1),
             rows = 1, gridExpand = TRUE)
    addStyle(wb, sht, twodp, cols = 22, rows = 2:(nrow(data)+1),
             gridExpand = TRUE)

    # Set column widths for event and proband - Currently not working 20211013
    #width_vec <- apply(data, 2, function(x) max(nchar(as.character(x)) + 2,
    #na.rm = TRUE))
    #setColWidths(wb, sht, cols = c(9,12), widths = width_vec[c(9,12)])

    # Add conditionalFormatting to difference and percentage columns
    conditionalFormatting(wb, sht, cols=c(14:(16+familymembers-1)),
                          rows = 2:(nrow(data)+1),
                          rule = NULL, style = c("#FCFCFF","#63BE7B"),
                          type = "colourScale")

    conditionalFormatting(wb, sht, cols=c(13), rows = 2:(nrow(data)+1),
                          rule = NULL, style=c("#F8696B","#FCFCFF","#F8696B"),
                          type = "colourScale")

    # Export report
    saveWorkbook(wb, paste(export,"/",sample,"_",gene,"_combined_dt_",
                           ".xlsx", sep=""),
                 overwrite = T)
}

#Simple code for the opposite of %in%
`%nin%` <- Negate(`%in%`)

GRanges.to.SAF <- function(gr, minAnchor=1){
    data.table(
        GeneID  = seq_along(gr),
        Chr     = as.factor(seqnames(gr)),
        Start   = start(gr) - (minAnchor - 1),
        End     = end(gr) + (minAnchor - 1),
        Strand  = as.factor(strand(gr))
    )
}

full_gene_figure <- function(gene,tx,sig_introns,strand){

    refseq <- fread("ref/refseq_introns_exons_hg38.tsv", sep = "\t")

    table <- refseq[gene_name == gene & tx_id %in% tx]

    table_length <- table[,.(max(tx_len))]
    table <- table[, ':=' (region_start = (region_start-min(region_start))/max(tx_len),
                           region_end = (region_end-min(region_start))/max(tx_len))]
    coding_table <- table[region_type=="cds"]
    non_coding_table <- table[region_type %in% c("5utr","3utr")]
    labels_table <- table[region_type == "intron"]
    sig_table <- labels_table[region_no %in% sig_introns]

    #View(mthfr)
    # SVG graphics device
    png("output/figs/full_gene_figure.png", width = 20, height = 3, units="cm", res = 150)

    par(mar =c(1,1,1,1))
    plot.new()
    height1 = 0.18
    height2 = 0.07
    height3 = 0.28

    if(strand == "+"){
        arrows(x0 = 0,y0 = 0.1,x1 = 0.05,y1 = 0.1, length=0.07, lwd = 2)
    }else{
        arrows(x0 = 1,y0 = 0.1,x1 = 0.95,y1 = 0.1, length=0.07, lwd = 2)
    }

    for(i in seq(1,nrow(labels_table))){
        start_label <- labels_table[canonical == 1, region_start][i]
        end_label <- labels_table[canonical == 1, region_end][i]
        start <- labels_table[, region_start][i]
        end <- labels_table[, region_end][i]
        label_pos <- (start_label + end_label)/2
        pos <- (start + end)/2
    }

    for(i in seq(1,nrow(labels_table))){
        start <- labels_table$region_start[i]
        end <- labels_table$region_end[i]
        segments(start, 0.5, end, 0.5, lwd=3)
    }

    for(i in seq(1,nrow(sig_table))){
        start <- sig_table$region_start[i]
        end <- sig_table$region_end[i]
        segments(start, 0.5, end, 0.5, lwd=3, col = "red")
        text(x = (start+end)/2, y = 0.9, sig_table$region_no[i])
    }

    for(i in seq(1,nrow(coding_table))){
        start <- coding_table$region_start[i]
        end <- coding_table$region_end[i]
        polygon(c(start, start, end, end), c((0.5-height1),(0.5+height1),(0.5+height1),(0.5-height1)), lwd=2, col="black")
    }

    for(i in seq(1,nrow(non_coding_table))){
        start <- non_coding_table$region_start[i]
        end <- non_coding_table$region_end[i]
        polygon(c(start, start, end, end), c((0.5-height2),(0.5+height2),(0.5+height2),(0.5-height2)), lwd=2, col="black")
    }
}

normalChangeBarPlot <- function(table){

    png("output/figs/normal_change_barplot.png", width = 6, height = 8, units="cm", res = 150)

    par(mar=c(5,5,2,2))
    plot.new()

    names_prop_ns <- c("controls","proband")
    mean_prop_ns <- c(as.numeric(table[1,controlavg]/table[1,controlavg]), as.numeric(table[1,probandpct]/table[1,controlavg]))
    sd_prop_ns <- c(table[1,controlsd],0)
    barplot(mean_prop_ns ~ names_prop_ns, las = 2, xlab = " ", ylab = " ", col = c("grey","white"), cex.lab = 1.7,
            cex.main = 1.5, cex.names = 0.7, axes = TRUE, ylim = c(0, 1+0.01+sd_prop_ns[[1]]), cex.axis = 1)

    arrows(x0=0.70, y0=1-sd_prop_ns[[1]], x1=0.70, y1=1+sd_prop_ns[[1]], code=3, angle=90, length=0.1, col="black", lwd=1)

    mtext("Proportion of Control Average", side = 2, line = 3, cex = 0.7, font = 2)
    dev.off()
}

splicingFrameConsequences <- function(table,intron){

    png("output/figs/splicing_frame_consequences.png", width = 10, height = 6, units="cm", res = 150)

    par(mar=c(3,4,2,1))

    table_ctrl_sj <- table[intron_no == intron & annotated != 'canonical' & SJ_IR == "SJ",sum(controlavg), by = frame_conserved]
    table_ctrl_ns <- table[intron_no == intron & annotated == 'canonical' & SJ_IR == "SJ",sum(controlavg), by = event]
    table_ctrl_ir <- table[intron_no == intron & annotated == 'canonical' & SJ_IR == "IR",sum(controlavg), by = event]

    table_pbd_sj <- table[intron_no == intron & annotated != 'canonical' & SJ_IR == "SJ",sum(probandpct), by = frame_conserved]
    table_pbd_ns <- table[intron_no == intron & annotated == 'canonical' & SJ_IR == "SJ",sum(probandpct), by = event]
    table_pbd_ir <- table[intron_no == intron & annotated == 'canonical' & SJ_IR == "IR",sum(probandpct), by = event]

    in_frame <- data.table(frame_conserved ="in-frame", V1 = formattable::percent(0, digits=1))

    if(nrow(table_ctrl_sj) == 1){
        table_ctrl_sj <- rbind(table_ctrl_sj,in_frame)
    }

    if(nrow(table_pbd_sj) == 1){
        table_pbd_sj <- rbind(table_pbd_sj,in_frame)
    }

    table_frame <- rbind(table_ctrl_ns ,table_ctrl_sj ,table_ctrl_ir, use.names=F)
    table_frame <- cbind(table_frame,rbind(table_pbd_ns[,2], table_pbd_sj[,2], table_pbd_ir[,2]))

    table_frame <- table_frame[,2:3]
    table_frame_mat <- as.matrix(table_frame)
    names(table_frame_mat) <- c("controlavg","proband")
    row.names(table_frame_mat) <- c("normal-splicing","in-frame","out-of-frame","intron-retention")

    barplot(table_frame_mat ~ names(table_frame_mat),
            col = c("white","red","yellow","black"),
            las = 2,
            legend = FALSE, beside = TRUE, space=c(0,1),
            ylim = c(0, 1))

    mtext("Proportion of Reads at Junction", side = 2, line = 3, cex = 0.8, font = 2)
    rm(table_frame)
    dev.off()
}


normalSpliceMap <- function(table, probands, genes){

    # table <- all_splicing_events_sample
    # probands <- familycols
    # genes <- testgenes

    #png(paste0("output/figs/",proband,"_normalSpliceMap.png"), width = 10, height = 6, units="cm", res = 150)
    pdf(paste0("output/figs/",proband,"_normalSpliceMap.pdf"), width = 5, height = 3)

    filtered_table <- as.data.frame(table[SJ_IR == "SJ" & annotated == 'canonical' & gene == genes])

    probpct <- as.vector(filtered_table[,probands])

    myplot1 <- ggplot() +
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
        scale_y_continuous(breaks=seq(0,1.0,0.1), limits = c(0,1.0)) +
        scale_x_continuous(breaks=seq(1,23,1), limits = c(1,23)) +
        ggtitle(paste0(proband)) + xlab("intron") + ylab("proportion of all splicing") +
        theme_minimal()

    print(myplot1)

    dev.off()

    #png(paste0("output/figs/",proband,"_normalSpliceMap_bar.png"), width = 10, height = 6, units="cm", res = 150)
    pdf(paste0("output/figs/",proband,"_normalSpliceMap_bar.pdf"), width = 5, height = 3)

    filtered_table <- as.data.frame(table[SJ_IR == "SJ" & annotated == 'canonical' & gene == genes])

    probpct <- as.vector(filtered_table[,probands])

    myplot1 <- ggplot() +
        geom_ribbon(aes(ymin = filtered_table$controlsd*-2,
                        ymax = filtered_table$controlsd*2,
                        x = filtered_table$intron_no), fill = "grey70", color = "grey70") +
        geom_col(aes(x=filtered_table$intron_no, y=filtered_table$difference/filtered_table$controlavg), fill = "maroon", width = 0.5) +
        geom_col(aes(x=filtered_table$intron_no, y=filtered_table$difference), fill = "black", width = 0.5) +
        #geom_line(aes(x=filtered_table$intron_no, y=filtered_table$controlsd*2), color = "black") +
        #geom_line(aes(x=filtered_table$intron_no, y=filtered_table$controlsd*-2), color = "black") +
        geom_hline(aes(yintercept=0.355),linetype=2, color = "grey70") +
        geom_hline(aes(yintercept=-0.355),linetype=2, color = "grey70") +
        geom_hline(aes(yintercept=0)) +
        scale_y_continuous(breaks=seq(-1.0,1.0,0.1), limits = c(-1.0,round(max(filtered_table$difference,filtered_table$controlsd*2),digits=1)+0.05)) +
        scale_x_continuous(breaks=seq(1,23,1), limits = c(1,23)) +
        ggtitle(paste0(proband)) + xlab("intron") + ylab("change in normal splicing prop.") +
        theme_minimal()

    print(myplot1)

    dev.off()
}


#
# splicingFrameConsequences <- function(table,intron){
#
    # table <- normaltable
    # intron <- sig_introns[[1]]
#
#     png("output/figs/splicing_frame_consequences.png", width = 10, height = 6, units="cm", res = 150)
#
#     par(mar=c(3,4,2,1))
#
#     table_ctrl_sj <<- table[intron_no == intron & annotated != 'canonical' & SJ_IR == "SJ",sum(controlavg), by = frame_conserved]
#     table_ctrl_ns <<- table[intron_no == intron & annotated == 'canonical' & SJ_IR == "SJ",sum(controlavg), by = event]
#     table_ctrl_ir <<- table[intron_no == intron & annotated == 'canonical' & SJ_IR == "IR",sum(controlavg), by = event]
#
#     table_pbd_sj <<- table[intron_no == intron & annotated != 'canonical' & SJ_IR == "SJ",sum(probandpct), by = frame_conserved]
#     table_pbd_ns <<- table[intron_no == intron & annotated == 'canonical' & SJ_IR == "SJ",sum(probandpct), by = event]
#     table_pbd_ir <<- table[intron_no == intron & annotated == 'canonical' & SJ_IR == "IR",sum(probandpct), by = event]
#
#     # for(count in list(table_ctrl_sj,table_ctrl_ns,table_ctrl_ir,table_pbd_sj,table_pbd_ns,table_pbd_ir)){
#     #     if(nrow(count) == 0){
#     #         count$event[1] <- "no event"
#     #         count$V1 <- formattable::percent(0, digits=1)
#     #     }
#     #     print(count)
#     #     print('\n')
#     # }
#
#
#     in_frame <<- data.table(frame_conserved ="in-frame", V1 = formattable::percent(0, digits=1))
#
#     if(nrow(table_ctrl_sj) == 1){
#         table_ctrl_sj <<- rbind(table_ctrl_sj,in_frame)
#     }
#
#     if(nrow(table_pbd_sj) == 1){
#         table_pbd_sj <<- rbind(table_pbd_sj,in_frame)
#     }
#
#     table_frame <<- rbind(table_ctrl_ns ,table_ctrl_sj ,table_ctrl_ir, use.names=F)
#     table_frame <<- cbind(table_frame,rbind(table_pbd_ns[,2], table_pbd_sj[,2], table_pbd_ir[,2]))
#
#     table_frame <<- table_frame[,2:3]
#     table_frame_mat <<- as.matrix(table_frame)
#     names(table_frame_mat) <- c("controlavg","proband")
#     row.names(table_frame_mat) <<- c("normal-splicing","in-frame","out-of-frame","intron-retention")[1]
#
#     barplot(table_frame_mat ~ names(table_frame_mat),
#             col = c("white","red","yellow","black"),
#             las = 2,
#             legend = FALSE, beside = TRUE, space=c(0,1),
#             ylim = c(0, 1))
#
#     mtext("Proportion of Reads at Junction", side = 2, line = 3, cex = 0.8, font = 2)
#     rm(table_frame)
#     dev.off()
# }













gen.exon.range <- function(strand, exon_range_start, exon_range_end) {
    #
    #     exon_range_start <- as.numeric(strsplit(exon_range_start,split=" ")[[1]][3])
    #     exon_range_end <- as.numeric(strsplit(exon_range_end,split=" ")[[1]][3])
    #     exon_ranges <- c(exon_range_start, exon_range_end)
    #
    #     if (!is.na(exon_range_start) & !is.na(exon_range_end)) {
    #         if(max(exon_ranges)-min(exon_ranges) == 1){
    #             return(paste("exon ",min(exon_ranges)," ~ ","exon ",
    #                          max(exon_ranges)+1, sep="", collapse=""))
    #         }else{
    #             return(paste("exon ",min(exon_ranges)," ~ ","exon ",
    #                          max(exon_ranges)+1, sep="", collapse=""))
    #         }
    #     }else if (!is.na(exon_range_start)){
    #         if (strand == "-"){
    #             return(paste("CD", " ~ ", "exon ",
    #                          exon_range_start+1, sep="", collapse=""))
    #         }
    #         if (strand == "+"){
    #             return(paste("exon ",exon_range_start, " ~ ",
    #                          "CA", sep="", collapse=""))
    #         }
    #
    #
    #     }else if (!is.na(exon_range_end)){
    #         if (strand == "-"){
    #             return(paste("exon ", exon_range_end, " ~ ",
    #                          "CA", sep="", collapse=""))
    #         }
    #         if (strand == "+"){
    #             return(paste("CD", " ~ ", "exon ",
    #                          exon_range_end+1, sep="", collapse=""))
    #         }
    #     }else{
    #         return("unannotated junctions")
    #     }
}
