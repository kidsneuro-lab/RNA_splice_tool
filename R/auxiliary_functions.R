
gen.exon.range <- function(strand, exon_range_start, exon_range_end) {

    exon_range_start <- as.numeric(strsplit(exon_range_start,split=" ")[[1]][3])
    exon_range_end <- as.numeric(strsplit(exon_range_end,split=" ")[[1]][3])
    exon_ranges <- c(exon_range_start, exon_range_end)

    if (!is.na(exon_range_start) & !is.na(exon_range_end)) {
        if(max(exon_ranges)-min(exon_ranges) == 1){
            return(paste("exon ",min(exon_ranges)," ~ ","exon ",
                         max(exon_ranges)+1, sep="", collapse=""))
        }else{
            return(paste("exon ",min(exon_ranges)," ~ ","exon ",
                         max(exon_ranges)+1, sep="", collapse=""))
        }
    }else if (!is.na(exon_range_start)){
        if (strand == "-"){
            return(paste("cryptic donor", " ~ ", "exon ",
                         exon_range_start+1, sep="", collapse=""))
        }
        if (strand == "+"){
            return(paste("exon ",exon_range_start, " ~ ",
                         "cryptic acceptor", sep="", collapse=""))
        }


    }else if (!is.na(exon_range_end)){
        if (strand == "-"){
            return(paste("exon ", exon_range_end, " ~ ",
                         "cryptic acceptor", sep="", collapse=""))
        }
        if (strand == "+"){
            return(paste("cryptic donor", " ~ ", "exon ",
                         exon_range_end+1, sep="", collapse=""))
        }
    }else{
        return("unannotated junctions")
    }
}

gen.normal <- function(exon_range_start, exon_range_end) {

    exon_range_start <- as.numeric(strsplit(exon_range_start,split=" ")[[1]][3])
    exon_range_end <- as.numeric(strsplit(exon_range_end,split=" ")[[1]][3])
    exon_ranges <- c(exon_range_start, exon_range_end)

    if (!is.na(exon_range_start) & !is.na(exon_range_end)) {
        if(max(exon_ranges)-min(exon_ranges) == 0){
            return(paste("Y"))
        }else{
            return("")
        }
    }else{
        return("")
    }
}

#Introns
gen.introns <- function(exon_range_start, exon_range_end) {

    exon_range_start <-as.numeric(strsplit(exon_range_start,split=" ")[[1]][3])
    exon_range_end <- as.numeric(strsplit(exon_range_end,split=" ")[[1]][3])
    exon_ranges <- c(exon_range_start, exon_range_end)

    if (!is.na(exon_range_start) & !is.na(exon_range_end)) {
        return(paste("intron ",min(exon_ranges, na.rm=T),"-intron ",
                     max(exon_ranges, na.rm=T),sep="", collapse=""))

    }else if (!is.na(exon_range_start) | !is.na(exon_range_end)) {
        return(paste("intron ",min(exon_ranges,na.rm=T), sep="", collapse=""))

    }else{
        return("")
    }
}

#Add gene names
gen.fetch <- function(exon_range_start, exon_range_end) {

    exon_range_start <- strsplit(exon_range_start,split=" ")[[1]][1]
    exon_range_end <- strsplit(exon_range_end,split=" ")[[1]][1]

    if (!is.na(exon_range_start) & !is.na(exon_range_end)) {
        return(exon_range_start)
    }else if (!is.na(exon_range_start)){
        return(exon_range_start)
    }else if (!is.na(exon_range_end)){
        return(exon_range_end)
    }else{
        return("no gene")
    }
}

#Calculate event frame
framed <- function(event, x,y){
    rfsq <- Refseq_Genes
    frame <- c()

    for(i in seq(1,length(x))){
        rg_event <- event[i]
        rg_start <- x[i]
        rg_end <- y[i]
        if(event[i] == "SJ"){
            #is the start/end annotated
            #both annotated
            if(rg_start %in% rfsq$region_start & rg_end %in% rfsq$region_end){
                pairstart <-rfsq$region_start[which(rfsq$region_end == rg_end)]
                pairend <-rfsq$region_end[which(rfsq$region_start == rg_start)]
                if(rg_end == pairend){
                    frame[i] <- TRUE
                } else{
                    dist2authentic <- abs(rg_end-pairend)
                    frame[i] <- dist2authentic%%3 == 0
                }
            }
            #start annotated
            else if(rg_start %in% rfsq$region_start){
                pairend <-rfsq$region_end[which(rfsq$region_start == rg_start)]
                dist2authentic <- abs(rg_end-pairend)
                frame[i] <- dist2authentic%%3 == 0
            }
            #end annotated
            else if(rg_end %in% rfsq$region_end){
                pairstart <-rfsq$region_start[which(rfsq$region_end == rg_end)]
                dist2authentic <- abs(rg_start-pairstart)
                frame[i] <- dist2authentic%%3 == 0
            }
            #unannotated junctions
            else{
                frame[i] <- NA
            }

        }else{
            frame[i] <- ""
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
    addStyle(wb, sht, style=pct, cols = c(14:(17+familymembers-1)),
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
    conditionalFormatting(wb, sht, cols=c(16:(18+familymembers-1)),
                          rows = 2:(nrow(data)+1),
                          rule = NULL, style = c("#FCFCFF","#63BE7B"),
                          type = "colourScale")

    conditionalFormatting(wb, sht, cols=c(15), rows = 2:(nrow(data)+1),
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
