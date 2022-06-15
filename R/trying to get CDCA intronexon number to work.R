11791000 %in% seq(11790899,11791206)
Sample_Refseq_Genes_IE <- Refseq_Genes[tx_id %in% Transcripts &
                                        gene_name %in% Genes$`Gene name`&
                                        region_type %in% c('intron','exon')]

testdt[genes == "MTHFR"]

testdt[, event := mapply(gen.exon.range, strand,
                             exon_range_start, exon_range_end, start, end)]

s

for(i in seq(1,nrow(Sample_Refseq_Genes_IE))){
    if(start %in% seq(Sample_Refseq_Genes_IE[i,region_start],Sample_Refseq_Genes_IE[i,region_end]) == TRUE){
        cat(Sample_Refseq_Genes_IE[i,c(region_type, region_no)])
    }
}

gen.exon.range <- function(strand, exon_range_start, exon_range_end, start, end) {

    exon_range_start <- as.numeric(strsplit(exon_range_start,split=" ")[[1]][3])
    exon_range_end <- as.numeric(strsplit(exon_range_end,split=" ")[[1]][3])
    exon_ranges <- c(exon_range_start, exon_range_end)

    if (!is.na(exon_range_start) & !is.na(exon_range_end)) {
        if(max(exon_ranges)-min(exon_ranges) == 1){
            return(paste("exon ",min(exon_ranges)," ~ ","exon ",
                         max(exon_ranges)+1, sep="", collapse=""))
            cat(paste("exon ",min(exon_ranges)," ~ ","exon ",
                      max(exon_ranges)+1, sep="", collapse=""))
        }else{
            return(paste("exon ",min(exon_ranges)," ~ ","exon ",
                         max(exon_ranges)+1, sep="", collapse=""))
            cat(paste("exon ",min(exon_ranges)," ~ ","exon ",
                      max(exon_ranges)+1, sep="", collapse=""))
        }
    }else if (!is.na(exon_range_start)){
        for(i in seq(1,nrow(Sample_Refseq_Genes_IE))){
            if(end %in% seq(Sample_Refseq_Genes_IE[i,region_start],Sample_Refseq_Genes_IE[i,region_end]) == TRUE){
                cryppos <- paste(Sample_Refseq_Genes_IE[i,c(region_type, region_no)])
            }else{cryppos <- "error"}
        }
        if (strand == "-"){
            return(paste(cryppos, "CD", " ~ ", "exon ",
                         exon_range_start+1, sep="", collapse=""))
        }
        if (strand == "+"){
            return(paste("exon ",exon_range_start, " ~ ",
                         cryppos, "CA", sep="", collapse=""))
        }


    }else if (!is.na(exon_range_end)){
        for(i in seq(1,nrow(Sample_Refseq_Genes_IE))){
            if(start %in% seq(Sample_Refseq_Genes_IE[i,region_start],Sample_Refseq_Genes_IE[i,region_end]) == TRUE){
                cryppos <- paste(Sample_Refseq_Genes_IE[i,c(region_type, region_no)])
            }else{cryppos <- "error"}
        }
        if (strand == "-"){
            return(paste("exon ", exon_range_end, " ~ ",
                         cryppos, "CA", sep="", collapse=""))
        }
        if (strand == "+"){
            return(paste(cryppos, "CD", " ~ ", "exon ",
                         exon_range_end+1, sep="", collapse=""))
        }
    }else{
        return("unannotated junctions")
    }
}




