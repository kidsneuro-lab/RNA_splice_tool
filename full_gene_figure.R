full_gene_figure <- function(gene,tx,sig_introns,strand){

    refseq <- fread("refseq_introns_exons_hg38.tsv", sep = "\t")

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
    png("full_gene_figure.png", width = 20, height = 2.5, units="cm", res = 150)

    par(mar =c(1,1,1,1))
    plot.new()
    height1 = 0.2
    height2 = 0.09
    height3 = 0.28

    if(strand == "+"){
        arrows(x0 = 0,y0 = 0.1,x1 = 0.05,y1 = 0.1, length=0.07, lwd = 2)
    }else{
        arrows(x0 = 1,y0 = 0.1,x1 = 0.95,y1 = 0.1, length=0.07, lwd = 2)
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
        polygon(c(start, start, end, end), c((0.5-height3),(0.5+height3),(0.5+height3),(0.5-height3)), lwd=2, border="red")
    }

}
