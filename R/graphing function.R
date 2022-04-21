
rnaseqgrapher <- function(combined_dt_intron_test, testgenes, sample, exportlocation){
gene_graph_data <- combined_dt_intron_test[which(combined_dt_intron_test$genes == testgenes)]
#mthfr_values <- fread("../../Reports/New Feature Development/281_MTHFR_P_MTHFR_combined_full.tsv")

gene_graph_data <- gene_graph_data[!which(gene_graph_data$SJ_IR == "IR" & gene_graph_data$annotated == "N")]
gene_graph_data <- gene_graph_data[!which(gene_graph_data$event == "unannotated junctions")]
gene_graph_data[normal == "NULL", normal := ""]
#mthfr_values <- mthfr_values[!which(mthfr_values$SJ_IR == "IR" & mthfr_values$annotated == "N")]
#mthfr_values <- mthfr_values[!which(mthfr_values$event == "unannotated junctions")]
#mthfr_values[normal == "NULL", normal := ""]

sig_introns <- gene_graph_data[,.(sum(two_sd)),by=introns][which(V1 != 0), introns]
#sig_introns <- gene_graph_data[,.(sum(two_sd)),by=introns][, introns]
#sig_introns <- mthfr_values[,.(sum(two_sd)),by=introns][which(V1 != 0), introns]

pctsample <- "probandpct"

gene_graph_data <- gene_graph_data[,c("introns","event", "normal", "SJ_IR", "frame_conserved", pctsample, "controlavg", "difference","two_sd", "unique"), with = F]
#mthfr_graph_data <- mthfr_values[,.(introns,event, normal, SJ_IR, frame_conserved, sj_pct_281_MTHFR_P, controlavg, difference,two_sd, unique)]

gene_graph_data[SJ_IR == "IR", ':=' (normal = "", frame_conserved = "")]
gene_graph_data[, orderintrons := sapply(gene_graph_data$introns,function(x) unlist(strsplit(x," "))[2])]
#mthfr_graph_data[SJ_IR == "IR", ':=' (normal = "", frame_conserved = "")]

gene_graph_data[frame_conserved == FALSE, colour := "red"]
gene_graph_data[frame_conserved == TRUE, colour := "yellow"]
gene_graph_data[SJ_IR == "IR", colour := "black"]
gene_graph_data[normal == "Y", colour := "lightblue"]

#mthfr_graph_data[frame_conserved == FALSE, colour := "red"]
#mthfr_graph_data[frame_conserved == TRUE, colour := "yellow"]
#mthfr_graph_data[SJ_IR == "IR", colour := "black"]
#mthfr_graph_data[normal == "Y", colour := "lightblue"]

gene_graph_data <- gene_graph_data[introns %in% sig_introns][order(as.numeric(orderintrons),normal,frame_conserved, decreasing = F)]
#mthfr_graph_data <- mthfr_graph_data[introns %in% sig_introns][order(introns,normal,frame_conserved, decreasing = F)]

#pdf(paste0(exportlocation,"/",sample,"_",testgenes,"_pie_charts",".pdf"))

if(length(sig_introns) == 0){
    print("no significant changes")
}

par(mfcol = c(2, length(sig_introns)))
par(mar = c(0, 0, 0, 0))

for(i in unique(gene_graph_data$introns)){
    gene_graph_data_working <- gene_graph_data[introns == i]
    pie(gene_graph_data[introns == i, ..pctsample][[1]], labels = "", col = gene_graph_data[introns == i, colour])
    text(c(0, 0), c(1.0,1.0), labels = strsplit(i,split=" ")[[1]][2], col="black", cex = 1.5)
    pie(gene_graph_data[introns == i, controlavg], labels = "", col = gene_graph_data[introns == i, colour])
    text(c(0, 0), c(1.0,1.0), labels = strsplit(i,split=" ")[[1]][2], col="black", cex = 1.5)
}

#dev.off()
}




mthfr <- refseq_introns_exons[gene_name == "MTHFR"]

mthfr <- mthfr[, figure_pos := ifelse(tx_id=="NM_001330358",0,1)]

mthfr_length <- mthfr[,.(max(tx_len))]
mthfr <- mthfr[, ':=' (region_start = (region_start-min(region_start))/max(tx_len),
                       region_end = (region_end-min(region_start))/max(tx_len))]
coding_mthfr <- mthfr[region_type=="cds"]
non_coding_mthfr <- mthfr[region_type %in% c("5utr","3utr")]
labels_mthfr <- mthfr[region_type == "intron"]


#View(mthfr)
par(mfrow =c(1,1))
plot.new()
plot.window(c(0,1.3), c(0,1))
#segments(0.01, 0.5, 0.99, 0.5, lwd=5)
#segments(0.01, 0.3, 0.99, 0.3, lwd=5)
gap = 0.15
height1 = 0.05
height2 = 0.02


for(i in seq(1,nrow(coding_mthfr))){
    start <- coding_mthfr$region_start[i]
    end <- coding_mthfr$region_end[i]
    cat(start,end,"\n")
    polygon(c(start, start, end, end), c((0.5-height1)-(gap*coding_mthfr$figure_pos[i]), (0.5+height1)-(gap*coding_mthfr$figure_pos[i]), (0.5+height1)-(gap*coding_mthfr$figure_pos[i]), (0.5-height1)-(gap*coding_mthfr$figure_pos[i])), lwd=2, col="black")
}

for(i in seq(1,nrow(non_coding_mthfr))){
    start <- non_coding_mthfr$region_start[i]
    end <- non_coding_mthfr$region_end[i]
    cat(start,end,"\n")
    polygon(c(start, start, end, end), c((0.5-height2)-(gap*non_coding_mthfr$figure_pos[i]), (0.5+height2)-(gap*non_coding_mthfr$figure_pos[i]), (0.5+height2)-(gap*non_coding_mthfr$figure_pos[i]), (0.5-height2)-(gap*non_coding_mthfr$figure_pos[i])), lwd=2,col="black")
}

for(i in seq(1,nrow(labels_mthfr))){
    start_label <- labels_mthfr[canonical == 1, region_start][i]
    end_label <- labels_mthfr[canonical == 1, region_end][i]
    start <- labels_mthfr[, region_start][i]
    end <- labels_mthfr[, region_end][i]
    label_pos <- (start_label + end_label)/2
    pos <- (start + end)/2
    cat(label_pos,"\n")
    text(c(label_pos, label_pos), c(0.60, 0.60), labels = labels_mthfr[canonical == 1, region_no][i], col="black")
    segments(pos-0.002, 0.5-(gap*labels_mthfr$figure_pos[i]), pos+0.002, 0.51-(gap*labels_mthfr$figure_pos[i]), lwd=3)
    segments(pos-0.002, 0.5-(gap*labels_mthfr$figure_pos[i]), pos+0.002, 0.49-(gap*labels_mthfr$figure_pos[i]), lwd=3)
}

for(i in seq(1,nrow(labels_mthfr))){
    start <- labels_mthfr$region_start[i]
    end <- labels_mthfr$region_end[i]
    segments(start, 0.5-(gap*labels_mthfr$figure_pos[i]), end, 0.5-(gap*labels_mthfr$figure_pos[i]), lwd=3)
}

for(i in seq(1,length(unique(mthfr$tx_id)))){
    pos <- i-1
    text(c(1.2, 1.2), c(0.50-(gap*pos), 0.50-(gap*pos)), labels = unique(mthfr$tx_id)[i], col="black", adj = c(1,0.5))
}

