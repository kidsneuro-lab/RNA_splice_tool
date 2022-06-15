library(data.table)
library(seqinr)
library(stringr)

irfasta <- read.fasta(file = "refseq_genes_introns.fasta", seqtype = c("DNA"), as.string = TRUE)
irfasta.dt <- t(as.data.table(irfasta))
irfasta.dt <- cbind(irfasta.dt, rownames(irfasta.dt))
rownames(irfasta.dt) <- NULL
irfasta.dt <- as.data.table(irfasta.dt)

irfasta.dt_results <- as.data.frame(matrix(nrow=0,ncol=14))

transint <- c()
transcrpt <- c()
intron <- c()
chrpos <- c()
chr <- c()
start <- c()
end <- c()
window1ptc_pos <- c()
window1ptc <- c()
window1_no_ptcs <- c()
window1_all_ptcs <- c()
window2ptc_pos <- c()
window2ptc <- c()
window2_no_ptcs <- c()
window2_all_ptcs <- c()
window3ptc_pos <- c()
window3ptc <- c()
window3_no_ptcs <- c()
window3_all_ptcs <- c()

for(i in seq(1,nrow(irfasta.dt))){
    transint[i] <- strsplit(as.character(irfasta.dt[i,2]),"::")[[1]][1]
    transcrpt[i] <- paste0(strsplit(transint[i],"_")[[1]][1],"_",strsplit(transint,"_")[[1]][2])
    intron[i] <- as.numeric(strsplit(transint[i],"_")[[1]][4])+1

    chrpos[i] <- strsplit(as.character(irfasta.dt[i,2]),"::")[[1]][2]
    chr[i] <- strsplit(chrpos[i],":")[[1]][1]
    start[i] <- strsplit(strsplit(chrpos[i],":")[[1]][2],"-")[[1]][1]
    end[i] <- strsplit(strsplit(chrpos[i],":")[[1]][2],"-")[[1]][2]

    window1ptc_pos <- which(substring(irfasta.dt[i,1],seq(1,str_length(irfasta.dt[i,1])-3,3),seq(3,str_length(irfasta.dt[i,1]),3)) %in% c("tag","taa","tga"))*3
    window1ptc[i] <- min(window1ptc_pos)
    window1_no_ptcs[i] <- length(window1ptc_pos)
    window1_all_ptcs[i] <- paste0(window1ptc_pos,collapse=",")

    window2ptc_pos <- which(substring(irfasta.dt[i,1],seq(2,str_length(irfasta.dt[i,1])-3,3),seq(4,str_length(irfasta.dt[i,1]),3)) %in% c("tag","taa","tga"))*3
    window2ptc[i] <- min(window2ptc_pos)
    window2_no_ptcs[i] <- length(window2ptc_pos)
    window2_all_ptcs[i] <- paste0(window2ptc_pos,collapse=",")

    window3ptc_pos <- which(substring(irfasta.dt[i,1],seq(3,str_length(irfasta.dt[i,1])-3,3),seq(5,str_length(irfasta.dt[i,1]),3)) %in% c("tag","taa","tga"))*3
    window3ptc[i] <- min(window3ptc_pos)
    window3_no_ptcs[i] <- length(window3ptc_pos)
    window3_all_ptcs[i] <- paste0(window3ptc_pos,collapse=",")
}

irfasta.dt_results <- as.data.table(list(chr,start,end,transcrpt,intron,window1ptc,window2ptc,window3ptc,window1_no_ptcs,window2_no_ptcs,window3_no_ptcs,window1_all_ptcs,window2_all_ptcs,window3_all_ptcs))

names(irfasta.dt_results) <- c("chr","start","end","transcript","intron","window1ptc","window2ptc","window3ptc","window1_no_ptcs","window2_no_ptcs","window3_no_ptcs","window1_all_ptcs","window2_all_ptcs","window3_all_ptcs")

window1_irfasta.dt_results <- irfasta.dt_results[,c(1,2,3,4,5,6,9,12)]

fwrite(irfasta.dt_results,"refseq_genes_introns_ptcs_all_windows.tsv",sep = "\t")
fwrite(window1_irfasta.dt_results,"refseq_genes_introns_ptcs_window_1.tsv",sep = "\t")

