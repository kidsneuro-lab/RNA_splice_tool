library(data.table)
library(seqinr)
library(stringr)
setwd("C:/Users/rhett/OneDrive/Documents/CAREER/CMRI/Bioinformatics/RNA-seq/ComparativeRNAseqTable/Development/v0.1.0 Cortar 02.22/cortar/R")

Refseq_Genes <- fread("../refseq_introns_exons_hg38.tsv", sep = "\t")
Ensembl_Genes <- fread("../hg38_mart_export_allgenes_chr1-Y.txt", sep = "\t")

Mendeliome_Genes <- fread("IR/clinically_relevant_genes_Sep2021.tsv.txt", sep = "\t")
Mendeliome_Introns <- Refseq_Genes_Introns[which(Refseq_Genes_Introns$gene_name %in% Mendeliome_Genes$Approved.Symbol)]

Refseq_Genes_Introns <- Refseq_Genes[region_type == "intron"]

Mendeliome_Introns$chrom <- paste0("chr",Mendeliome_Introns$chrom)

mend_intr_bed <- Mendeliome_Introns[,.(chrom,
                      region_start-1,
                      region_end,
                      paste0(chrom,"_",region_start,"_",region_end,"_",gene_name,"_",strand),
                      0,
                      strand
                      )]

mend_intr_bed_uniq <- unique(mend_intr_bed)

#bedtools getfasta -fi ref/hg38.fasta -bed all_mendelian_gene_introns.bed -fo all_mendelian_gene_introns_sequences.tsv -name -tab -s

mendelian_intron_sequences <- fread("IR/all_mendelian_gene_introns_sequences.tsv", header = F, sep = "\t")

frame <- c()
intron_length <- c()
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


for(sequence in 1:nrow(mendelian_intron_sequences)){
    intron_length[sequence] <- str_length(mendelian_intron_sequences[[sequence,2]])
    if(intron_length[sequence] < 3){
        message(sequence)
    }

    else{
    frame[sequence] <- if(intron_length %% 3 == 0) TRUE else FALSE

    window1ptc_pos <- which(substring(mendelian_intron_sequences[[sequence,2]],seq(1,str_length(mendelian_intron_sequences[[sequence,2]])-3,3),seq(3,str_length(mendelian_intron_sequences[[sequence,2]]),3)) %in% c("TAG","TAA","TGA"))*3
    window1ptc[sequence] <- min(window1ptc_pos)
    window1_no_ptcs[sequence] <- length(window1ptc_pos)
    window1_all_ptcs[sequence] <- paste0(window1ptc_pos,collapse=",")

    window2ptc_pos <- which(substring(mendelian_intron_sequences[[sequence,2]],seq(2,str_length(mendelian_intron_sequences[[sequence,2]])-3,3),seq(4,str_length(mendelian_intron_sequences[[sequence,2]]),3)) %in% c("TAG","TAA","TGA"))*3
    window2ptc[sequence] <- min(window2ptc_pos)
    window2_no_ptcs[sequence] <- length(window2ptc_pos)
    window2_all_ptcs[sequence] <- paste0(window2ptc_pos,collapse=",")

    window3ptc_pos <- which(substring(mendelian_intron_sequences[[sequence,2]],seq(3,str_length(mendelian_intron_sequences[[sequence,2]])-3,3),seq(5,str_length(mendelian_intron_sequences[[sequence,2]]),3)) %in% c("TAG","TAA","TGA"))*3
    window3ptc[sequence] <- min(window3ptc_pos)
    window3_no_ptcs[sequence] <- length(window3ptc_pos)
    window3_all_ptcs[sequence] <- paste0(window3ptc_pos,collapse=",")
    }
}

















fwrite(mend_intr_bed_uniq, "mend_intron.bed", sep="\t", col.names = F)


irfasta <- read.fasta(file = "IR/cortar_genes.fasta", seqtype = c("DNA"), as.string = TRUE)
mendelian_intron_sequences <- t(as.data.table(irfasta))
mendelian_intron_sequences <- cbind(mendelian_intron_sequences, rownames(mendelian_intron_sequences))
rownames(mendelian_intron_sequences) <- NULL
mendelian_intron_sequences <- as.data.table(mendelian_intron_sequences)

mendelian_intron_sequences_results <- as.data.frame(matrix(nrow=0,ncol=14))



for(i in seq(1,nrow(mendelian_intron_sequences))){
    transint[i] <- strsplit(as.character(mendelian_intron_sequences[i,2]),"::")[[1]][1]
    transcrpt[i] <- paste0(strsplit(transint[i],"_")[[1]][1],"_",strsplit(transint,"_")[[1]][2])
    intron[i] <- as.numeric(strsplit(transint[i],"_")[[1]][4])+1

    chrpos[i] <- strsplit(as.character(mendelian_intron_sequences[i,2]),"::")[[1]][2]
    chr[i] <- strsplit(chrpos[i],":")[[1]][1]
    start[i] <- strsplit(strsplit(chrpos[i],":")[[1]][2],"-")[[1]][1]
    end[i] <- strsplit(strsplit(chrpos[i],":")[[1]][2],"-")[[1]][2]

    window1ptc_pos <- which(substring(mendelian_intron_sequences[i,1],seq(1,str_length(mendelian_intron_sequences[i,1])-3,3),seq(3,str_length(mendelian_intron_sequences[i,1]),3)) %in% c("tag","taa","tga"))*3
    window1ptc[i] <- min(window1ptc_pos)
    window1_no_ptcs[i] <- length(window1ptc_pos)
    window1_all_ptcs[i] <- paste0(window1ptc_pos,collapse=",")

    window2ptc_pos <- which(substring(mendelian_intron_sequences[i,1],seq(2,str_length(mendelian_intron_sequences[i,1])-3,3),seq(4,str_length(mendelian_intron_sequences[i,1]),3)) %in% c("tag","taa","tga"))*3
    window2ptc[i] <- min(window2ptc_pos)
    window2_no_ptcs[i] <- length(window2ptc_pos)
    window2_all_ptcs[i] <- paste0(window2ptc_pos,collapse=",")

    window3ptc_pos <- which(substring(mendelian_intron_sequences[i,1],seq(3,str_length(mendelian_intron_sequences[i,1])-3,3),seq(5,str_length(mendelian_intron_sequences[i,1]),3)) %in% c("tag","taa","tga"))*3
    window3ptc[i] <- min(window3ptc_pos)
    window3_no_ptcs[i] <- length(window3ptc_pos)
    window3_all_ptcs[i] <- paste0(window3ptc_pos,collapse=",")
}

mendelian_intron_sequences_results <- as.data.table(list(mendelian_intron_sequences[,1],intron_length,frame,window1ptc,window2ptc,window3ptc,window1_no_ptcs,window2_no_ptcs,window3_no_ptcs,window1_all_ptcs,window2_all_ptcs,window3_all_ptcs))



names(mendelian_intron_sequences_results) <- c("intron","length","frame","window1ptc","window2ptc","window3ptc","window1_no_ptcs","window2_no_ptcs","window3_no_ptcs","window1_all_ptcs","window2_all_ptcs","window3_all_ptcs")
mendelian_intron_sequences_results$frame <- mendelian_intron_sequences_results$length %% 3 == 0

View(mendelian_intron_sequences_results)

window1_mendelian_intron_sequences_results <- mendelian_intron_sequences_results[,c(1,2,3,4,5,6,9,12)]

fwrite(mendelian_intron_sequences_results,"refseq_genes_introns_ptcs_all_windows.tsv",sep = "\t")
fwrite(window1_mendelian_intron_sequences_results,"refseq_genes_introns_ptcs_window_1.tsv",sep = "\t")


`%nin%` <- Negate(`%in%`)
