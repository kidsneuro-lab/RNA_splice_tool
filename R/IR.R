library(Rsamtools)
    # read in the bam file
    bam <- scanBam(sample_list$bamfile[1])[[1]] # the result comes in nested lists
    # filter reads without match position
    ind <- ! is.na(bam$pos)
    ## remove non-matches, they are not relevant to us
    bam <- lapply(bam, function(x) x[ind])


    introns_of_interest <- GRanges(seqnames = Sample_Refseq_Genes$chrom,
                                   IRanges(start = Sample_Refseq_Genes$region_start,
                                           end = Sample_Refseq_Genes$region_end),
                                   strand = as.factor(Sample_Refseq_Genes$strand))

    introns_of_interest <- IRanges(start = Sample_Refseq_Genes$region_start[1],
                                           end = Sample_Refseq_Genes$region_end[1])

    seqlevelsStyle(introns_of_interest) <- 'UCSC'

    ranges
    sj[[1]]
    introns_of_interest

    coverage(subsetByOverlaps(ranges,introns_of_interest,type = "within"))

    for(i in seq(1,nrow(Sample_Refseq_Genes))){

        intron <- IRanges(start = Sample_Refseq_Genes$region_start[i]-Sample_Refseq_Genes$region_start[i],
                            end = Sample_Refseq_Genes$region_end[i]-Sample_Refseq_Genes$region_start[i])
        ranges <- IRanges(start=bam$pos-Sample_Refseq_Genes$region_start[i],
                          width=bam$qwidth)
        overlap <- subsetByOverlaps(ranges,intron,type = "within")
        cat(Sample_Refseq_Genes$gene_name[i],
              Sample_Refseq_Genes$region_type[i],
              Sample_Refseq_Genes$region_no[i],
              Sample_Refseq_Genes$tx_id[i],
              ":")
        coverage(overlap)
    }

    Sample_Refseq_Genes_Introns <- Refseq_Genes[gene_name %in% Genes$`Gene name` &
                                            region_type == 'intron']
    Sample_Refseq_Genes_Exons <- Refseq_Genes[gene_name %in% Genes$`Gene name` &
                                                    region_type == 'exon']

    MTHFR_exons <- Sample_Refseq_Genes_Exons[gene_name == "MTHFR"]
    MTHFR_introns <- Sample_Refseq_Genes_Introns[gene_name == "MTHFR"]

    #library(ggplot2)

    #MTHFR_exons %>% ggplot() + geom_point(aes(region_start,region_end))

    tx_bump=0

    plot(NULL, xlim=c(min(MTHFR_exons$region_start), max(MTHFR_exons$region_end)),ylim=c(0,1))
    for(i in seq(1,nrow(MTHFR_exons))){
        if(i>1){
          if(MTHFR_exons$tx_id[i-1] != MTHFR_exons$tx_id[i]){
            #tx_bump = tx_bump+0.1
          }
        }
        if(sum(sapply(MTHFR_exons$region_start[i],between,MTHFR_introns$region_start,MTHFR_introns$region_end)) == 1 | sum(sapply(MTHFR_exons$region_end[i],between,MTHFR_introns$region_start,MTHFR_introns$region_end)) == 1){
            #segments(MTHFR_exons$region_start[i],0.20+tx_bump,MTHFR_exons$region_end[i],0.20+tx_bump,lwd=4, col = "red")
        }else{
            segments(MTHFR_exons$region_start[i],0.20+tx_bump,MTHFR_exons$region_end[i],0.20+tx_bump,lwd=4)
        }
            segments(MTHFR_introns$region_start[i],0.20+tx_bump,MTHFR_introns$region_end[i],0.20+tx_bump,lwd=1)
            tx_bump = tx_bump+0.01
    }


    sum(sapply(MTHFR_exons$region_start,between,MTHFR_introns$region_start,MTHFR_introns$region_end)) | sum(sapply(MTHFR_exons$region_end,between,MTHFR_introns$region_start,MTHFR_introns$region_end))


