# ########################################################################### #

#░▄▀▀░▄▀▄▒█▀▄░▀█▀▒▄▀▄▒█▀▄ v0.3.0                                      Mar 2022
#░▀▄▄░▀▄▀░█▀▄░▒█▒░█▀█░█▀▄                        e: rmar4592@uni.sydney.edu.au


#-Clinically-focused relative quantification of aberrant pre-mRNA splicing.


# ########################################################################### #

#Parameters
SampleFile <- "../../../Reports/v0.3 Cortar 03.22/CAPN3/subset/CAPN3_Perth_RNAseq_CAPN3.tsv" #"cortar v030 testing samplefile.tsv"
Assembly <- list("hg19","UCSC", FALSE, 0) # hg38/hg19 | UCSC/1000genomes | paired? | stranded (0,1,2)
Export <- "../../../Reports/v0.3 Cortar 03.22/CAPN3/"

#--Load Environment------------------------------------------------------------
sapply(c("data.table",
         "biomaRt",
         "GenomicFeatures",
         "GenomicAlignments",
         "GenomicRanges",
         "BSgenome.Hsapiens.UCSC.hg38",
         "BSgenome.Hsapiens.UCSC.hg19",
         "BSgenome.Hsapiens.1000genomes.hs37d5",
         "Rsubread",
         "openxlsx"),
       require,warn.conflicts=F,
       quietly=T,
       character.only = TRUE)

source("R/auxiliary_functions.R")


#--Load Data - samples, refseq/ensembl, UCSC/1000 genomes, hg38/hg19-----------
sample_list <- fread(SampleFile, sep ="\t", fill=T)

Transcripts <- unique(sample_list$transcript)

if(Assembly[1] == "hg38"){
    Refseq_Genes <- fread("refseq_introns_exons_hg38.tsv", sep = "\t")
    Ensembl_Genes <- fread("hg38_mart_export_allgenes_chr1-Y.txt", sep = "\t")
    Genome_Assembly = BSgenome.Hsapiens.UCSC.hg38
}

if(Assembly[1] == "hg19"){
    Refseq_Genes <- fread("refseq_introns_exons_hg37.tsv")
    Ensembl_Genes <- fread("hg19_mart_export_allgenes_chr1-Y.txt", sep = "\t")
    if(Assembly[2] == "UCSC") Genome_Assembly <- BSgenome.Hsapiens.UCSC.hg19
    if(Assembly[2] == "1000genomes") Genome_Assembly <- BSgenome.Hsapiens.1000genomes.hs37d5
}


if(sum(sample_list$transcript %in% Refseq_Genes$tx_id) == length(sample_list$genes)){
    message("")
}else{stop("Gene(s) not found.")}

if(sum(sapply(sample_list$bamfile,file.exists)) == nrow(sample_list)){
    message("")
}else{stop("File(s) not found.")}



#--Extract Genes---------------------------------------------------------------
Genes <- Ensembl_Genes[`Gene name` %in% sample_list$genes]

genes.GRanges <- GRanges(seqnames = Genes$`Chromosome/scaffold name`,
                         IRanges(start = Genes$`Gene start (bp)`,
                                 end = Genes$`Gene end (bp)`),
                         strand = Genes$Strand)

#if(Assembly[2] == "UCSC")
    seqlevelsStyle(genes.GRanges) <- 'UCSC'

if(nrow(Genes) == length(unique(sample_list$genes))){
    message("")
}else{stop("Gene(s) not found.")}

#--Extract Reads---------------------------------------------------------------
sj = list()

param <- ScanBamParam(which=genes.GRanges,
                      flag=scanBamFlag(isDuplicate = FALSE,
                                       isSecondaryAlignment = F,
                                       isPaired = Assembly[[3]]))

for (i in 1:nrow(sample_list)) {
    sample <- sample_list[i, sampleID]
    message("\t", sample)
    if(Assembly[4] == 0){
        gal <- readGAlignments(file = sample_list[i, bamfile],
                                   param = param)
    }else{
        gal <- readGAlignmentPairs(file = sample_list[i, bamfile],
                               param = param, strandMode = Assembly[[4]])
    }
    sj[[sample]] <- summarizeJunctions(gal,genome = Genome_Assembly)
    strand(sj[[sample]]) <- mcols(sj[[sample]])[,"intron_strand"]
}

#--Count Reads-----------------------------------------------------------------
combined_sj <- unique(unlist(GRangesList(unlist(sj))))

mcols(combined_sj)[c('score','plus_score','minus_score','intron_motif',
                     'intron_strand')] <- NULL

for (sample in sample_list$sampleID) {
    message("\t",sample)
    mcols(combined_sj)[paste0("sj_", sample)] <- 0

    qryhits <- findOverlaps(sj[[sample]], combined_sj, type = "equal")
    mcols(combined_sj[subjectHits(qryhits)])[paste0("sj_", sample)] <-
        mcols(sj[[sample]][queryHits(qryhits)])[,'score']

}

#--Annotate exons--------------------------------------------------------------
Sample_Refseq_Genes <- Refseq_Genes[tx_id %in% Transcripts &
                                    gene_name %in% Genes$`Gene name` &
                                    region_type == 'intron']

introns_of_interest <- GRanges(seqnames = Sample_Refseq_Genes$chrom,
                               IRanges(start = Sample_Refseq_Genes$region_start,
                                       end = Sample_Refseq_Genes$region_end),
                               strand = as.factor(Sample_Refseq_Genes$strand))

if (Assembly[2] == "UCSC") seqlevelsStyle(introns_of_interest) <- 'UCSC'

    mcols(introns_of_interest)['exon_no'] <- paste(
        Sample_Refseq_Genes$gene_name,
        " In ", Sample_Refseq_Genes$region_no, "", sep="")

    mcols(combined_sj)['exon_range_start'] <- NA
    mcols(combined_sj)['exon_range_end'] <- NA
    mcols(combined_sj)['annotated'] <- 'N'

    qryhits <- findOverlaps(introns_of_interest, combined_sj, type = "start")
    mcols(combined_sj[subjectHits(qryhits)])['exon_range_start'] <- mcols(
        introns_of_interest[queryHits(qryhits)])[,'exon_no']

    qryhits <- findOverlaps(introns_of_interest, combined_sj, type = "end")
    mcols(combined_sj[subjectHits(qryhits)])['exon_range_end'] <- mcols(
        introns_of_interest[queryHits(qryhits)])[,'exon_no']

    qryhits <- findOverlaps(introns_of_interest,combined_sj, type = "equal")
    mcols(combined_sj[subjectHits(qryhits)])['annotated'] <- "Y"

#--Extract Splice Junctions from Samples---------------------------------------

combined.dt <- sortSeqlevels(combined_sj)
combined.dt <- sort(combined_sj)

for (sample_name in names(sj)){
    message("\t",sample_name)

    lhs <- as.data.table(combined.dt[,paste0("sj_", sample_name)])
    rhs <- as.data.table(combined.dt[,paste0("sj_", sample_name)])

    lhs_rhs_merged <- merge(lhs, rhs, by="seqnames", allow.cartesian = T,
                             suffixes = c("","_overlap"))

    lhs_rhs_merged[, overlap_sj := ifelse(
        (start == start_overlap & end == end_overlap &
                strand == strand_overlap) |
            (start == start_overlap & end != end_overlap &
                    strand == strand_overlap) |
            (start != start_overlap & end == end_overlap &
                    strand == strand_overlap), 1, 0)]

    sj_overlap_column_name <- paste0("sj_", sample_name, "_overlap")

    sj_overlap <- lhs_rhs_merged[, .(sj_count = sum(ifelse(overlap_sj == 1,
                                                           get(sj_overlap_column_name), 0))),
                                 by = .(seqnames, start, end, strand)]

    setnames(sj_overlap,'sj_count', paste0("sj_", sample_name, "_overlap"))

    sj_overlap_gr <- GRanges(seqnames = sj_overlap$seqnames,
                                ranges = IRanges(start = sj_overlap$start,
                                                end = sj_overlap$end),
                                strand = sj_overlap$strand)

    mcols(sj_overlap_gr)[sj_overlap_column_name] <- sj_overlap[,..sj_overlap_column_name]

    sj_overlap_gr <- sortSeqlevels(sj_overlap_gr)
    sj_overlap_gr <- sort(sj_overlap_gr)

    mcols(combined.dt)[sj_overlap_column_name] <- mcols(sj_overlap_gr)[,sj_overlap_column_name]
}

#--Calculate splice-junction usage---------------------------------------------

for (sample_name in names(sj)) {
    message("\t",sample_name)

    mcols(combined.dt)[paste0("sj_pct_", sample_name)] <-
            mcols(combined.dt)[,paste0("sj_", sample_name)] /
            mcols(combined.dt)[,paste0("sj_", sample_name, "_overlap")]

    mcols(combined.dt)[paste0("sj_pct_", sample_name)] <-
            nafill(mcols(combined.dt)[,paste0("sj_pct_", sample_name)],
                   type = "const", fill = 0)
}

combineddt <- as.data.table(combined.dt)

combineddt[, event := mapply(gen.exon.range, strand,
                                 exon_range_start, exon_range_end)]
combineddt[, introns := mapply(gen.introns,
                                   exon_range_start, exon_range_end)]
combineddt[, normal := mapply(gen.normal,
                                  exon_range_start, exon_range_end)]
combineddt[, genes := mapply(gen.fetch,
                                 exon_range_start, exon_range_end)]
combineddt[, `:=`(exon_range_start = NULL, exon_range_end = NULL)]

#--Extract Intron Retention from Samples---------------------------------------
rightparams <- c(end,-2,-1,"RHS")
leftparams <- c(start,+1,+2,"LHS")
IRparams <- list(leftparams,rightparams)
for(i in seq(1,2)){
        sj_start_pos <<- GRanges(seqnames = seqnames(combined_sj),
                                ranges = IRanges(
                                    start = IRparams[[i]][[1]](combined_sj) + IRparams[[i]][[2]],
                                    end = IRparams[[i]][[1]](combined_sj) + IRparams[[i]][[3]]),
                                strand = strand(combined_sj))

        rsubreadCounts <- featureCounts(files = sample_list$bamfile,
                                        annot.ext=GRanges.to.SAF(sj_start_pos, 5),
                                        minOverlap=5*2,
                                        allowMultiOverlap=TRUE,
                                        checkFragLength=FALSE,
                                        strandSpecific=Assembly[[4]],

                                        # activating long read mode
                                        isLongRead=F,

                                        # multi-mapping reads
                                        countMultiMappingReads=FALSE,

                                        # unstranded case: for counting only
                                        # non-spliced reads we skip this
                                        #information
                                        isPairedEnd=Assembly[[3]],
                                        nonSplitOnly = T,

                                        # sorting only needed for paired-end reads
                                        autosort=T,
                                        nthreads=1,
                                        tmpDir="temp")

        if(i == 1){
            infocols <<- c("seqnames","start","end","width","strand","annotated",
                           "genes","event","introns","normal")
            combinedintron <- cbind(combineddt[,infocols, with = F],
                                    rsubreadCounts$counts)
            setnames(combinedintron, basename(sample_list$bamfile),
                     paste0("left_", sample_list$sampleID))

        }else if(i == 2){
            combinedintron <- cbind(combinedintron, rsubreadCounts$counts)
            setnames(combinedintron, basename(sample_list$bamfile),
                     paste0("right_", sample_list$sampleID))
        }
}

#--Calculate intron retention--------------------------------------------------

    #combined.dt <- combined.backup.dts
    #combinedintron <- combinedIntrons.backup.dt
    #combined.backup.dts <<- combineddt
    #combinedIntrons.backup.dt <<- combinedintron

    skipping_events_exons <- combineddt[which(event != "unannotated junctions")]
    skipping_events_introns <- combinedintron[which(event != "unannotated junctions")]

    skipping_events_exons[, introns :=  sapply(skipping_events_exons$introns,function(x) unlist(strsplit(x,"-"))[2])]
    combineddt <- rbind(combineddt, skipping_events_exons[!is.na(introns) & normal != "Y"])
    combineddt[, introns :=  sapply(combineddt$introns,function(x) unlist(strsplit(x,"-"))[1])]

    skipping_events_introns[, introns :=  sapply(skipping_events_introns$introns,function(x) unlist(strsplit(x,"-"))[2])]
    combinedintron <- rbind(combinedintron, skipping_events_introns[!is.na(introns) & normal != "Y"])
    combinedintron[, introns :=  sapply(combinedintron$introns,function(x) unlist(strsplit(x,"-"))[1])]



    for (sample_name in names(sj)) {

        message("\t",sample_name)

        overlaps_exons <- combineddt[normal == 'Y']
        overlaps_introns <- combinedintron[normal == 'Y']

        ir_column <- paste0("ir_", sample_name)
        sj_column <- paste0("sj_", sample_name)
        sj_pct_column <- paste0("sj_pct_", sample_name)
        sj_overlap_column <- paste0("sj_", sample_name, "_overlap")
        ns_left_column <- paste0("left_", sample_name)
        ns_right_column <- paste0("right_", sample_name)
        total_ir_column <- paste0("total_ir_", sample_name)
        total_ir_column_1 <- paste0("total_ir_1", sample_name)
        total_sj_column <- paste0("total_sj_", sample_name)
        total_sj_column_1 <- paste0("total_sj_1", sample_name)

        listofoverlaps_exons <- list()
        listofoverlaps_right <- list()
        listofoverlaps_left <- list()

        for(i in seq(1,nrow(combineddt))){
            x <- overlaps_exons[which(overlaps_exons$genes == combineddt$genes[i]&
                                          overlaps_exons$introns == combineddt$introns[i]),
                                sj_overlap_column, with = F]
            y <- overlaps_introns[which(overlaps_introns$genes == combinedintron$genes[i]&
                                            overlaps_introns$introns == combinedintron$introns[i]),
                                  c(ns_left_column), with = F]
            z <- overlaps_introns[which(overlaps_introns$genes == combinedintron$genes[i]&
                                            overlaps_introns$introns == combinedintron$introns[i]),
                                  c(ns_right_column), with = F]
            listofoverlaps_exons[i] <- x
            listofoverlaps_right[i] <- y
            listofoverlaps_left[i] <- z
        }

        listofoverlaps_exons <<- unlist(ifelse(sapply(listofoverlaps_exons, length) == 0, 0, listofoverlaps_exons))
        listofoverlaps_right <<- unlist(ifelse(sapply(listofoverlaps_right, length) == 0, 0, listofoverlaps_right))
        listofoverlaps_left <<- unlist(ifelse(sapply(listofoverlaps_left, length) == 0, 0, listofoverlaps_left))

        combineddt[, c(total_sj_column_1) := listofoverlaps_exons]
        #print(unique(combineddt[,total_sj_column_1, with = F]))
        combinedintron[, c(total_ir_column_1) := (as.numeric(listofoverlaps_right) + as.numeric(listofoverlaps_left))/2]
        #print(unique(combinedintron[,total_ir_column_1, with = F]))

        #print(unique(sapply(c(listofoverlaps_exons,listofoverlaps_left, listofoverlaps_right), class)))
        #print("\n")
        overlaps <- (listofoverlaps_exons + (listofoverlaps_right + listofoverlaps_left)/2)

        #print(unique(sapply(c(listofoverlaps_exons,listofoverlaps_left, listofoverlaps_right), class)))

        combineddt[, c(total_sj_column) := overlaps]
        #print(sapply(c(listofoverlaps_exons,listofoverlaps_left, listofoverlaps_right), class))
        #print(unique(combineddt[,total_sj_column, with = F]))
        combinedintron[, c(total_ir_column) := combineddt[[total_sj_column]]]
        #print(unique(combinedintron[,total_ir_column, with = F]))


        combinedintron[, c(ir_column) :=
                               nafill((get(ns_left_column) + get(ns_right_column))/2 /
                                          get(total_ir_column),
                                      "const", 0)]

        combinedintron[, c(ns_left_column) :=
                               nafill((get(ns_left_column) + get(ns_right_column))/2,
                                      "const", 0)]

        combineddt[, c(sj_pct_column) :=
                         nafill(get(sj_column) / get(total_sj_column),
                                "const", 0)]


        #combinedintron[, c(total_ir_column) :=
        #                nafill((get(ns_left_column) +
        #                      combineddt[[sj_overlap_column]]) +
        #                      (get(ns_right_column) +
        #                      combineddt[[sj_overlap_column]])  / 2,
        #                    "const", 0)]

        #combinedintron[, c(total_ir_column) :=
        #                 nafill(get(total_ir_column) +
        #                          combineddt[[total_sj_column]],
        #                        "const", 0)]



        # combineddt[, c(total_sj_column) := combinedintron[[total_ir_column]]]


        #combineddt[, c(total_sj_column) :=
        #             nafill((get(sj_overlap_column) +
        #                        ((combinedintron[[ns_left_column]] +
        #                         combinedintron[[ns_right_column]])/2)),
        #                        "const", 0)]





    }


#--Combining IR and SJ Data----------------------------------------------------

    #Labelling IR and SJ prior to combining
    combinedintron$SJ_IR <- "IR"
    combinedintron$frame_conserved <- framed(combinedintron$SJ_IR,
                                             combinedintron$start,
                                             combinedintron$end)
    combineddt$SJ_IR <- "SJ"
    combineddt$frame_conserved <- framed(combineddt$SJ_IR,
                                         combineddt$start,
                                         combineddt$end)

    #Combining intron/exon data in a dt and combining read location in a 2nd
    combined_intron_final <- combinedintron[,c(paste0("ir_", names(sj)),
                                               paste0("left_", names(sj)),
                                               paste0("total_ir_", names(sj))),
                                            with=F]
    combined_intron_info <- combinedintron[,c(infocols,
                                              "SJ_IR","frame_conserved"),
                                           with = F]
    combined_dt_final <- combineddt[,c(paste0("sj_pct_", names(sj)),
                                       paste0("sj_", names(sj)),
                                       paste0("total_sj_", names(sj))),
                                    with=F]
    combined_dt_info <- combineddt[,c(infocols,
                                      "SJ_IR","frame_conserved"),
                                   with = F]
    combined_dt_intron_data <- as.data.table(mapply(c,combined_dt_final,
                                                    combined_intron_final))
    combined_dt_intron_info <- rbind(combined_dt_info, combined_intron_info)

    #Combining intron/exon data and read locations into a single dataframe
    combined_dt_intron <- cbind(combined_dt_intron_info,
                                combined_dt_intron_data)
    combined_dt_intron$event[which(combined_dt_intron$introns != "" &
                                       combined_dt_intron$SJ_IR == "IR")] <-
    combined_dt_intron$introns[which(combined_dt_intron$introns != "" &
                                             combined_dt_intron$SJ_IR == "IR")]
    combined_dt_intron$assembly <- Assembly[1]

#--Calculate true IR-----------------------------------------------------------

# get all IR counts for introns and intronic cryptics
# annotate data.table with whether each intron has an intronic cryptic
# use annotation to select counts which will be duplicated - for averaging.
# continue workflow with adjusted IR counts.
# what if there are intronic cryptics on both sides of the intron?

#cryptics.df <- combinedIntronsExons.dt[which(SJ_IR == "SJ" & event != "unannotated junctions" &
#                                annotated == "N")]

#for(i in seq(1,nrow(cryptics.df))){

#  if((cryptics.df$end[i] %in% intronsOfInterest.df$end & cryptics.df$start[i] %in% intronsOfInterest.df$start) == FALSE){
#    if(cryptics.df$end[i] %in% intronsOfInterest.df$end){
#      if(cryptics.df$start[i] > intronsOfInterest.df$start[which(intronsOfInterest.df$end == cryptics.df$end[i])]){
#        cat(c(i," start", "\n"))
#      }
#    }else if(cryptics.df$start[i] %in% intronsOfInterest.df$start){
#      if(cryptics.df$end[i] < intronsOfInterest.df$end[which(intronsOfInterest.df$start == cryptics.df$start[i])]){
#        cat(c(i," end", "\n"))
#      }
#    }
#  }



#}


#cryptics.df[1,2] < intronsOfInterest.df$start[which(intronsOfInterest.df$end == cryptics.df[1,3])]



#--Compare splicing between test and controls----------------------------------
    for(sample in seq(1,nrow(sample_list))){
        combined_dt_intron_test <- combined_dt_intron
        if(sample_list$sampletype[sample] == 'test'){
            family <- sample_list$sampleID[which(
                sample_list$family == sample_list$family[sample])]
            message("\t",family)
            familycols <- paste0("sj_pct_", family)
            familyreadcols <- paste0("sj_", family)

            ctrls <- sample_list$sampleID[which(
                sample_list$family != sample_list$family[sample])]# &
                #sample_list$gene != sample_list$gene[sample])]
            print(ctrls)
            ctrlscols <- paste0("sj_pct_", ctrls)
            print(ctrlscols)
            ctrlsreadcols <- paste0("sj_", ctrls)

            proband <- sample_list$sampleID[sample]

            combined_dt_intron_test$proband <- sample_list$sampleID[sample]

            combined_dt_intron_test$controlavg <- rowMeans(
                combined_dt_intron_test[,ctrlscols, with=F])
            combined_dt_intron_test$controlavgreads <- rowMeans(
                combined_dt_intron_test[,ctrlsreadcols, with=F])
            combined_dt_intron_test <- cbind(combined_dt_intron_test,
                                             controlsd = apply(combined_dt_intron_test[,ctrlscols,
                                                                                       with=F], 1, sd))

            combined_dt_intron_test$difference <- combined_dt_intron_test[,
                                                                          paste0("sj_pct_", proband), with=F] -
                combined_dt_intron_test$controlavg
            combined_dt_intron_test$controln <- length(ctrlscols)
            combined_dt_intron_test$two_sd <- abs(
                combined_dt_intron_test$difference) >
                combined_dt_intron_test$controlsd*2
            combined_dt_intron_test$three_sd <- abs(
                combined_dt_intron_test$difference) >
                combined_dt_intron_test$controlsd*3
            combined_dt_intron_test$four_sd <- abs(
                combined_dt_intron_test$difference) >
                combined_dt_intron_test$controlsd*4

            for (i in seq(1,nrow(combined_dt_intron_test))){
                event_unique_count <- 0
                for (member in familycols){
                    if (combined_dt_intron_test$controlavg[i] == 0 &
                        combined_dt_intron_test[,member, with=F][i] != 0){
                        event_unique_count <- event_unique_count + 1
                    }
                }
                if (event_unique_count >= 1){
                    combined_dt_intron_test$unique[i] <- paste(
                        event_unique_count, "/", length(familycols), sep="")
                }else{
                    combined_dt_intron_test$unique[i] <- ""
                }
            }

            #Extracting the columns in the required order
            combined_dt_intron_final <- combined_dt_intron_test[,c("assembly",
                                                                   infocols[-which(infocols %in% c("normal"))],
                                                                   "SJ_IR","frame_conserved","unique",
                                                                   "proband","difference",familycols,"controlavg",
                                                                   "controlsd","controln", familyreadcols,
                                                                   "controlavgreads","two_sd","three_sd",
                                                                   "four_sd"), with=F]

            #combined_dt_intron_final <- combined_dt_intron_test[,c("assembly",
            #                                    infocols,"SJ_IR",
            #                                    "frame_conserved","unique",
            #                                    "proband","difference",
            #                                    familycols,"controlavg",
            #                                    "controlsd","controln",
            #                                    familyreadcols,"controlavgreads",
            #                                    paste0("total_sj_", family),
            #                                    "two_sd","three_sd",
            #                                    "four_sd"), with=F]

            #Sort by the greatest difference in percentage
            combined_dt_intron_final <- combined_dt_intron_final[order(
                abs(combined_dt_intron_final$difference), decreasing = T),]

            #Exclude Novel Junction and any non-annotated IR reads
            combined_dt_intron_final <- combined_dt_intron_final[!which(
                combined_dt_intron_final$SJ_IR == "IR" &
                    combined_dt_intron_final$annotated == "N")]
            combined_dt_intron_final <- combined_dt_intron_final[!which(
                combined_dt_intron_final$event == "unannotated junctions")]

            report <- T

            testgenes <- unique(sample_list$genes[sample])

            if(report == TRUE){
                #for(i in seq(0,length(testgenes))){
                #combined_dt_intron_final <- combined_dt_intron_final[which(
                #combined_dt_intron_final$genes == testgenes[i])]
                #combined_dt_intron_final <- combined_dt_intron_final[which(
                #combined_dt_intron_final$two_sd == TRUE)]
                generate.report(combined_dt_intron_final[which(
                    combined_dt_intron_final$genes == testgenes)],
                    length(familycols), testgenes, Export,
                    sample_list$sampleID[sample])

                combined_dt_intron_test$normal <- as.character(combined_dt_intron_test$normal)
                fwrite(combined_dt_intron_test[which(combined_dt_intron_test$genes == testgenes)],
                       file = paste0(Export,"/",sample_list$sampleID[sample],"_",testgenes,"_combined_full",
                                    ".tsv"),sep="\t")

                #topdf <- combined_dt_intron_test[genes == testgenes & two_sd == TRUE,
                                        #c(event, sj_pct_332_RNASEH2B_P,controlavg,difference,unique)]

                    rmarkdown::render(
                        input = "test.Rmd",
                        output_file = paste(Export,"/",sample_list$sampleID[sample],"_",testgenes,"_combined_full",
                                            ".pdf", sep=""),
                        params = list("table" = paste(Export,"/",sample_list$sampleID[sample],"_",testgenes,"_combined_full",
                                                      ".tsv", sep=""),
                                      "testgenes" = testgenes,
                                      "control_no" = length(ctrls),
                                      "sample" = proband,
                                      "proband" = paste0("sj_pct_", proband),
                                      "transcripts" = Transcripts,
                                      "refseq" = Refseq_Genes,
                                      "strand" = Refseq_Genes[gene_name == testgenes & canonical == 1,strand][1]))

                #rnaseqgrapher(combined_dt_intron_test, testgenes, sample_list$sampleID[sample], "../../Reports/CDK5RAP3/")

                #}
            }else{
                #Exporting the combineddt dataframe to an excel spreadsheet
                openxlsx::write.xlsx(combined_dt_intron_final,
                                     paste(Export,sample_list$sampleID[sample],
                                           "_combined_dt_",sample_list$assembly[1],
                                           ".xlsx", sep=""),
                                     asTable = T,
                                     overwrite = T)
            }
        }
    }
