# *************************************************************************** #
#░▄▀▀░▄▀▄▒█▀▄░▀█▀▒▄▀▄▒█▀▄ v0.5.0                                      Jun 2022
#░▀▄▄░▀▄▀░█▀▄░▒█▒░█▀█░█▀▄                        e: rmar4592@uni.sydney.edu.au
# ___________________________________________________________________________ #


#--Load Environment------------------------------------------------------------

#Load Parameters
SampleFile <- "src/AGRF_all_blood_subset_samplefile.tsv"
Assembly <- list("hg38","UCSC", TRUE, 2) # hg38/hg19 | UCSC/1000genomes | paired? | stranded (0,1,2)
Export <- "output/"


#Load Dependencies
sapply(c("data.table", "GenomicFeatures", "GenomicAlignments", "formattable",
         "GenomicRanges", "BSgenome.Hsapiens.UCSC.hg38", "openxlsx",
         "BSgenome.Hsapiens.UCSC.hg19", "BSgenome.Hsapiens.1000genomes.hs37d5",
         "magrittr","ggplot2"),
       require,warn.conflicts=F,
       quietly=T,
       character.only = TRUE)

source("R/additional_functions.R")

Refseq_Genes <- fread("ref/refseq_introns_exons_hg38.tsv", sep = "\t")
Ensembl_Genes <- fread("ref/hg38_mart_export_allgenes_chr1-Y.txt", sep = "\t")
Genome_Assembly = BSgenome.Hsapiens.UCSC.hg38
Sample_File <- fread(SampleFile, sep ="\t", fill=T)




#-Select Genes and Introns of Interest-----------------------------------------

#Select genes of interest and create a GRanges object
Genes <- Ensembl_Genes[`Gene name` %in% Sample_File$genes]
genes <- unique(Refseq_Genes[tx_id %in% Sample_File$transcript &
                                 region_type == c('intron'),"gene_name"])
genes.GRanges <- GRanges(seqnames = Genes$`Chromosome/scaffold name`,
                         IRanges(start = Genes$`Gene start (bp)`,
                                 end = Genes$`Gene end (bp)`),
                         strand = Genes$Strand)

seqlevelsStyle(genes.GRanges) <- 'UCSC'


#Select introns of interest for chosen transcript and create a GRanges object
Introns <- Refseq_Genes[tx_id %in% Sample_File$transcript &
                            region_type == c('intron')]

introns.GRanges <- GRanges(seqnames = Introns$chrom,
                           IRanges(start = Introns$region_start,
                                   end = Introns$region_end),
                           strand = Introns$strand)

mcols(introns.GRanges)["intron_no"] <- Introns$region_no
mcols(introns.GRanges)["gene"] <- Introns$gene_name

seqlevelsStyle(introns.GRanges) <- 'UCSC'


#Select introns of interest for other transcripts and create a GRanges object
Introns_Other_Tx <- Refseq_Genes[gene_name %in% unlist(genes) &
                            tx_id %nin% Sample_File$transcript &
                            region_type == c('intron')]

introns_other_tx.GRanges <- GRanges(seqnames = Introns_Other_Tx$chrom,
                           IRanges(start = Introns_Other_Tx$region_start,
                                   end = Introns_Other_Tx$region_end),
                           strand = Introns_Other_Tx$strand)

mcols(introns_other_tx.GRanges)["tx_id"] <- Introns_Other_Tx$tx_version_id
mcols(introns_other_tx.GRanges)["intron_no"] <- Introns_Other_Tx$region_no
mcols(introns_other_tx.GRanges)["gene"] <- Introns_Other_Tx$gene_name

seqlevelsStyle(introns_other_tx.GRanges) <- 'UCSC'

#Create a GRanges object for intron-exon junctions
intron_starts.GRanges <- GRanges(seqnames = Introns$chrom,
                                 IRanges(start = Introns$region_start-4,
                                         end = Introns$region_start+3),
                                 strand = Introns$strand)

intron_ends.GRanges <- GRanges(seqnames = Introns$chrom,
                               IRanges(start = Introns$region_end-3,
                                       end = Introns$region_end+4),
                               strand = Introns$strand)

seqlevelsStyle(intron_starts.GRanges) <- 'UCSC'
seqlevelsStyle(intron_ends.GRanges) <- 'UCSC'


#--Extract all junctions and intron retention (non-split >4bp intron/>4bp exon) reads----
sj = list()
ir = list()

param <- ScanBamParam(
            which=genes.GRanges,flag=scanBamFlag(
                                        isDuplicate = FALSE,
                                        isSecondaryAlignment = FALSE,
                                        isPaired = TRUE))

for (sample_number in 1:nrow(Sample_File)) {
    sample_name <- Sample_File[sample_number, sampleID]
    message("\t", sample_name)

    if(Assembly[4] == 0){
        alignment <- readGAlignments(
                            file = Sample_File[sample_number, bamfile],
                            param = param)
    }else{
        alignment <- readGAlignmentPairs(
                            file = Sample_File[sample_number, bamfile],
                            param = param,
                            strandMode = Assembly[[4]])
    }

    overlaps_intron_starts <- countOverlaps(
                                    intron_starts.GRanges,
                                    alignment,
                                    minoverlap = 8,
                                    type = "any")

    overlaps_intron_ends <- countOverlaps(
                                    intron_ends.GRanges,
                                    alignment,
                                    minoverlap = 8,
                                    type = "any")

    mcols(introns.GRanges)["ir_s"] <- overlaps_intron_starts
    mcols(introns.GRanges)["ir_e"] <- overlaps_intron_ends
    mcols(introns.GRanges)["ir"] <- (overlaps_intron_starts + overlaps_intron_ends)/2
    ir[[sample_name]] <- introns.GRanges
    sj[[sample_name]] <- summarizeJunctions(alignment,genome = Genome_Assembly)
    strand(sj[[sample_name]]) <- mcols(sj[[sample_name]])[,"intron_strand"]
}


#--Aggregate and count junctional reads for each sample------------------------
combined_sj <- unique(unlist(GRangesList(unlist(sj))))
combined_ir <- unique(unlist(GRangesList(unlist(ir))))

mcols(combined_sj)[c('score','plus_score','minus_score','intron_motif',
                     'intron_strand')] <- NULL
mcols(combined_ir)[c('ir','intron_no','genes')] <- NULL

for (sample_name in Sample_File$sampleID) {
    message("\t",sample_name)

    mcols(combined_sj)$SJ_IR <- "SJ"
    mcols(combined_sj)[paste0("count_", sample_name)] <- 0
    qryhits <- findOverlaps(sj[[sample_name]], combined_sj, type = "equal")
    mcols(combined_sj[subjectHits(qryhits)])[paste0("count_", sample_name)] <-
        mcols(sj[[sample_name]][queryHits(qryhits)])[,'score']

    mcols(combined_ir)$SJ_IR <- "IR"
    mcols(combined_ir)[paste0("count_", sample_name)] <- 0
    qryhits <- findOverlaps(ir[[sample_name]], combined_ir, type = "equal")
    mcols(combined_ir[subjectHits(qryhits)])[paste0("count_", sample_name)] <-
        mcols(ir[[sample_name]][queryHits(qryhits)])[,'ir']
}
combined_sj <- c(combined_sj,combined_ir)



#--Extract, annotate, and quantify all events at canonical junctions-----------
combined_sj_sorted <- sortSeqlevels(combined_sj)
combined_sj_sorted <- sort(combined_sj)

events_by_intron <- list()

#Extract and Annotate All Events at Canonical Junctions
for (query_intron in seq(1,nrow(Introns))){
    intron_name <- paste0(Introns$gene_name[query_intron]," intron ",Introns$region_no[query_intron])
    message(paste0("\t",intron_name))

    #For all samples extract all the reads which overlap the query intron jxns
    qryhits_start <- findOverlaps(introns.GRanges[query_intron], combined_sj_sorted, type = "start")
    qryhits_end <- findOverlaps(introns.GRanges[query_intron], combined_sj_sorted, type = "end")

    #Combine the overlapping reads for the query intron
    query_intron.GRanges <- c(combined_sj_sorted[subjectHits(qryhits_start)],
                             combined_sj_sorted[subjectHits(qryhits_end)])

    #Identify the non-canonical splicing junctions
    mcols(query_intron.GRanges)$"intron_jxn_start" <- NA
    mcols(query_intron.GRanges)$"intron_jxn_end" <- NA

    qryhits_start_allintrons <- findOverlaps(introns.GRanges, query_intron.GRanges, type = "start")
    qryhits_end_allintrons <- findOverlaps(introns.GRanges, query_intron.GRanges, type = "end")

    mcols(query_intron.GRanges[subjectHits(qryhits_start_allintrons)])$"intron_jxn_start" <-
        mcols(introns.GRanges[queryHits(qryhits_start_allintrons)])[,'intron_no']
    mcols(query_intron.GRanges[subjectHits(qryhits_end_allintrons)])$"intron_jxn_end" <-
        mcols(introns.GRanges[queryHits(qryhits_end_allintrons)])[,'intron_no']

    #Annotate cryptics that are actually other isoforms
    mcols(query_intron.GRanges)$"annotated" <- ""

    qryhits_equal_othertx <- findOverlaps(introns_other_tx.GRanges, query_intron.GRanges, type = "equal")

    if(length(qryhits_equal_othertx) > 0){
        mcols(query_intron.GRanges[subjectHits(qryhits_equal_othertx)])$"annotated" <- "alternative"
    }

    #Annotate canonical splice junctions
    qryhits_equal_allintrons <- findOverlaps(introns.GRanges, query_intron.GRanges, type = "equal")

    mcols(query_intron.GRanges[subjectHits(qryhits_equal_allintrons)])$"annotated" <- "canonical"

    #Convert intron GRange into data.table
    query_intron.dt <- as.data.table(query_intron.GRanges)
    query_intron.dt <- unique(query_intron.dt)

    #Calculate the proportion of splicing at the intron each event represents
    for (sample_name in Sample_File$sampleID) {
        query_intron.dt[,paste0("pct_", sample_name)] <-
            nafill((query_intron.dt[[paste0("count_", sample_name)]]/
                    sum(query_intron.dt[[paste0("count_", sample_name)]])),
                    type="const",
                    fill=0,
                    nan=NA)
    }

    #Annotate with assembly, gene, intron_no
    query_intron.dt$assembly <- Assembly[[1]]
    query_intron.dt$gene <- Introns$gene_name[query_intron]
    query_intron.dt$intron_no <- Introns$region_no[query_intron]

    #Annotate with event
    query_intron.dt$event <- eventAnnotation(query_intron.dt)

    #Annotate with frame
    query_intron.dt$frame_conserved <- framed(query_intron.dt)

    #Append intron analysis to all events by intron
    events_by_intron[[intron_name]] <- query_intron.dt

}


#Merge all introns into a single data.table
all_splicing_events <- rbindlist(events_by_intron)


#--Compare splicing between test and controls and Generate Report--------------
for(sample_number in seq(1,nrow(Sample_File))){

#Initialise new copy of the all_splicing_events dataset
    all_splicing_events_sample <- all_splicing_events
    if(Sample_File$sampletype[sample_number] == 'test'){

    #Identify columns for the proband and family members
        proband <- Sample_File$sampleID[sample_number]
        family <- Sample_File$sampleID[which(
            Sample_File$family == Sample_File$family[sample_number])]
        message("\t",family)
        familycols <- paste0("pct_", family)
        familyreadcols <- paste0("count_", family)

    #Identify columns for the controls
        ctrls <- Sample_File$sampleID[which(
            Sample_File$family != Sample_File$family[sample_number])]
        ctrlscols <- paste0("pct_", ctrls)
        ctrlsreadcols <- paste0("count_", ctrls)

#Initialise various columns
    #Proband
        all_splicing_events_sample$proband <- Sample_File$sampleID[sample_number]

    #Control average pct, read count, sd, and n
        all_splicing_events_sample$controlavg <- rowMeans(
            all_splicing_events_sample[,..ctrlscols])
        all_splicing_events_sample$controlavgreads <- rowMeans(
            all_splicing_events_sample[,..ctrlsreadcols])
        all_splicing_events_sample <- cbind(all_splicing_events_sample,
                                         controlsd = apply(all_splicing_events_sample[,..ctrlscols], 1, sd))
        all_splicing_events_sample$controln <- length(ctrlscols)

    #Difference between proband and control average
        all_splicing_events_sample$difference <- all_splicing_events_sample[,paste0("pct_", proband),with=F] -
            all_splicing_events_sample$controlavg

    #Compute the standard deviation thresholds
        all_splicing_events_sample$two_sd <- abs(
            all_splicing_events_sample$difference) > all_splicing_events_sample$controlsd*2
        all_splicing_events_sample$three_sd <- abs(
            all_splicing_events_sample$difference) > all_splicing_events_sample$controlsd*3
        all_splicing_events_sample$four_sd <- abs(
            all_splicing_events_sample$difference) >
            all_splicing_events_sample$controlsd*4


    #Identify unique events
        for (i in seq(1,nrow(all_splicing_events_sample))){
            event_unique_count <- 0
            for (member in familycols){
                if (all_splicing_events_sample$controlavg[i] == 0 &
                    all_splicing_events_sample[,..member][i] != 0){
                    event_unique_count <- event_unique_count + 1
                }
            }
            if (event_unique_count >= 1){
                all_splicing_events_sample$unique[i] <- paste(
                    event_unique_count, "/", length(familycols), sep="")
            }else{
                all_splicing_events_sample$unique[i] <- ""
            }
        }

    #Order columns and sort by the greatest difference
        all_splicing_events_sample <- all_splicing_events_sample[order(
            abs(all_splicing_events_sample$difference), decreasing = T),
            c("assembly",
              "proband",
              "seqnames",
              "start",
              "end",
              "width",
              "strand",
              "gene",
              "event",
              "annotated",
              "frame_conserved",
              "unique",
              "difference",
              familycols,
              "controlavg",
              "controlsd",
              "controln",
              familyreadcols,
              "controlavgreads",
              "two_sd",
              "three_sd",
              "four_sd",
              "intron_no",
              "SJ_IR"),with=F]

    #Set report outputs and parameters
        report <- T
        splicing_diagnostics_report <- F
        full_all_genes_report <- F
        testgenes <- unique(Sample_File$genes[sample_number])
        normalSpliceMap(all_splicing_events_sample, familycols[1], testgenes)

    #Generate filtered excel spreadsheet +/- summary html
        #Excel spreadsheet
        if(report == TRUE){
            generate.report(all_splicing_events_sample[gene %in% testgenes],
                            length(familycols),
                            gene = testgenes,
                            Export,
                            Sample_File$sampleID[sample_number])

            fwrite(all_splicing_events_sample,#[which(combined_dt_intron_test$genes == testgenes)],
                   file = paste0(Export,"/",Sample_File$sampleID[sample_number],"_",testgenes,"_combined_full",
                                 ".tsv"),sep="\t")

        #--HTML summary report---- You are up to here! Nearly done! -----------
            if(splicing_diagnostics_report == TRUE){

                report_table <- all_splicing_events_sample

                names(report_table)[which(names(report_table) == paste0("pct_",proband))] <- "probandpct"

                report_table$seqnames <- paste0(report_table$seqnames,":",report_table$start,"-",report_table$end)
                report_table$difference <- formattable::percent(report_table$difference, digits=1)
                report_table$probandpct <- formattable::percent(report_table$probandpct, digits=1)
                report_table$controlavg <- formattable::percent(report_table$controlavg, digits=1)


                sig_introns <- unique(report_table[two_sd == TRUE & abs(difference) > 0.05, intron_no])

                no_introns <- length(sig_introns)
                strand <- Refseq_Genes[tx_id %in% unique(Sample_File$transcript), strand][1]
                full_gene_figure(testgenes, unique(Sample_File$transcript), sig_introns, strand)
                # Close the graphics device
                dev.off()



                normaltable <- report_table[two_sd == TRUE & intron_no %in% sig_introns & SJ_IR == "SJ" & annotated == "canonical"]

                normalevents <- normaltable[,.(seqnames,event,difference,probandpct,controlavg,frame_conserved)][order(abs(difference), decreasing = T)]

                names(normalevents) <- c("splice-junction","event","difference","proband","controls","frame")

                normalChangeBarPlot(normaltable)

                splicingFrameConsequences(normaltable, sig_introns[[1]])

                sig_introns_list <- list()

                for(intron in sig_introns){

                    message("## Intron ",intron)

                    introntable <- report_table[two_sd == TRUE & introns == intron]
                    insigrow <- colSums(table[two_sd == FALSE & introns == intron,.(difference,probandpct,controlavg)])
                    insigrow <- as.data.table(t(c("","NS events < 2sd",insigrow,"")))
                    insigrow$difference <- formattable::percent(insigrow$difference, digits=1)
                    insigrow$probandpct <- formattable::percent(insigrow$probandpct, digits=1)
                    insigrow$controlavg <- formattable::percent(insigrow$controlavg, digits=1)


                    events <- introntable[,.(seqnames,event,difference,probandpct,controlavg,frame_conserved)][order(abs(difference), decreasing = T)]

                    names(events) <- c("splice-junction","event","difference","proband","controls","frame")

                    sig_introns_list[[intron]] <- rbind(events,insigrow, use.names=FALSE)


                sds <- table[,lapply(.(two_sd,three_sd,four_sd),sum)]


                rmarkdown::render(
                     input = "R/splicing_diagnostics_report.Rmd",
                     output_file = paste(Export,"/",Sample_File$sampleID[sample_number],"_",testgenes,"_combined_full",
                                         ".html", sep=""),
                     params = list("normal" = normalevents,
                                   "aberrant" = sig_introns_list,
                                   "testgenes" = testgenes,
                                   "control_no" = length(ctrls),
                                   "sample" = Sample_File$sampleID[[sample_number]],
                                   "transcript" = Sample_File$transcript[[sample_number]],
                                   "refseq" = Refseq_Genes,
                                   "strand" = strand,
                                   "sds" = sds,
                                   "details" = "../../../Reports/v0.4 Cortar 05.22/AGRF Batch 3/Blood/AGRF_blood_batch_3_details.txt"
                                   )
                   )
            }

    #Full dataset export - no report option
        }else if(full_all_genes_report == T){
            #Exporting the combineddt dataframe to an excel spreadsheet
            openxlsx::write.xlsx(all_splicing_events_sample,
                                 paste(Export,Sample_File$sampleID[sample_number],
                                       "_combined_dt_",Sample_File$assembly[1],
                                       ".xlsx", sep=""),
                                 asTable = T,
                                 overwrite = T)
        }
    }
    }
}
