############################################################################# #

#░▄▀▀░▄▀▄▒█▀▄░▀█▀▒▄▀▄▒█▀▄ v0.1.0                                      Nov 2021
#░▀▄▄░▀▄▀░█▀▄░▒█▒░█▀█░█▀▄                        e: rmar4592@uni.sydney.edu.au


## Clinically-focused relative quantification of aberrant pre-mRNA splicing.

## TO DO:

  #== Errors ==================================================================
        # Fetching ensembl object does not load/wait for the code
            # Use a biomart download?

  #== Features ================================================================
        # Don’t compare SJs of two samples testing the same gene
        # Call proximal variants from RNAseq to measure allele bias
        # Compare VCFs with RNAseq BAMs for AGRF cases
        # How can we distinguish between IR and intronic cryptics
            # Split reads! - don’t we already do this?

  #== Optimisation ============================================================
        # Optimise splice junction analysis for multiple samples with same gene
        # Optimise a gene panel(s) to be run on the tool

  #== Reporting ===============================================================
        # Create code to automate creation of the summary excel sheet
        # Enable reports for multiple genes to be analysed for a single sample
        # Combine all of the excel sheets into a single file with summary
        # Clean up code – especially spreadsheet processing

  #== Other ===================================================================
        # Re-run all old cases including Broad on newer versions of the tool
        # Slideshow for podcast
        # Time running the code for short and long genes
        # Release as an R package on GitHub (paper?)


###############################################################################


#--Run Cortar------------------------------------------------------------------

# Load samples for analysis (will be depreciated)
sample_file.tsv <-
    "../../Reports/New Feature Development/NFD_Dataset_Broad_26062013.tsv"

# Runs the full pipeline
cortar <- function(sampledata,exportlocation,sampleID="",genes="",alt_genes="",
                   assembly="",variants="",VCF="", controls=F, chx = ""){
    load.environment()
    sampleList.dt <<- load.samples(sample_file.tsv,sampleID,genes,alt_genes,
                                   assembly,variants,VCF)
    sampleList.dt <<- load.controls(sampleList.dt,controls,chx)
    validate.samples(sampleList.dt)
    set.assembly(sampleList.dt)
    sampleList.dt <<- validate.genes(sampleList.dt)
    genes.GRanges <<- extract.gene.ranges(sampleList.dt)
    sjs.GRanges <<- obtain.SJs(sampleList.dt)
    combined.dt <<- combine.SJs(sampleList.dt,sjs.GRanges)
    combined.dt <<- annotate.exons(sampleList.dt,combined.dt)
    combined.dt <<- extract.SJs(combined.dt,sjs.GRanges)
    combined.dt <<- calculate.SJs(combined.dt,sjs.GRanges)
    combinedIntrons.dt <<- extract.IR(sampleList.dt, combined.dt)
    combinedIntrons.dt <<- calculate.IR(combined.dt,combinedIntrons.dt,
                                        sjs.GRanges)
    combinedIntronsExons.dt <<- combine.IR.SJ(combinedIntrons.dt,combined.dt,
                                              sampleList.dt,sjs.GRanges)
    compare.samples(sampleList.dt, combinedIntronsExons.dt, exportlocation)
}


#--Load Environment------------------------------------------------------------

# Loads the required packages and variables required for validating input
load.environment <- function(){
    message("Loading environment...")
    required_packages <- c("data.table","biomaRt","GenomicFeatures",
                           "GenomicAlignments","GenomicRanges",
                           "BSgenome.Hsapiens.UCSC.hg38",
                           "BSgenome.Hsapiens.1000genomes.hs37d5",
                           "Rsubread","openxlsx",
                           "BSgenome.Hsapiens.UCSC.hg19")

    sapply(required_packages,require,warn.conflicts=F,quietly=T,
           character.only = TRUE)
    sampleFileColumns <<- c("assembly","RNAseq","tissue","sampleID","family",
                            "sampletype","genes","alt_genes","bamfile",
                            "variants","VCF")
    infocols <<- c("seqnames","start","end","width","strand","annotated",
                   "genes","event","introns")
    message("\t","Environment loaded.","\n")
}

#--Load samples----------------------------------------------------------------

# Loads the .tsv file containing multiple RNAseq files to be analysed
load.samples <- function(samples, sampleID="",genes="",alt_genes="",
                         assembly="",variants="",VCF=""){
    message("Loading samples...")
    if(grepl("\\.tsv$", samples)){
        sampleList.dt <- fread(samples, sep ="\t", fill=T)
        message("\t","Fetching ", samples)
    }else{
        sampleList.dt <- data.table(assembly,sampleID,1,"test",genes,alt_genes,
                                    samples,variants,VCF)
        names(sampleList.dt) <- sampleFileColumns
    }
    message("\t","\t","Samples loaded.","\n")
    return(sampleList.dt)
}

#--Generate control dataset----------------------------------------------------

# Adds additional controls based on chosen parameters (age, sex, tissue)
load.controls <- function(samples,controls,chx=""){
    message("Loading controls...")

    if(controls == T){
        controlRNAseq <- fread("RNAseq_sample_database.tsv", sep ="\t", fill=T)
        controlRNAseq <- controlRNAseq[which(
                                    tissue %in% unique(samples$tissue) &
                                    gene %nin% unique(samples$genes) &
                                    chx_dmso %in% c("",chx) &
                                    sampleID %nin% unique(samples$sampleID))]
        controlRNAseq <- controlRNAseq[,c("assembly","RNAseq","tissue",
                                          "sampleID","familyID","sampletype",
                                          "gene","alt_genes",
                                          "bamfile_batch_subset","variants",
                                          "VCF_affected_gene"),with=F]
        samplecontrolRNAseq <- rbind(samples,controlRNAseq, use.names=F)
        message("\t","Additional controls included.","\n")
    }else{
        samplecontrolRNAseq <- samples
        message("\t","Additional controls excluded.","\n")
    }
    return(samplecontrolRNAseq)
}

#--Load Assembly Specific Parameters-------------------------------------------

# Pre-loads variables to either hg38 or hg19 depending on desired assembly
set.assembly <- function(samples){

    message("Setting assembly parameters...")

    if(identical(unique(samples$assembly),"hg38")){

        refseq_introns_exons <<- fread("refseq_introns_exons_hg38.tsv.gz",
                                       sep = "\t")
        ensembl_gene_list <<- fread("hg38_mart_export_allgenes_chr1-Y.txt",
                                    sep = "\t")
        #ensembl <<- useEnsembl(biomart = "genes", host = "www.ensembl.org",
        #                       version = "104",
        #                       dataset = "hsapiens_gene_ensembl")
        message("\t","Parameters set for hg38 assembly.","\n")

    }else if (identical(unique(samples$assembly),"hg19")){

        refseq_introns_exons <<- fread("refseq_introns_exons_hg37.tsv.gz")
        ensembl_gene_list <<- fread("hg19_mart_export_allgenes_chr1-Y.txt",
                                    sep = "\t")
        #ensembl <<- useEnsembl(biomart = "genes", host = "grch37.ensembl.org",
        #                       dataset = "hsapiens_gene_ensembl")
        message("\t","Parameters set for hg19 assembly.","\n")

    }else{

        stop("Single assembly not specified. Provide either hg19 or hg38.")
    }
}


#--Validate Files and Genes----------------------------------------------------

#Checks the format of the multiple sample input file is correct
validate.samples <- function(samples){
    message("Validating sample file...")
    if(identical(names(samples),sampleFileColumns) == T){
        message("\t","Sample file format valid.","\n")
        rm(sampleFileColumns, pos = 1)
    }else{
        stop('Invalid file format. Ensure column headers are correct and file
             is tab-delimited')
    }
}

#Checks that genes are in RefSeq and Ensembl and finds matches in alt_genes
validate.genes <- function(samples){
    message("Checking gene names...")

    sampleList.dt$refseq_genes <- sampleList.dt$genes
    sampleList.dt$ensembl_genes <- sampleList.dt$genes
    if (length(sampleList.dt$genes[which(sampleList.dt$genes %nin%
                                         refseq_introns_exons$gene_name)])==0 &
        length(sampleList.dt$genes[which(sampleList.dt$genes %nin%
                                         ensembl_gene_list$`Gene name`)])==0){
        message("\t",
                "All gene names recognised by RefSeq & Ensembl databases.",
                "\n")
        return(sampleList.dt)
    }else{
        message("\t","Genes not in RefSeq and/or Ensembl: ",length(
            sampleList.dt$genes[which(sampleList.dt$genes %nin%
                                          refseq_introns_exons$gene_name)]))
    }

    message("\t","Querying alt_genes against RefSeq database...")

    for(i in sampleList.dt$genes[which(sampleList.dt$genes %nin%
                                       refseq_introns_exons$gene_name)]){
        alt_genes <- unlist(strsplit(sampleList.dt$alt_genes[which(
            sampleList.dt$genes == i)],split=","))
        if(length(alt_genes) == 0){
            stop("Gene ","\"",i,"\""," not in RefSeq database.",
                 "Try providing \"alt_genes\" or changing the gene name.")
            break
        }
        for(j in seq(1,length(alt_genes))){
            if(alt_genes[j] %in% refseq_introns_exons$gene_name){
                sampleList.dt$refseq_genes[which(
                                    sampleList.dt$genes== i)] <- alt_genes[j]
                break
            }else if(j == length(alt_genes)){
                stop("Gene ","\"",i,"\"",
                     " and alt_genes not in RefSeq database.")
            }
        }
    }
    if(identical(sampleList.dt$genes,sampleList.dt$refseq_genes)){
        message("\t","\t","All gene names recognised by RefSeq database.","\n")

    }else{
        message("\t","\t","RefSeq gene names found: ",sum(sampleList.dt$genes != sampleList.dt$refseq_genes),"\n")
    }

    message("\t","Querying alt_genes against Ensembl database...")

    for(i in sampleList.dt$genes[which(sampleList.dt$genes %nin% ensembl_gene_list$`Gene name`)]){
        alt_genes <- unlist(strsplit(sampleList.dt$alt_genes[which(sampleList.dt$genes == i)],split=","))

        if(length(alt_genes) == 0){
            stop("Gene ","\"",i,"\""," not in Ensembl database. Try providing \"alt_genes\" or changing the gene name.")
            break
        }
        for(j in seq(1,length(alt_genes))){

            if(alt_genes[j] %in% ensembl_gene_list$`Gene name`){
                sampleList.dt$ensembl_genes[which(sampleList.dt$genes== i)] <- alt_genes[j]
                break

            }else if(j == length(alt_genes)){
                stop("Gene ","\"",i,"\""," and alt_genes not in RefSeq database.")
            }
        }
    }

    if(identical(sampleList.dt$genes,sampleList.dt$ensembl_genes)){
        message("\t","\t","All gene names recognised by Ensembl database.","\n")

    }else{
        message("\t","\t","Ensembl gene names found: ",sum(sampleList.dt$genes != sampleList.dt$ensembl_genes),"\n")
    }

    sampleList.dt$genes <- sampleList.dt$ensembl_genes
    return(sampleList.dt)
}


#--Extract Gene Names----------------------------------------------------------

#Creates a granges object with the coordinates of the genes of interest
extract.gene.ranges <- function(samples){
    message("Extracting gene ranges...")
    genes <- getBM(attributes = c('external_gene_name','chromosome_name','start_position','end_position','strand'),
             filters = c('chromosome_name','external_gene_name'),
             values = list(c(as.character(seq(1,22)), "X", "Y"), unlist(sapply(unique(samples$genes), strsplit, split = ","))),
             mart = ensembl,
             uniqueRows = T)

    genes.GRanges <- GRanges(seqnames = genes$chromosome_name, IRanges(start = genes$start_position, end = genes$end_position), strand = genes$strand)
    #seqlevelsStyle(genes.GRanges) <- 'UCSC'

    if(unique(samples$assembly) == "hg38"){
        seqlevelsStyle(genes.GRanges) <- 'UCSC'
    }

    if(nrow(genes) != length(unlist(sapply(unique(samples$genes), strsplit, split = ",")))){
        warning("Gene ranges extracted for only ",nrow(genes)," of ",length(unlist(sapply(unique(samples$genes), strsplit, split = ",")))," genes.","\n")

    }else{
        message("\t","All gene ranges extracted successfully.","\n")
    }
    return(genes.GRanges)
}


#--Obtain Splice Junctions-----------------------------------------------------

# Obtain splice junctions for all samples for selected genes
obtain.SJs <- function(samples){
    message("Obtaining splice junctions...")
    sj = list()

    param <- ScanBamParam(which = genes.GRanges,
                      flag=scanBamFlag(isDuplicate = FALSE,
                                       isSecondaryAlignment = F,
                                       isPaired = T))

    for (row_number in 1:nrow(samples)) {
        message("\t",samples[row_number, sampleID])

        sampleID <- samples[row_number, sampleID]
        bamfile <- samples[row_number, bamfile]

        gal <- readGAlignmentPairs(file = bamfile, param = param, strandMode = 2)

        #hg38
        if (unique(samples$assembly) == "hg38"){
            sj[[sampleID]] <- summarizeJunctions(gal, genome = BSgenome.Hsapiens.UCSC.hg38)

        #hg19
        }else if (unique(samples$assembly) == "hg19"){
            sj[[sampleID]] <- summarizeJunctions(gal, genome = BSgenome.Hsapiens.1000genomes.hs37d5)
            #sj[[sampleID]] <- summarizeJunctions(gal, genome = BSgenome.Hsapiens.UCSC.hg19)
        }
        strand(sj[[sampleID]]) <- mcols(sj[[sampleID]])[,"intron_strand"]
    }
    message("\t","\t","Splice junctions obtained.","\n")
    return(sj)
}


#--Combine Splice Junctions----------------------------------------------------

# Obtain all reads for combined.dt splice junctions
combine.SJs <- function(samples,sj){
    message("Combining splice junctions...")

    combined_grlist <- GRangesList(unlist(sj))
    combined.dt <- unique(unlist(combined_grlist))

    mcols(combined.dt)[c('score','plus_score','minus_score','intron_motif','intron_strand')] <- NULL



    for (row_number in 1:nrow(samples)) {
        sampleID <- samples[row_number, sampleID]
        message("\t",sampleID)
        mcols(combined.dt)[paste0("sj_", sampleID)] <- NA

        qryhits <- findOverlaps(sj[[sampleID]], combined.dt, type = "equal")
        mcols(combined.dt[subjectHits(qryhits)])[paste0("sj_", sampleID)] <- mcols(sj[[sampleID]][queryHits(qryhits)])[,'score']

        mcols(combined.dt)[paste0("sj_", sampleID)] <- nafill(mcols(combined.dt)[, paste0("sj_", sampleID)], type = "const", fill = 0)
    }
    message("\t","\t","Splice junctions combined.dt.","\n")
    return(combined.dt)
}


#--Annotate Exons--------------------------------------------------------------
annotate.exons <- function(samples,combinedsamples){
    message("Annotating exons...")
    refseq_introns_exons_of_interest <<- refseq_introns_exons[canonical == 1 & gene_name %in% unique(samples$refseq_genes) & region_type == 'intron']

    introns_of_interest <<- GRanges(
        seqnames = refseq_introns_exons_of_interest$chrom,
        IRanges(start = refseq_introns_exons_of_interest$region_start, end = refseq_introns_exons_of_interest$region_end),
        strand = as.factor(refseq_introns_exons_of_interest$strand)
    )

    if (unique(samples$assembly) == "hg38"){
        seqlevelsStyle(introns_of_interest) <- 'UCSC'
    }

    mcols(introns_of_interest)['exon_no'] <- paste(refseq_introns_exons_of_interest$gene_name, " In ", refseq_introns_exons_of_interest$region_no, "", sep="")

    mcols(combinedsamples)['exon_range_start'] <- NA
    mcols(combinedsamples)['exon_range_end'] <- NA
    mcols(combinedsamples)['annotated'] <- 'N'

    qryhits <<- findOverlaps(introns_of_interest, combinedsamples,  type = "start")
    mcols(combinedsamples[subjectHits(qryhits)])['exon_range_start'] <- mcols(introns_of_interest[queryHits(qryhits)])[,'exon_no']

    qryhits <<- findOverlaps(introns_of_interest, combinedsamples,  type = "end")
    mcols(combinedsamples[subjectHits(qryhits)])['exon_range_end'] <- mcols(introns_of_interest[queryHits(qryhits)])[,'exon_no']

    qryhits <<- findOverlaps(introns_of_interest, combinedsamples,  type = "equal")
    mcols(combinedsamples[subjectHits(qryhits)])['annotated'] <- "Y"
    message("\t","Exons annotated.","\n")
    return(combinedsamples)
}


#--Extract Splice Junctions from Samples---------------------------------------
extract.SJs <- function(combinedsamples,sj){

    message("Extracting annotated splice junctions...")
    combined.dt <- sortSeqlevels(combinedsamples)
    combined.dt <- sort(combinedsamples)

    for (sample_name in names(sj)){
        message("\t",sample_name)

        lhs <<- as.data.table(combined.dt[,paste0("sj_", sample_name)])
        rhs <<- as.data.table(combined.dt[,paste0("sj_", sample_name)])

        lhs_rhs_merged <<- merge(lhs, rhs, by="seqnames", allow.cartesian = T, suffixes = c("","_overlap"))

        lhs_rhs_merged[, overlap_sj := ifelse(
            (start == start_overlap & end == end_overlap & strand == strand_overlap) |
            (start == start_overlap & end != end_overlap & strand == strand_overlap) |
            (start != start_overlap & end == end_overlap & strand == strand_overlap), 1, 0)]

        sj_overlap_column_name <- paste0("sj_", sample_name, "_overlap")

        sj_overlap <- lhs_rhs_merged[, .(sj_count = sum(ifelse(overlap_sj == 1, get(sj_overlap_column_name), 0))), by = .(seqnames, start, end, strand)]
        setnames(sj_overlap, 'sj_count', paste0("sj_", sample_name, "_overlap"))

        sj_overlap_gr <- GRanges(seqnames = sj_overlap$seqnames,
                                ranges = IRanges(start = sj_overlap$start, end = sj_overlap$end),
                                strand = sj_overlap$strand)

        mcols(sj_overlap_gr)[sj_overlap_column_name] <- sj_overlap[,..sj_overlap_column_name]

        sj_overlap_gr <- sortSeqlevels(sj_overlap_gr)
        sj_overlap_gr <- sort(sj_overlap_gr)

        mcols(combined.dt)[sj_overlap_column_name] <- mcols(sj_overlap_gr)[,sj_overlap_column_name]
    }
    message("\t","\t","Annotated splice junctions extracted.","\n")
    return(combined.dt)
}


#--Calculate usage of splice junctions-----------------------------------------
calculate.SJs <- function(combinedsamples, sj){
    message("Calculating splice junction usage...")
    for (sample_name in names(sj)) {
        message("\t",sample_name)
        mcols(combinedsamples)[paste0("sj_pct_", sample_name)] <- mcols(combinedsamples)[,paste0("sj_", sample_name)] / mcols(combinedsamples)[,paste0("sj_", sample_name, "_overlap")]
        mcols(combinedsamples)[paste0("sj_pct_", sample_name)] <- nafill(mcols(combinedsamples)[,paste0("sj_pct_", sample_name)], type = "const", fill = 0)
    }

    combineddt <- as.data.table(combinedsamples)

    combineddt[, event := mapply(gen.exon.range, strand, exon_range_start, exon_range_end)]
    combineddt[, introns := mapply(gen.introns, exon_range_start, exon_range_end)]
    combineddt[, genes := mapply(gen.fetch, exon_range_start, exon_range_end)]
    combineddt[, `:=`(exon_range_start = NULL, exon_range_end = NULL)]

    infocols <- c("seqnames","start","end","width","strand","annotated","genes","event","introns")

    message("\t","\t","Splice junction usage calculated.","\n")
    return(combineddt)
}


#--Extract Intron Retention from Samples---------------------------------------
extract.IR <- function(samples, combineddt){

    rightparams <- c(end,-2,-1,"RHS")
    leftparams <- c(start,+1,+2,"LHS")
    IRparams <- list(leftparams,rightparams)
    message("Extracting intron retention reads...")
    for(i in seq(1,2)){
        sj_start_pos <- GRanges(seqnames = seqnames(combineddt),
                            ranges = IRanges(start = IRparams[[i]][[1]](combineddt) + IRparams[[i]][[2]],
                                             end = IRparams[[i]][[1]](combineddt) + IRparams[[i]][[3]]),
                            strand = strand(combineddt))

        rsubreadCounts <- featureCounts(files = samples$bamfile,
                                    annot.ext=GRanges.to.SAF(sj_start_pos, 5),
                                    minOverlap=5*2,
                                    allowMultiOverlap=TRUE,
                                    checkFragLength=FALSE,
                                    strandSpecific=2,

                                    # activating long read mode
                                    isLongRead=F,

                                    # multi-mapping reads
                                    countMultiMappingReads=TRUE,

                                    # unstranded case: for counting only non spliced reads we
                                    # skip this information
                                    isPairedEnd=T,
                                    nonSplitOnly = T,

                                    # sorting only needed for paired-end reads
                                    autosort=T,
                                    nthreads=1,
                                    tmpDir="temp")

        if(i == 1){
            combinedintron <- cbind(combineddt[,infocols, with = F], rsubreadCounts$counts)
            setnames(combinedintron, basename(samples$bamfile), paste0("left_", samples$sampleID))

        }else if(i == 2){
            combinedintron <- cbind(combinedintron, rsubreadCounts$counts)
            setnames(combinedintron, basename(samples$bamfile), paste0("right_", samples$sampleID))
        }
        message("\t","Extracted ",IRparams[[i]][[4]]," IR reads.")
    }
    message("\t","\t","Intron Retention reads extracted.","\n")
    return(combinedintron)
}


#--Calculate intron retention--------------------------------------------------
calculate.IR <- function(combined.dt,combinedintron,sj){
    message("Calculating intron retention...")
    for (sample_name in names(sj)) {
        ir_column <- paste0("ir_", sample_name)
        sj_column <- paste0("sj_", sample_name)
        ns_left_column <- paste0("left_", sample_name)
        ns_right_column <- paste0("right_", sample_name)

        combinedintron[, c(ir_column) := nafill(( (get(ns_left_column)/(get(ns_left_column) + combineddt[[sj_column]])) +
                                                   (get(ns_right_column)/(get(ns_right_column) + combineddt[[sj_column]])) ) / 2, "const", 0)]
    }
    message("\t","Intron Retention calculated.","\n")
    return(combinedintron)
}


#--Combining IR and SJ Data----------------------------------------------------
combine.IR.SJ <- function(combinedintron, combineddt, samples, sj){
    message("Combining intron retention and splice junction data...")
    #Labelling IR and SJ prior to combining
    combinedintron$SJ_IR <- "IR"
    combinedintron$frame_conserved <- framed(combinedintron$SJ_IR,combinedintron$start,combinedintron$end)
    combineddt$SJ_IR <- "SJ"
    combineddt$frame_conserved <- framed(combineddt$SJ_IR,combineddt$start,combineddt$end)

    #Combining intron/exon data in one data frame and combining read location in a second
    combined_intron_final <- combinedintron[,c(paste0("ir_", names(sj)),paste0("left_", names(sj))), with=F]
    combined_intron_info <- combinedintron[,c(infocols,"SJ_IR","frame_conserved"), with = F]
    combined_dt_final <- combineddt[,c(paste0("sj_pct_", names(sj)),paste0("sj_", names(sj))), with=F]
    combined_dt_info <- combineddt[,c(infocols,"SJ_IR","frame_conserved"), with = F]
    combined_dt_intron_data <- as.data.table(mapply(c,combined_dt_final,combined_intron_final))
    combined_dt_intron_info <- rbind(combined_dt_info, combined_intron_info)

    #Combining intron/exon data and read locations into a single dataframe
    combined_dt_intron <- cbind(combined_dt_intron_info, combined_dt_intron_data)
    combined_dt_intron$event[which(combined_dt_intron$introns != "" & combined_dt_intron$SJ_IR == "IR")] <- combined_dt_intron$introns[which(combined_dt_intron$introns != "" & combined_dt_intron$SJ_IR == "IR")]
    combined_dt_intron$assembly <- unique(samples$assembly)

    message("\t","Data combined.dt.","\n")
    return(combined_dt_intron)
}


#The code from here onwards needs to be incorporated into a for loop for each proband/family
compare.samples <- function(samples, combined_dt_intron, exportlocation){
    message("Comparing samples...")
    for(sample in seq(1,nrow(samples))){
        combined_dt_intron_test <- combined_dt_intron
        if(samples$sampletype[sample] == 'test'){
            family <- samples$sampleID[which(samples$family == samples$family[sample])]
            message("\t",family)
            familycols <- paste0("sj_pct_", family)
            familyreadcols <- paste0("sj_", family)

            ctrls <- samples$sampleID[which(samples$family != samples$family[sample])]
            ctrlscols <- paste0("sj_pct_", ctrls)
            ctrlsreadcols <- paste0("sj_", ctrls)

            proband <- samples$sampleID[sample]

            combined_dt_intron_test$proband <- samples$sampleID[sample]

            combined_dt_intron_test$controlavg <- rowMeans(combined_dt_intron_test[,ctrlscols, with=F])
            combined_dt_intron_test$controlavgreads <- rowMeans(combined_dt_intron_test[,ctrlsreadcols, with=F])
            combined_dt_intron_test <- cbind(combined_dt_intron_test, controlsd = apply(combined_dt_intron_test[,ctrlscols, with=F], 1, sd))

            combined_dt_intron_test$difference <- combined_dt_intron_test[,paste0("sj_pct_", proband), with=F] - combined_dt_intron_test$controlavg
            combined_dt_intron_test$controln <- length(ctrlscols)
            combined_dt_intron_test$two_sd <- abs(combined_dt_intron_test$difference) > combined_dt_intron_test$controlsd*2
            combined_dt_intron_test$three_sd <- abs(combined_dt_intron_test$difference) > combined_dt_intron_test$controlsd*3
            combined_dt_intron_test$four_sd <- abs(combined_dt_intron_test$difference) > combined_dt_intron_test$controlsd*4

            for (i in seq(1,nrow(combined_dt_intron_test))){
                event_unique_count <- 0
                for (member in familycols){
                    if (combined_dt_intron_test$controlavg[i] == 0 & combined_dt_intron_test[,member, with=F][i] != 0){
                        event_unique_count <- event_unique_count + 1
                    }
                }
                if (event_unique_count >= 1){
                    combined_dt_intron_test$unique[i] <- paste(event_unique_count, "/", length(familycols), sep="")
                }else{
                    combined_dt_intron_test$unique[i] <- ""
                }
            }

        #Extracting the columns in the required order
        combined_dt_intron_final <- combined_dt_intron_test[,c("assembly",infocols[-which(infocols=="introns")],"SJ_IR","frame_conserved","unique","proband","difference",familycols,"controlavg","controlsd","controln", familyreadcols, "controlavgreads","two_sd","three_sd","four_sd"), with=F]

        #Sort by the greatest difference in percentage
        combined_dt_intron_final <- combined_dt_intron_final[order(abs(combined_dt_intron_final$difference), decreasing = T),]

        #Exclude Novel Junction not mapping to a known SJ and any not-annotated IR reads
        combined_dt_intron_final <- combined_dt_intron_final[!which(combined_dt_intron_final$SJ_IR == "IR" & combined_dt_intron_final$annotated == "N")]
        combined_dt_intron_final <- combined_dt_intron_final[!which(combined_dt_intron_final$event == "unannotated junctions")]

        report <- T

        testgenes <- unique(samples$genes[sample])

        if(report == TRUE){
            #for(i in seq(0,length(testgenes))){
                #combined_dt_intron_final <- combined_dt_intron_final[which(combined_dt_intron_final$genes == testgenes[i])]
                #combined_dt_intron_final <- combined_dt_intron_final[which(combined_dt_intron_final$two_sd == TRUE)]
                generate.report(combined_dt_intron_final[which(combined_dt_intron_final$genes == testgenes)], length(familycols), testgenes, exportlocation, samples$sampleID[sample])
            #}
        }else{
            #Exporting the combined.dt dataframe to an excel spreadsheet
            openxlsx::write.xlsx(combined_dt_intron_final, paste(exportlocation,samples$sampleID[sample],"_combined_dt_",samples$assembly[1],".xlsx", sep=""), asTable = T, overwrite = T)
        }
    }
}
    message("\t","\t","Reports generated...")
}


#--Optional - add variants and allele bias-------------------------------------
add.VCF <- function(samples){
    for (row_number in 1:nrow(samples)) {
        sampleID <- samples[row_number, sampleID]
        vcffile <- samples[row_number, VCF]
        message(sampleID)
        vcf <- fread(file=vcffile, sep ="\t", fill=T)
        names(vcf) <- c("chr","start","end","alt","genotype","AB")
        vcf$end <- vcf$start+1
        vcfGR <- makeGRangesFromDataFrame(vcf[,c("chr","start","end","genotype","AB"), with=F], keep.extra.columns = T)
        head(vcfGR)
    }
}







#--Other-----------------------------------------------------------------------
GRanges.to.SAF <- function(gr, minAnchor=1){
    data.table(
        GeneID  = seq_along(gr),
        Chr     = as.factor(seqnames(gr)),
        Start   = start(gr) - (minAnchor - 1),
        End     = end(gr) + (minAnchor - 1),
        Strand  = as.factor(strand(gr))
    )
}

#Splicing junctions
gen.exon.range <- function(strand, exon_range_start, exon_range_end) {

    exon_range_start <- as.numeric(strsplit(exon_range_start,split=" ")[[1]][3])
    exon_range_end <- as.numeric(strsplit(exon_range_end,split=" ")[[1]][3])
    exon_ranges <- c(exon_range_start, exon_range_end)

    if (!is.na(exon_range_start) & !is.na(exon_range_end)) {
        return(paste("exon ",min(exon_ranges)," ~ ","exon ", max(exon_ranges)+1, sep="", collapse=""))

    }else if (!is.na(exon_range_start)){
        if (strand == "-"){
            return(paste("cryptic donor", " ~ ", "exon ", exon_range_start+1, sep="", collapse=""))
        }
        if (strand == "+"){
            return(paste("exon ",exon_range_start, " ~ ", "cryptic acceptor", sep="", collapse=""))
        }


    }else if (!is.na(exon_range_end)){
        if (strand == "-"){
            return(paste("exon ", exon_range_end, " ~ ", "cryptic acceptor", sep="", collapse=""))
        }
        if (strand == "+"){
            return(paste("cryptic donor", " ~ ", "exon ", exon_range_end+1, sep="", collapse=""))
        }
    }else{
        return("unannotated junctions")
    }
}

#Introns
gen.introns <- function(exon_range_start, exon_range_end) {

    exon_range_start <- as.numeric(strsplit(exon_range_start,split=" ")[[1]][3])
    exon_range_end <- as.numeric(strsplit(exon_range_end,split=" ")[[1]][3])
    exon_ranges <- c(exon_range_start, exon_range_end)

    if (!is.na(exon_range_start) & !is.na(exon_range_end)) {
        return(paste("intron ",min(exon_ranges), sep="", collapse=""))

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
    rfsq <- refseq_introns_exons_of_interest
    frame <- c()

    for(i in seq(1,length(x))){
        rg_event <- event[i]
        rg_start <- x[i]
        rg_end <- y[i]
        if(event[i] == "SJ"){
            #is the start/end annotated
            #both annotated
            if(rg_start %in% rfsq$region_start & rg_end %in% rfsq$region_end){
                pairstart <- rfsq$region_start[which(rfsq$region_end == rg_end)]
                pairend <- rfsq$region_end[which(rfsq$region_start == rg_start)]
                if(rg_end == pairend){
                    frame[i] <- TRUE
                } else{
                    dist2authentic <- abs(rg_end-pairend)
                    frame[i] <- dist2authentic%%3 == 0
                }
            }
            #start annotated
            else if(rg_start %in% rfsq$region_start){
                pairend <- rfsq$region_end[which(rfsq$region_start == rg_start)]
                dist2authentic <- abs(rg_end-pairend)
                frame[i] <- dist2authentic%%3 == 0
            }
            #end annotated
            else if(rg_end %in% rfsq$region_end){
                pairstart <- rfsq$region_start[which(rfsq$region_end == rg_end)]
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
    headerStyle <- createStyle(textRotation = 45, fgFill = "#4F81BD", textDecoration = "bold", fontColour = "white")

    # Add data to the worksheet
    writeDataTable(wb, sht, data)

    # Add the percent, centre, and header styles to the desired cells
    addStyle(wb, sht, style=pct, cols = c(14:(17+familymembers-1)), rows = 2:(nrow(data)+1), gridExpand=TRUE) #note this may change with additional family members
    addStyle(wb, sht, style=centre, cols = c(1,6:8,10:12), rows = 2:(nrow(data)+1), gridExpand=TRUE)
    addStyle(wb, sht, headerStyle, cols = 1:(ncol(data)+1), rows = 1, gridExpand = TRUE)
    addStyle(wb, sht, twodp, cols = 22, rows = 2:(nrow(data)+1), gridExpand = TRUE)

    # Set column widths for event and proband - Currently not working 20211013
    #width_vec <- apply(data, 2, function(x) max(nchar(as.character(x)) + 2, na.rm = TRUE))
    #setColWidths(wb, sht, cols = c(9,12), widths = width_vec[c(9,12)])

    # Add conditionalFormatting to difference and percentage columns
    conditionalFormatting(wb, sht, cols=c(15:(17+familymembers-1)), rows = 2:(nrow(data)+1),
                          rule = NULL, style = c("#FCFCFF","#63BE7B"), type = "colourScale")

    conditionalFormatting(wb, sht, cols=c(14), rows = 2:(nrow(data)+1),
                          rule = NULL, style = c("#F8696B","#FCFCFF","#F8696B"), type = "colourScale")

    # Export report
    saveWorkbook(wb, paste(export,"/",sample,"_",gene,"_combined_dt_",".xlsx", sep=""), overwrite = T)
}

#Simple code for the opposite of %in%
`%nin%` <- Negate(`%in%`)
