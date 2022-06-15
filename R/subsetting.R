library(data.table)

sampledatabase <- fread("RNAseq_sample_database.csv")

forsubset <- sampledatabase[tissue == "muscle" & batch == "2021-2"]

forsubset$bamfile_full_hpc

subsetgenes <- Ensembl_Genes[`Gene name` %in% forsubset$gene]
#forsubset$gene %nin% subsetgenes$`Gene name`

subsetgenes2 <- paste0("\'","\'","chr",subsetgenes$`Chromosome/scaffold name`,":",subsetgenes$`Gene start (bp)`-1000,"-",
       subsetgenes$`Gene end (bp)`+1000,"\'","\'")

subsetgenes3 <- paste(subsetgenes2[1],subsetgenes2[2])

fwrite(as.list(forsubset$bamfile_full_hpc), "subset_bamfiles.tsv", sep = "\t")
fwrite(as.list(subsetgenes2), "subset_genes.tsv", sep = "\t")

fs_samplefile <- forsubset[,c(1,2,3,4,5,6,7)]
names(fs_samplefile)[7] <- "bamfile"
fwrite(fs_samplefile,"batch_samplefile.tsv", sep = "\t")


fs_samplefile <- forsubset[,c(1,2,3,4,5,6,20)]
names(fs_samplefile)[7] <- "bamfile"
fwrite(fs_samplefile,"tissue_samplefile.tsv", sep = "\t")
