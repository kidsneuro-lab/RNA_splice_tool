#' Extract gene coordinates for RNA-seq subsetting
#'
#' Returns the gene coordinates for a given entrez_gene_symbol or
#' RefSeq/Ensembl transcript to be used for RNA-seq subsetting
#'
#' @param genes A character vector with the entrez_gene_symbols,
#' RefSeq transcript ids, or Ensembl gene or transcript ids for the genes for
#' analysis.
#' @param hg Either 19 or 38
#' @param overhang Number of nucleotides flanking the gene to be included
#' (default: 1000nt)
#'
#' @return A character vector.
#' @export
#'
#' @examples
#' == COMING SOON==
#'

subsetBamfiles <- function(genes, hg, overhang = 1000){

    if (hg == 38) {
        Refseq_Genes <- refseq_introns_exons_hg38
        Ensembl_Genes <- ensembl_allgenes_chr1_Y_hg38
    } else if (hg == 19) {
        Refseq_Genes <- refseq_introns_exons_hg19
        Ensembl_Genes <- ensembl_allgenes_chr1_Y_hg19
    }

    genes_tx <- data.table("gene_name" = character(), "tx" = character())

    for (gene in seq(1, length(genes))) {
        if (length(grep("NM_[0-9]+\\.[0-9]+", genes[gene])) > 0) {
            tx <- stringr::str_extract_all(genes[gene], "NM_[0-9]+\\.[0-9]+")[[1]][1]
            gene_name <- unique(Refseq_Genes[tx_version_id == tx & region_type == c("intron"), "gene_name"])[[1]]
        } else {
            gene_name <- genes[gene]
            tx <- unique(Refseq_Genes[gene_name == genes[gene] & canonical == 1, "tx_version_id"])[[1]]
        }
        genes_tx <- rbind(genes_tx, data.table("gene_name" = gene_name, "tx" = tx))
    }

    subsetgenes <- Ensembl_Genes[`Gene name` %in% genes_tx$gene_name]
    forsubset <- paste0("\'","\'","chr",subsetgenes$`Chromosome/scaffold name`,":",subsetgenes$`Gene start (bp)`-overhang,"-",
                           subsetgenes$`Gene end (bp)`+overhang,"\'","\'")
    forsubset <- paste(forsubset,collapse=" ")
    return(forsubset)
}
