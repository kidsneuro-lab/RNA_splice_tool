#' Split a string
#'
#' @param string A character vector with, at most, one element.
#' @inheritParams stringr::str_split
#'
#' @return A character vector.
#' @export
#'
#' @examples
#' x <- "alfa,bravo,charlie,delta"
#' str_split_one(x, pattern = ",")
#' str_split_one(x, pattern = ",", n = 2)
#'
#' y <- "192.168.0.1"
#' str_split_one(y, pattern = stringr::fixed("."))
#'
# INPUT = Vector of gene names and/or transcripts
# OUTPUT = Subsetted Refseq table and GRanges object
#'
selectGenesTranscripts <- function(genes, assembly, annotation) {


message("Selecting genes and transcripts...")
message("")

  if (assembly == "hg38") {
      Refseq_Genes <- refseq_introns_exons_hg38
      Ensembl_Genes <- ensembl_allgenes_chr1_Y_hg38
  } else if (assembly == "hg19") {
      Refseq_Genes <- refseq_introns_exons_hg19.tsv
      Ensembl_Genes <- ensembl_allgenes_chr1_Y_hg19
  }

  gene_tx <- tx_extraction(genes, Refseq_Genes)
  genes.GRanges <- gene_to_GRange(gene_tx, assembly, annotation, Refseq_Genes, Ensembl_Genes)
  introns.GRanges <- introns_to_GRange(gene_tx, assembly, annotation, Refseq_Genes)
  introns_other_tx.GRanges <- introns_other_tx_to_GRange(genes, gene_tx, assembly, annotation, Refseq_Genes)
  introns_jx.GRanges <- introns_jx_to_GRange(gene_tx, assembly, annotation, Refseq_Genes)

  return(list(
    genes.GRanges,
    introns.GRanges,
    introns_other_tx.GRanges,
    unlist(introns_jx.GRanges)
  ))
}

# add error - transcript/gene could not be found and then return suggestions using grep
tx_extraction <- function(genes, Refseq_Genes) {
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
  return(as.data.table(genes_tx))
}

gene_to_GRange <- function(gene_tx, assembly, annotation, Refseq_Genes, Ensembl_Genes) {
  genes <- Ensembl_Genes[`Gene name` %in% gene_tx$gene_name]

  genes.GRanges <- GenomicRanges::GRanges(
    seqnames = genes$`Chromosome/scaffold name`,
    IRanges::IRanges(
      start = genes$`Gene start (bp)`,
      end = genes$`Gene end (bp)`
    ),
    strand = genes$Strand
  )

  if (annotation == "UCSC") {
    GenomeInfoDb::seqlevelsStyle(genes.GRanges) <- "UCSC"
  }
  return(genes.GRanges)
}

introns_to_GRange <- function(gene_tx, assembly, annotation, Refseq_Genes) {

  introns <- Refseq_Genes[tx_version_id %in% gene_tx$tx & region_type == c("intron")]

  introns.GRanges <- GenomicRanges::GRanges(
    seqnames = introns$chrom,
    IRanges::IRanges(
      start = introns$region_start,
      end = introns$region_end
    ),
    strand = introns$strand
  )

  GenomicRanges::mcols(introns.GRanges)["intron_no"] <- introns$region_no
  GenomicRanges::mcols(introns.GRanges)["gene"] <- introns$gene_name

  if (annotation == "UCSC") {
    GenomeInfoDb::seqlevelsStyle(introns.GRanges) <- "UCSC"
  }
  return(list(introns.GRanges, introns))
}

introns_other_tx_to_GRange <- function(genes, gene_tx, assembly, annotation, Refseq_Genes) {
  introns <- Refseq_Genes[
    gene_name %in% unlist(genes) &
      tx_version_id %nin% gene_tx$tx &
      region_type == c("intron")
  ]

  introns.GRanges <- GenomicRanges::GRanges(
    seqnames = introns$chrom,
    IRanges::IRanges(
      start = introns$region_start,
      end = introns$region_end
    ),
    strand = introns$strand
  )

  GenomicRanges::mcols(introns.GRanges)["tx_id"] <- introns$tx_version_id
  GenomicRanges::mcols(introns.GRanges)["intron_no"] <- introns$region_no
  GenomicRanges::mcols(introns.GRanges)["gene"] <- introns$gene_name

  if (annotation == "UCSC") {
    GenomeInfoDb::seqlevelsStyle(introns.GRanges) <- "UCSC"
  }
  return(introns.GRanges)
}

introns_jx_to_GRange <- function(gene_tx, assembly, annotation, Refseq_Genes) {
  introns <- Refseq_Genes[tx_version_id %in% gene_tx$tx & region_type == c("intron")]

  intron_starts.GRanges <- GenomicRanges::GRanges(
    seqnames = introns$chrom,
    IRanges::IRanges(
      start = introns$region_start - 4,
      end = introns$region_start + 3
    ),
    strand = introns$strand
  )

  intron_ends.GRanges <- GenomicRanges::GRanges(
    seqnames = introns$chrom,
    IRanges::IRanges(
      start = introns$region_end - 3,
      end = introns$region_end + 4
    ),
    strand = introns$strand
  )

  if (annotation == "UCSC") {
    GenomeInfoDb::seqlevelsStyle(intron_starts.GRanges) <- "UCSC"
    GenomeInfoDb::seqlevelsStyle(intron_ends.GRanges) <- "UCSC"
  }

  return(list(intron_starts.GRanges, intron_ends.GRanges))
}

`%nin%` <- Negate(`%in%`)
