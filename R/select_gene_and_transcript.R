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
tx_extraction <- function(genes) {
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

gene_to_GRange <- function(gene_tx, assembly) {
  genes <- Ensembl_Genes[`Gene name` %in% gene_tx$gene_name]

  genes.GRanges <- GenomicRanges::GRanges(
    seqnames = genes$`Chromosome/scaffold name`,
    IRanges::IRanges(
      start = genes$`Gene start (bp)`,
      end = genes$`Gene end (bp)`
    ),
    strand = genes$Strand
  )

  if (assembly == "hg38") {
    GenomeInfoDb::seqlevelsStyle(genes.GRanges) <- "UCSC"
  }
  return(genes.GRanges)
}


introns_to_GRange <- function(gene_tx, assembly) {
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

  if (assembly == "hg38") {
      GenomeInfoDb::seqlevelsStyle(introns.GRanges) <- "UCSC"
  }
  return(introns.GRanges)
}

# Select introns of interest for other transcripts and create a GRanges object
Introns_Other_Tx <- Refseq_Genes[gene_name %in% unlist(genes) &
  tx_id %nin% Sample_File$transcript &
  region_type == c("intron")]

introns_other_tx.GRanges <- GRanges(
  seqnames = Introns_Other_Tx$chrom,
  IRanges(
    start = Introns_Other_Tx$region_start,
    end = Introns_Other_Tx$region_end
  ),
  strand = Introns_Other_Tx$strand
)

mcols(introns_other_tx.GRanges)["tx_id"] <- Introns_Other_Tx$tx_version_id
mcols(introns_other_tx.GRanges)["intron_no"] <- Introns_Other_Tx$region_no
mcols(introns_other_tx.GRanges)["gene"] <- Introns_Other_Tx$gene_name

seqlevelsStyle(introns_other_tx.GRanges) <- "UCSC"

# Create a GRanges object for intron-exon junctions
intron_starts.GRanges <- GRanges(
  seqnames = Introns$chrom,
  IRanges(
    start = Introns$region_start - 4,
    end = Introns$region_start + 3
  ),
  strand = Introns$strand
)

intron_ends.GRanges <- GRanges(
  seqnames = Introns$chrom,
  IRanges(
    start = Introns$region_end - 3,
    end = Introns$region_end + 4
  ),
  strand = Introns$strand
)

seqlevelsStyle(intron_starts.GRanges) <- "UCSC"
seqlevelsStyle(intron_ends.GRanges) <- "UCSC"
