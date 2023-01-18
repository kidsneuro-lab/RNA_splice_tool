#' Select Genes & Transcripts
#'
#' `selectGenesTranscripts` extracts genes and/or transcripts from the cortar
#' samplefile and returns gene and intron GRange objects required for further
#' analysis.
#'
#' @param genes A character vector of gene names and/or RefSeq transcript IDs
#'     with version numbers (e.g. NM_004006.2)
#' @param assembly Assembly used for alignment: either `"hg38"` or `"hg19"`
#' @param annotation  Annotation used for alignment: either `"UCSC"`or
#'     `"1000genomes"`
#'
#' @returns A list of GRanges for:
#' * genes of interest
#' * introns in genes of interest
#' * introns in genes of interest (other transcripts)
#' * exon-intron junctions of genes of interest
#' @export
#'
#' @examples
#' #### == COMING SOON == ####
#'
selectGenesTranscripts <- function(genes,
                                   assembly,
                                   annotation) {
  # Initialisation message
  message("Selecting genes and transcripts...")
  message("")

  if (assembly == "hg38") {
    Refseq_Genes <- refseq_introns_exons_hg38
    Ensembl_Genes <- ensembl_allgenes_chr1_Y_hg38
  } else if (assembly == "hg19") {
    Refseq_Genes <- refseq_introns_exons_hg19
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

#' Extract Genes & Transcripts
#'
#' `tx_extraction` extracts genes and/or transcripts from a list and fetches
#' the corresponding gene and/or transcript as required.
#'
#' @param genes A character vector of gene names and/or RefSeq transcript IDs
#' @param refseq_assembly A data.table of RefSeq coding genes. Either:
#'     `refseq_introns_exons_hg38` (default) or `refseq_introns_exons_hg19`
#'
#' @returns A data.table with two columns `gene_name` and `tx`, with a
#' row for each gene/transcript provided, containing:
#' * the query gene name and the canonical RefSeq transcript - if a gene name
#' was provided
#' * the gene name and the query transcript - if a RefSeq transcript was
#' provided
#' @export
#'
#' @examples
#' #### == COMING SOON == ####
#'
tx_extraction <- function(genes,
                          refseq_assembly = refseq_introns_exons_hg38) {
  genes_tx <- data.table("gene_name" = character(), "tx" = character())

  for (gene in seq(1, length(genes))) {
    if (length(grep("NM_[0-9]+\\.[0-9]+", genes[gene])) > 0) {
      tx <- stringr::str_extract_all(genes[gene], "NM_[0-9]+\\.[0-9]+")[[1]][1]
      gene_name <- unique(refseq_assembly[tx_version_id == tx & region_type == c("intron"), "gene_name"])[[1]]
    } else {
      gene_name <- genes[gene]
      tx <- unique(refseq_assembly[gene_name == genes[gene] & canonical == 1, "tx_version_id"])[[1]]
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
