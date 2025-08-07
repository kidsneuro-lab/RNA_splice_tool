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
                                   annotation,
                                   debug = "") {
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

  gene_tx <- tx_extraction(genes, Refseq_Genes, debug)
  genes.GRanges <- gene_to_GRange(gene_tx, assembly, annotation, Refseq_Genes, Ensembl_Genes, debug)
  introns.GRanges <- introns_to_GRange(gene_tx, assembly, annotation, Refseq_Genes, debug)
  introns_other_tx.GRanges <- introns_other_tx_to_GRange(genes, gene_tx, assembly, annotation, Refseq_Genes, debug)
  introns_jx.GRanges <- introns_jx_to_GRange(gene_tx, assembly, annotation, Refseq_Genes, debug)

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
                          refseq_assembly = refseq_introns_exons_hg38,
                          debug = "") {
  genes_tx <- data.table("gene_name" = character(), "tx" = character())

  for (gene in seq(1, length(genes))) {
    if (length(grep("NM_[0-9]+\\.[0-9]+", genes[gene])) > 0) {
      tx <- stringr::str_extract_all(genes[gene], "NM_[0-9]+\\.[0-9]+")[[1]][1]
      if(tx %in% refseq_assembly[,tx_version_id]){
        gene_name <- unique(refseq_assembly[tx_version_id == tx & region_type == c("intron"), "gene_name"])[[1]]
      } else {
        stop(paste0("Transcript identifier `",tx,"` is invalid"))
      }
    } else if (length(grep("[A-Za-z0-9]{3,10}", genes[gene])) > 0){
      gene_name <- genes[gene]
      if(genes[gene] %in% refseq_assembly[,gene_name]){
        tx <- unique(refseq_assembly[gene_name == genes[gene] & canonical == 1, "tx_version_id"])[[1]]
      } else {
        stop(paste0("Gene name `",genes[gene],"` is invalid"))
      }
    } else {
      gene_name <- NULL
      tx <- NULL
    }
    genes_tx <- rbind(genes_tx, data.table("gene_name" = gene_name, "tx" = tx))
  }

  if(debug != "" | debug == FALSE){
    fwrite(as.data.table(genes_tx),paste0(debug,"/","1_tx_extraction.tsv"), sep = "\t")
  }

  return(as.data.table(genes_tx))
}

gene_to_GRange <- function(gene_tx, assembly, annotation, Refseq_Genes, Ensembl_Genes, debug = "") {
  genes <- Ensembl_Genes[`Gene name` %in% gene_tx$gene_name]
  rs_genes <- unique(Refseq_Genes[gene_name %in% gene_tx$gene_name, .(gene_name)])

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

  # output_dir
  if (nrow(genes) != nrow(rs_genes)) {
    message("Gene name invalid.
    Would you like to continue? (Note: This is unlikely to work)
    \t 1. Yes
    \t 2. No")
    print(genes)
    print(rs_genes)
    selection <- readline(prompt = "Selection: ")
    if (selection %in% c("1", "Yes", "Y", "yes", "y")) {
    } else {
      stop(
        "Gene name invalid"
      )
    }
  }

  if(debug != "" | debug == FALSE){
    fwrite(as.data.table(genes.GRanges),paste0(debug,"/","2_gene_to_GRange.tsv"), sep = "\t")
  }

  return(genes.GRanges)
}

introns_to_GRange <- function(gene_tx, assembly, annotation, Refseq_Genes, debug = "") {
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

  if(debug != "" | debug == FALSE){
    fwrite(as.data.table(introns.GRanges),paste0(debug,"/","3_introns_to_GRange_GRanges.tsv"), sep = "\t")
    fwrite(as.data.table(introns),paste0(debug,"/","3_introns_to_GRange_introns.tsv"), sep = "\t")
  }

  return(list(introns.GRanges, introns))
}

introns_other_tx_to_GRange <- function(genes, gene_tx, assembly, annotation, Refseq_Genes, debug = "") {
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

  if(debug != "" | debug == FALSE){
    fwrite(as.data.table(introns.GRanges),paste0(debug,"/","4_introns_other_tx_to_GRange.tsv"), sep = "\t")
  }

  return(introns.GRanges)
}

introns_jx_to_GRange <- function(gene_tx, assembly, annotation, Refseq_Genes, debug = "") {
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

  if(debug != "" | debug == FALSE){
    fwrite(as.data.table(intron_starts.GRanges),paste0(debug,"/","5_introns_jx_to_GRange_starts.tsv"), sep = "\t")
    fwrite(as.data.table(intron_ends.GRanges),paste0(debug,"/","5_introns_jx_to_GRange_ends.tsv"), sep = "\t")
  }

  return(list(intron_starts.GRanges, intron_ends.GRanges))
}

`%nin%` <- Negate(`%in%`)
