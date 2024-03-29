% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_select_gene_and_transcript.R
\name{tx_extraction}
\alias{tx_extraction}
\title{Extract Genes & Transcripts}
\usage{
tx_extraction(genes, refseq_assembly = refseq_introns_exons_hg38)
}
\arguments{
\item{genes}{A character vector of gene names and/or RefSeq transcript IDs}

\item{refseq_assembly}{A data.table of RefSeq coding genes. Either:
\code{refseq_introns_exons_hg38} (default) or \code{refseq_introns_exons_hg19}}
}
\value{
A data.table with two columns \code{gene_name} and \code{tx}, with a
row for each gene/transcript provided, containing:
\itemize{
\item the query gene name and the canonical RefSeq transcript - if a gene name
was provided
\item the gene name and the query transcript - if a RefSeq transcript was
provided
}
}
\description{
\code{tx_extraction} extracts genes and/or transcripts from a list and fetches
the corresponding gene and/or transcript as required.
}
\examples{
#### == COMING SOON == ####

}
