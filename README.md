
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cortar 

<!-- badges: start -->
<!-- badges: end -->

The goal of cortar is to extract reads at the exon-intron junction from RNA-seq
and report the proportion of splicing events across each intron. By comparing
these values to controls, deviations from normal splicing can be characterised
in test samples.

## Installation

You can install the development version of cortar from
[GitHub](https://github.com/) with: 

``` r
# install.packages("devtools")
devtools::install_github("kidsneuro-lab/RNA_splice_tool")
```

## Usage
To use cortar, a samplefile needs to be created for each run. This file
contains six columns and a row for each sample (as shown below):
* sampleID: A unique identifier for the sample
* familyID: A unique identifier for related samples that should not be compared
to one another
* sampletype: Whether analysis and report is desired for this sample ("test")
* genes: Gene symbol for gene under investigation (can be blank for
panel/research mode)
* transcript: Transcript under investigation (blank transcript will default to
the canonical RefSeq transcript)
* bamfile: Absolute or relative file path to the processed bamfile of the sample

``` r
#> sampleID    familyID   sampletype   genes      transcript       bamfile
#> ---------   --------   ----------   ------     --------------   --------------------------------
#> proband_1   1	  test         DMD        NM_004006.3	   Z:/path/to/bamfile/proband_1.bam
#> proband_2   2	  test         TTN	  NM_001267550.2   Z:/path/to/bamfile/proband_2.bam
#> proband_3   3	  test         CFTR	  NM_000492.4	   Z:/path/to/bamfile/proband_3.bam
#> proband_4   4	  test         NF1	  NM_001042492.3   Z:/path/to/bamfile/proband_4.bam
#> mother_4    4		       NF1	  NM_001042492.4   Z:/path/to/bamfile/mother_4.bam
#> proband_5   5	  test         COL2A1	  NM_001844.5	   Z:/path/to/bamfile/proband_5.bam
```

