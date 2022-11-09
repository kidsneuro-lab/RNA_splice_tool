#' Run cortar
#'
#' `cortar()` runs the entire cortar pipeline with specified parameters,
#' and returns excel & pdf reports.
#'
#' @param file A file path, pointing to the cortar samplefile - see readme for
#'     more information
#' @param mode One of `"default"` (default),`"panel"`, or `"research"` - see
#'     readme for more information
#' @param assembly Assembly used for alignment: either `"hg38"` (default) or
#'     `"hg19"`
#' @param annotation  Annotation used for alignment: either `"UCSC"` or
#'     `"1000genomes"`
#' @param paired Is the RNA-seq paired-end?: `TRUE`/`FALSE`
#' @param stranded Strandedness of the RNA-seq: `0` for unstranded, `1`
#'     for forward stranded or `2` for reverse stranded
#' @param subset Does the RNA-Seq need to be subsetted to the genes
#'     of interest? `TRUE`/`FALSE` (Optional, but improves speed of subsequent
#'     runs.)
#' @param output_dir A directory path, pointing to the desired location for
#'     export of cortar results (e.g. `"output/"`)
#' @param genelist A character vector with genes/RefSeq transcripts of interest
#'     (Only for panel or research modes; default = NULL)
#' @param testsamples A character vector with sampleIDs of test samples for
#'     gene panel testing (Only for panel mode; default = NULL)
#'
#' @export
#'
#' @examples
#' #### == COMING SOON == ####
#'

cortar <- function(file,
                   mode = "default",
                   assembly = "hg38",
                   annotation = "UCSC",
                   paired = TRUE,
                   stranded = 2,
                   subset = FALSE,
                   output_dir = "output/test",
                   genelist = NULL,
                   testsamples = NULL) {

 #Error catching
  if(!file.exists(file)){
    stop("File '", file, "' does not exist or is non-readable. getwd()=='",
         getwd(), "'")
  }

  if(!dir.exists(output_dir)){
    message("Directory '", output_dir, "' does not exist.
    Would you like to create this directory?
    \t 1. Yes
    \t 2. No"
            )
    selection <- readline(prompt="Selection: ")
    if(selection == "1"){
      dir.create(output_dir)
      message("Directory '", output_dir, "' created.")
      message("")
    }else{
      stop("Directory '", output_dir, "' does not exist. getwd()=='",
         getwd(), "'")
    }
  }

  if(mode != "default"){
        if(mode == "panel" | mode == "research"){
        }else{
            stop("'",mode,"' is not an available cortar mode ('default','panel','research')")
        }
  }

  if(assembly %nin% c("hg19","hg38")){
    stop("Assembly '", assembly, "' is not an available assembly ('hg38','hg19')")
  }

  if(annotation %nin% c("1000genomes","UCSC")){
    stop("Annotation '", annotation, "' is not an available annotation ('1000genomes','UCSC')")
  }

  if(class(paired) != "logical"){
    stop("Paired must be a logical: TRUE, FALSE.
         Supplied:", paired,"
         Note: You may need to remove quotation marks.")
  }


  message(paste0("Running cortar "))
  message(paste0("        file: ", file))
  message(paste0("        mode: ", mode))
  message(paste0("    assembly: ", assembly))
  message(paste0("  annotation: ", annotation))
  message(paste0("      paired: ", paired))
  message(paste0("    stranded: ", stranded))
  message(paste0("      output: ", output_dir))
  message("")


  file <- data.table::fread(file)

  if(mode == "panel" | mode == "research"){
      genes_tx <- selectGenesTranscripts(
        genes = genelist,
        assembly = assembly,
        annotation = annotation
      )
    }

  if(mode == "default"){
      genes_tx <- selectGenesTranscripts(
        genes = file$genes,
        assembly = assembly,
        annotation = annotation
      )
    }

  if (subset == T) {
    subsetBamfiles
  }

  reads <- extractCountReads(
    genes.GRanges = genes_tx[[1]],
    introns.GRanges = genes_tx[[2]][[1]],
    intron_starts.GRanges = genes_tx[[4]][[1]],
    intron_ends.GRanges = genes_tx[[4]][[2]],
    bamfiles = file$bamfile,
    sample_names = file$sampleID,
    assembly = assembly,
    paired = T,
    stranded = 2
  )

  events <- annotateQuantifyEvents(
    ids = file$sampleID,
    combined_sj = reads,
    introns.GRanges = genes_tx[[2]][[1]],
    introns_other_tx.GRanges = genes_tx[[3]],
    introns = genes_tx[[2]][[2]],
    assembly = assembly
  )

  comparisons <- compareSplicing(
    all_splicing_events = events,
    Sample_File = file,
    mode = mode
  )

  generateReport(
    comparisons = comparisons,
    Sample_File = file,
    Export = output_dir,
    mode = mode
  )


  message("")
  message(paste("Done! Reports saved in:", output_dir))

}
