#' Run the noisyR pipeline for the transcript approach
#' @description Calls the functions to run each of the three steps of the pipeline
#' (similarity calculation, noise quantification, noise removal), with the specified parameters.
#' See the individual function documentation for more details and required arguments.
#' Required steps: \code{\link{calculate_expression_similarity_transcript}},
#' \code{\link{calculate_noise_threshold}}. \code{\link{remove_noise_from_bams}}.
#' Optional steps: \code{\link{calculate_noise_threshold_method_statistics}}
#' @param bams,path.bams either a path to the directory where the BAM files are
#' or a vector of paths to each individual file; if a path is specified,
#' it extracts all files that end in .bam; looks in the working directory by default
#' @param genes a tibble of the exons extracted from the gtf file;
#' this is meant for speed if the output of \code{\link{cast_gtf_to_genes}} is already generated,
#' or if the user wants to only calculate similarity for a subset of exons
#' @param path.gtf the path to the gtf/gff annotation file (only used if genes is not
#' provided); if unspecified, looks for one in the working directory
#'@param ncores Number of cores for parallel computation; defaults to sequential computation,
#' but parallelisation is highly encouraged; it is set to detectCores() if higher
#' @param similarity.threshold,method.chosen parameters passed on to
#' \code{\link{calculate_noise_threshold}}; they can be single values or vectors;
#' if they are vectors optimal values are computed by calling
#' \code{\link{calculate_noise_threshold_method_statistics}} and
#' minimising the coefficient of variation across samples; all possible values for
#' method.chosen can be viewed by \code{\link{get_methods_calculate_noise_threshold}};
#' only boxplot based methods are accepted for the transcript approach due to the number of observations
#' and high variance
#' @param ... arguments to be passed on to individual pipeline steps; see their documentation
#' for more details and required arguments
#' @return The denoised BAM files are created, as specified by the destination.files argument
#' of remove_noise_from_bams()
#' @seealso \code{\link{noisyr}}, \code{\link{noisyr_counts}}
#' @export
#' @examples
#' bams <- rep(system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE), 2)
#' genes <- data.frame("id" = 1:2,
#'                     "gene_id" = c("gene1", "gene2"),
#'                     "seqid" = c("seq1", "seq2"),
#'                     "start" = 1,
#'                     "end" = 1600)
#' noisyr_transcript(
#'   bams = bams,
#'   genes = genes,
#'   destination.files = paste0(tempdir(), "/", basename(bams), ".noisefiltered.bam")
#' )
noisyr_transcript = function(
  bams = NULL,
  path.bams = ".",
  genes = NULL,
  path.gtf = list.files(".", pattern="\\.g[tf]f$"),
  ncores = 1,
  similarity.threshold = 0.25,
  method.chosen = "Boxplot-IQR",
  ...
){
  base::message(">>> noisyR transcript approach pipeline <<<")

  if(base::is.null(bams)){
    bams <- base::list.files(path.bams, pattern="\\.bam$", full.names=TRUE)
  }
  if(base::length(bams)<2) base::stop("Please provide at least 2 BAM files")

  if(base::is.null(genes)){
    base::message("Creating gene table from gtf file...")
    genes <- noisyr::cast_gtf_to_genes(path.gtf, ...)
  }
  if(base::nrow(genes) < 2) base::stop("Please provide at least 2 genes")

  expression.summary <- noisyr::calculate_expression_similarity_transcript(
    bams = bams,
    genes = genes,
    ncores = ncores,
    ...
  )

  accepted.methods <- noisyr::get_methods_calculate_noise_threshold()[8:10]
  if(!base::any(method.chosen %in% accepted.methods)){
    base::message("Method(s) chosen not suitable for the transcript approach, defaulting to Boxplot-IQR...")
    method.chosen <- "Boxplot-IQR"
  }else if(!base::all(method.chosen %in% accepted.methods)){
    base::message("Ignoring non-Boxplot methods for the transcript approach...")
    method.chosen <- method.chosen[method.chosen %in% accepted.methods]
  }

  if(base::length(similarity.threshold) > 1 | base::length(method.chosen) > 1){
    base::message("Selecting parameters that minimise the coefficient of variation...")
    stats.table <- noisyr::calculate_noise_threshold_method_statistics(
      expression = expression.summary,
      similarity.threshold.sequence = similarity.threshold,
      method.chosen.sequence = method.chosen,
      ...
    )
    row.min.coef.var <- base::which.min(stats.table$noise.threshold.coefficient.of.variation)
    similarity.threshold <- stats.table$similarity.threshold[row.min.coef.var]
    method.chosen <- paste(stats.table$approach[row.min.coef.var],
                           stats.table$method[row.min.coef.var],
                           sep="-")
  }

  noise.thresholds <- noisyr::calculate_noise_threshold(
    expression = expression.summary,
    similarity.threshold = similarity.threshold,
    method.chosen = method.chosen,
    ...
  )

  noisyr::remove_noise_from_bams(
    bams = bams,
    genes = genes,
    noise.thresholds = noise.thresholds,
    expression = expression.summary,
    ...
  )

  base::message(">>> Done! <<<")
  return(invisible(NULL))
}
