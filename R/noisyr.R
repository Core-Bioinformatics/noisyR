#' Run the noisyR pipeline
#' @description Calls one of noisyr_counts or noisyr_transcript, with the specified parameters.
#' See the individual function documentation for more details and required arguments:
#' \code{\link{noisyr_counts}}, \code{\link{noisyr_transcript}}
#' @param approach.for.similarity.calculation which approach to use for the similarity calculation;
#' defaults to counts
#' @param ... arguments to be passed on to noisyr_counts or noisyr_transcript; see their documentation
#' for more details and required arguments
#' @return For the counts approach, the denoised expression matrix. For the transcript approach,
#' the numeric vector of noise thresholds per sample. For more details, see their respective documentation.
#' @export
#' @examples
#' noisyr(approach.for.similarity.calculation = "counts",
#'        expression.matrix = matrix(1:100, ncol = 5))
noisyr = function(approach.for.similarity.calculation = c("counts", "transcript"), ...){
  if(approach.for.similarity.calculation[1] == "counts"){
    noisyr::noisyr_counts(...)
  }else if(approach.for.similarity.calculation[1] == "transcript"){
    noisyr::noisyr_transcript(...)
  }else{
    base::stop("Invalid approach, please select counts or transcript.")
  }
}
