#' Run the noisyR pipeline
#' @description Calls one of noisyr_counts or noisyr_transcript, with the specified parameters.
#' @param approach_for_similarity_calculation which approach to use for the similarity calculation;
#' defaults to counts
#' @param ... arguments to be passed on to noisyr_counts or noisyr_transcript; see their documentation
#' for more details and required arguments
#' @return For the counts approach, the denoised expression matrix. For the transcript approach,
#' the numeric vector of noise thresholds per sample. For more details, see their respective documentation.
#' @seealso [noisyr_counts()], [noisyr_transcript()]
#' @export
#' @examples
#' noisyr(approach_for_similarity_calculation = "counts",
#'        expression.matrix = matrix(1:100, ncol = 5))
noisyr = function(approach_for_similarity_calculation = c("counts", "transcript"), ...){
  if(approach_for_similarity_calculation[1] == "counts"){
    noisyr::noisyr_counts(...)
  }else if(approach_for_similarity_calculation[1] == "transcript"){
    noisyr::noisyr_transcript(...)
  }else{
    base::stop("Invalid approach, please select counts or transcript.")
  }
}
