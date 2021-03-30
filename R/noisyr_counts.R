#' Run the noisyR pipeline for the count matrix approach
#' @description Calls the functions to run each of the three steps of the pipeline
#' (similarity calculation, noise quantification, noise removal), with the specified parameters.
#' @param expression.matrix the expression matrix used as input for the similarity calculation;
#' this argument is required
#' @param n.elements.per.window number of elements to have in a window passed to
#' calculate_expression_similarity_counts(); default 10\% of the number of rows
#' @param optimise.window.length.logical whether to call optimise_window_length to try and
#' optimise the value of n.elements.per.window
#' @param similarity.threshold,method.chosen parameters passed on to calculate_noise_threshold_base()
#' @param minimise.coefficient.of.variation whether to compute optimal values for similarity.threshold
#' and method.chosen by minimising the coefficient of variation across samples
#' @param ... arguments to be passed on to individual pipeline steps; see their documentation
#' for more details and required arguments
#' @return For the counts approach, the denoised expression matrix. For the transcript approach,
#' the numeric vector of noise thresholds per sample. For more details, see their respective documentation.
#' @seealso [noisyr()], [noisyr_transcript()]
#' @export
#' @examples
#' noisyr_counts(
#'     expression.matrix = matrix(1:100, ncol = 5),
#'     similarity.measure = "correlation_pearson",
#'     n.elements.per.window = 3)
noisyr_counts = function(
  expression.matrix,
  n.elements.per.window = NULL,
  optimise.window.length.logical = FALSE,
  similarity.threshold = 0.25,
  method.chosen = "Boxplot-IQR",
  minimise.coefficient.of.variation = FALSE,
  ...
){
  base::message(">>> noisyR counts approach pipeline <<<")
  expression.matrix <- noisyr::cast_matrix_to_numeric(expression.matrix)

  if(optimise.window.length.logical){
    n.elements.per.window <- noisyr::optimise_window_length(
      expression.matrix = expression.matrix,
      ...
    )
  }

  expression.summary <- noisyr::calculate_expression_similarity_counts(
    expression.matrix = expression.matrix,
    n.elements.per.window = n.elements.per.window,
    ...
  )

  if(minimise.coefficient.of.variation){
    base::message("Selecting parameters that minimise the coefficient of variation...")
    stats.table <- noisyr::calculate_noise_threshold_method_statistics(
      expression = expression.summary,
      ...
    )
    row.min.coef.var <- base::which.min(stats.table$noise.threshold.coefficient.of.variation)
    similarity.threshold <- stats.table$similarity.threshold[row.min.coef.var]
    method.chosen <- paste(stats.table$approach[row.min.coef.var],
                           stats.table$method[row.min.coef.var],
                           sep="-")
  }

  noise.thresholds <- noisyr::calculate_noise_threshold_base(
    expression = expression.summary,
    similarity.threshold = similarity.threshold,
    method.chosen = method.chosen,
    ...
  )

  expression.matrix.denoised <- noisyr::remove_noise_from_matrix(
    expression.matrix = expression.matrix,
    noise.thresholds = noise.thresholds,
    ...
  )
  base::message(">>> Done! <<<")
  return(expression.matrix.denoised)
}
