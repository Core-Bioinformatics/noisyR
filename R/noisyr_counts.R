#' Run the noisyR pipeline for the count matrix approach
#' @description Calls the functions to run each of the three steps of the pipeline
#' (similarity calculation, noise quantification, noise removal), with the specified parameters.
#' See the individual function documentation for more details and required arguments.
#' Required steps: \code{\link{calculate_expression_similarity_counts}},
#' \code{\link{calculate_noise_threshold}}. \code{\link{remove_noise_from_matrix}}.
#' Optional steps: \code{\link{optimise_window_length}},
#' \code{\link{calculate_noise_threshold_method_statistics}}
#' @param expression.matrix the expression matrix used as input for the similarity calculation;
#' this argument is required
#' @param n.elements.per.window number of elements to have in a window passed to
#' calculate_expression_similarity_counts(); default 10\% of the number of rows
#' @param optimise.window.length.logical whether to call optimise_window_length to try and
#' optimise the value of n.elements.per.window
#' @param similarity.threshold,method.chosen parameters passed on to
#' \code{\link{calculate_noise_threshold}}; they can be single values or vectors;
#' if they are vectors optimal values are computed by calling
#' \code{\link{calculate_noise_threshold_method_statistics}} and
#' minimising the coefficient of variation across samples; all possible values for
#' method.chosen can be viewed by \code{\link{get_methods_calculate_noise_threshold}}
#' @param ... arguments to be passed on to individual pipeline steps
#' @return The denoised expression matrix.
#' @seealso \code{\link{noisyr}}, \code{\link{noisyr_transcript}}
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

  expression.matrix.denoised <- noisyr::remove_noise_from_matrix(
    expression.matrix = expression.matrix,
    noise.thresholds = noise.thresholds,
    ...
  )
  base::message(">>> Done! <<<")
  return(expression.matrix.denoised)
}
