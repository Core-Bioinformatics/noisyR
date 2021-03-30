#' Function to tabulate statistics for different methods of calculating the noise threshold
#' @description This function is used to tabulate and compare different combinations of similarity
#' threshold and method to calculate the noise threshold for a given expression matrix.
#' @param expression either an expression summary (as calculated by calculate_expression_similarity_*()),
#' which should be a list with 3 slots: expression.matrix, expression.levels, expression.levels.similarity;
#' alternatively, just an expression matrix; only density based methods are available for the latter case
#' @param similarity.threshold.sequence similarity (correlation or inverse distance) threshold(s) to be used
#' to find corresponding noise threshold; can be a single value or a numeric vector;
#'  the default, 0.25 is usually suitable for the Pearson correlation (the default similarity measure)
#' @param method.chosen.sequence methods to use to calculate the noise thresholds,
#' must be a subset of get_methods_calculate_noise_threshold(); defaults to all
#' @param dump.stats name of csv to export different thresholds calculated (optional)
#' @param ... other arguments (for the boxplot methods) passed to calculate_noise_threshold_base()
#' @return A tibble containing information on noise thresholds calculated using the input
#' similarity thresholds and methods (optionally written in a csv file).
#' The columns list the chosen method and similarity threshold, the minimum, mean,
#' coefficient of variation, and maximum of the noise thresholds, and all the noise thresholds
#' concatenated as a string.
#' @export
#' @examples
#' expression.summary <- calculate_expression_similarity_counts(
#'     expression.matrix = matrix(1:100, ncol=5),
#'     method = "correlation_pearson",
#'     n.elements.per.window = 3)
#' calculate_noise_threshold_method_statistics(expression.summary)
calculate_noise_threshold_method_statistics <- function(
  expression,
  similarity.threshold.sequence = 0.25,
  method.chosen.sequence = noisyr::get_methods_calculate_noise_threshold(),
  dump.stats=NULL,
  ...
)
{
  if(base::is.matrix(expression)){
    expression.matrix <- expression
    expression.levels <- NULL
    expression.levels.similarity <- NULL
  }else if(base::is.list(expression) &
           identical(names(expression),
                     c("expression.matrix",
                       "expression.levels",
                       "expression.levels.similarity"))){
    expression.matrix <- expression$expression.matrix
    expression.levels <- expression$expression.levels
    expression.levels.similarity <- expression$expression.levels.similarity
  }else{
    stop("Please provide an expression.matrix or an expression.summary list")
  }

  if(!base::any(method.chosen.sequence %in% noisyr::get_methods_calculate_noise_threshold())){
    stop("Please provide a valid method from get_methods_calculate_noise_threshold()")
  }else if(!base::all(method.chosen.sequence %in% noisyr::get_methods_calculate_noise_threshold())){
    method.chosen.sequence <-
      method.chosen.sequence[method.chosen.sequence %in% noisyr::get_methods_calculate_noise_threshold()]
    base::warning("Dropped invalid method(s)")
  }

  if(base::is.null(expression.levels)){
    if(!base::any(method.chosen.sequence %in% noisyr::get_methods_calculate_noise_threshold()[1:3])){
      stop("Only density based methods are available for a simple matrix input")
    }else if(!base::all(method.chosen.sequence %in% noisyr::get_methods_calculate_noise_threshold()[1:3])){
      method.chosen.sequence <-
        method.chosen.sequence[method.chosen.sequence %in% noisyr::get_methods_calculate_noise_threshold()[1:3]]
      base::warning("Only density based methods are available for a simple matrix input - dropped invalid method(s)")
    }
  }

  stats.table <- tibble::tibble()

  for(similarity.threshold in similarity.threshold.sequence){
    for(method.chosen in method.chosen.sequence){
      noise.thresholds <- base::suppressMessages(
        noisyr::calculate_noise_threshold_base(
          expression = expression,
          similarity.threshold = similarity.threshold,
          method.chosen = method.chosen,
          ...
        )
      )
      stats.row <- tibble::tibble(approach=strsplit(method.chosen, "-")[[1]][1],
                                  method=strsplit(method.chosen, "-")[[1]][2],
                                  similarity.threshold=similarity.threshold,
                                  noise.threshold.min=base::min(noise.thresholds),
                                  noise.threshold.mean=base::mean(noise.thresholds),
                                  noise.threshold.coefficient.of.variation=
                                    stats::sd(noise.thresholds) / base::mean(noise.thresholds),
                                  noise.threshold.max=base::max(noise.thresholds),
                                  noise.thresholds.all=base::paste(noise.thresholds,
                                                                   collapse=","))
      stats.table <- base::rbind(stats.table, stats.row)
    }
  }

  if(!base::is.null(dump.stats)){
    utils::write.table(stats.table, file=dump.stats, sep=",", append=FALSE,
                       row.names=FALSE, col.names = TRUE)
  }

  return(stats.table)

}
