#' Show the methods for calculating a noise threshold
#' @description This function outputs the methods available for the calculation of
#' the noise threshold. To be used as input in \code{\link{calculate_noise_threshold}}.
#' @return A character vector of options for the method.chosen arguement of
#' \code{\link{calculate_noise_threshold}}
#' @export
#' @examples get_methods_calculate_noise_threshold()
get_methods_calculate_noise_threshold = function(){
  approaches <- base::c(base::rep("Density_based", 3),
                        base::rep("Line_plot", 4),
                        base::rep("Boxplot", 3))
  methods <- base::c("No_normalisation",
                     "RPM_normalisation",
                     "Quantile_normalisation",
                     "No_smoothing",
                     "loess10_smoothing",
                     "loess25_smoothing",
                     "loess50_smoothing",
                     "Median",
                     "IQR",
                     "Quant5")
  vec <- base::paste(approaches, methods, sep="-")
  return(vec)
}
