#' Function to remove the noisy reads from the expression matrix using a specified method
#' @description This function is used to remove the noisy reads from the expression matrix.
#' It can use thresholds precalculated by calculate_threshold_noise(), if available.
#' If a vector of thresholds needs to be input manually, use remove_noise_matrix() instead.
#' @param expression.matrix the expression matrix
#' @param method.chosen the method to be used for the noise removal;
#' must be one of get_methods_calculate_noise_threshold()
#' @param stats.df a tibble, as output by calculate_threshold_noise();
#' if supplied and the chosen method appears in it, the corresponding threshold is used directly.
#' @param abn.matrix,dist.matrix the input distance and abundance matrices as calculated by
#' calculate_distance_matrices(); only needed if stats.df is not supplied;
#' if either is not supplied, only a fixed threshold is
#' calculated based on the density
#' @param dist.thresh a distance threshold to be used if the noise thresholds are not
#' pre-calculated; the default 0.25 is suitable for correlation measures
#' @param binsize size of each bin in the boxplot methods; defaults to 0.1 (on a log-scale)
#' @param add.thresh whether to add the noise threshold to all values in the expression matrix
#' (default), or set entries below the threshold to the threshold
#' @param average.thresh if TRUE (default), uses tthe average of the vector of thresholds across all samples;
#' if FALSE, uses the thresholds as supplied
#' @param remove.noisy.features logical, whether rows of the expression matrix that are
#' fully under the noise threshold should be removed (default TRUE)
#' @param export.csv option to write the matrix into a csv after the noise removal;
#' should be NULL or the name of the output file.
#' @return Returns a matrix of the same dims as the expression matrix, with the noise removed
#' using the specified method.
#' All entries below the noise threshold are replaced with the noise threshold.
#' @export
#' @examples obj <- calculate_distance_matrices_counts(
#'     expression.matrix = matrix(1:100, ncol=5),
#'     method="correlation_pearson",
#'     n.elements.per.window=3)
#' remove_noise_method(
#'     expression.matrix=obj$exp,
#'     abn.matrix=obj$abn,
#'     dist.matrix=obj$dist)
remove_noise_method = function(expression.matrix,
                               method.chosen="Boxplot-IQR",
                               stats.df=NULL,
                               abn.matrix=NULL, dist.matrix=NULL,
                               dist.thresh=0.25, binsize=0.1,
                               add.thresh=TRUE, average.thresh=TRUE,
                               remove.noisy.features=TRUE, export.csv=NULL){
  approach <- base::strsplit(method.chosen, split="-")[[1]][1]
  method <- base::strsplit(method.chosen, split="-")[[1]][2]
  abn.thresh <- NULL
  if(!base::is.null(stats.df)){
    methods.preran <- base::paste(stats.df$approach, stats.df$method, sep="-")
    if(method.chosen %in% methods.preran){
      abn.thresh <- base::as.numeric(base::strsplit(
        stats.df$abn.thresh.all[match(method.chosen, methods.preran)],
        split=",")[[1]])

    }
  }
  if(base::is.null(abn.thresh)){
    if(approach!="Density_based_fixed_threshold" &
       (base::is.null(dist.matrix) | base::is.null(abn.matrix))){
      stop(base::message("Method", method.chosen, "requires a distance and abundance matrix"))
    }
    abn.thresh <- noisyr::calculate_threshold_noise(
      expression.matrix=expression.matrix,
      dist.matrix=dist.matrix,
      abn.matrix=abn.matrix,
      dist.thresh=dist.thresh,
      binsize=binsize,
      dump.stats=NULL,
      method.chosen=method.chosen)
  }
  expression.matrix.noNoise <- noisyr::remove_noise_matrix(
    expression.matrix=expression.matrix,
    abn.thresh=abn.thresh,
    add.thresh = add.thresh,
    average.thresh = average.thresh,
    remove.noisy.features=remove.noisy.features,
    export.csv=export.csv)
  return(expression.matrix.noNoise)
}
