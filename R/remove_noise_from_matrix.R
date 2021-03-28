#' Function to remove the noisy reads from the expression matrix
#' @description This function is used to remove the noisy reads from the expression matrix.
#' It uses as input a vector of abundance thresholds;
#' all entries below the noise threshold are replaced with the noise threshold.
#' @param expression.matrix the expression matrix
#' @param noise.thresholds a vector of expression thresholds by sample;
#' must be the same length as the number of columns of the expression matrix,
#' or a singular value to be used as a fixed noise threshold
#' @param add.threshold whether to add the noise threshold to all values in the expression matrix
#' (default), or set entries below the threshold to the threshold
#' @param average.threshold if TRUE (default), uses tthe average of the vector of thresholds across all samples;
#' if FALSE, uses the thresholds as supplied
#' @param remove.noisy.features logical, whether rows of the expression matrix that are
#' fully under the noise threshold should be removed (default TRUE)
#' @param export.csv option to write the matrix into a csv after the noise removal;
#' should be NULL or the name of the output file.
#' @return Returns a matrix of the same dims as the expression matrix, with the noise removed.
#' This matrix has no entries remaining below the noise threshold.
#' @export
#' @examples
#' expression.matrix <- matrix(1:100, ncol=5)
#' noise.thresholds <- c(5,30,45,62,83)
#' expression.matrix.denoised <- remove_noise_from_matrix(
#'     expression.matrix = expression.matrix,
#'     noise.thresholds = noise.thresholds
#' )
#'
remove_noise_from_matrix <- function(expression.matrix, noise.thresholds,
                                     add.threshold=TRUE, average.threshold=TRUE,
                                     remove.noisy.features=TRUE, export.csv=NULL){
  if(base::length(noise.thresholds) == 1){
    base::message("noise.thresholds only has 1 value, using a fixed threshold...")
    noise.thresholds <- base::rep(noise.thresholds, base::ncol(expression.matrix))
  }else if(base::length(noise.thresholds) != base::ncol(expression.matrix)){
    base::stop("noise.thresholds needs to be length 1 or ncol(expression.matrix)")
  }
  expression.matrix.denoised <- expression.matrix
  threshold.matrix <- base::matrix(base::rep(noise.thresholds, base::nrow(expression.matrix)),
                                   ncol=base::ncol(expression.matrix), byrow = TRUE)
  if(remove.noisy.features){
    above.noise.threshold <- base::as.vector(
      base::rowSums(expression.matrix >= threshold.matrix) > 0)
    expression.matrix.denoised <- expression.matrix.denoised[above.noise.threshold,]
    threshold.matrix <- threshold.matrix[above.noise.threshold,]
  }
  if(average.threshold){
    noise.thresholds.mean <- base::mean(noise.thresholds)
    threshold.matrix[] <- noise.thresholds.mean
  }

  threshold.matrix <- base::round(threshold.matrix)

  if(!add.threshold){
    expression.matrix.denoised <- base::pmax(expression.matrix.denoised, threshold.matrix)
  }else{
    expression.matrix.denoised <- expression.matrix.denoised + threshold.matrix
  }
  if(!base::is.null(export.csv)){
    utils::write.csv(expression.matrix.denoised,
                     file=export.csv)
  }
  return(expression.matrix.denoised)
}
