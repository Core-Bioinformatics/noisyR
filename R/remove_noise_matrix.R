#' Function to remove the noisy reads from the expression matrix
#' @description This function is used to remove the noisy reads from the expression matrix.
#' It uses as input a vector of abundance thresholds;
#' all entries below the noise threshold are replaced with the noise threshold.
#' @param expression.matrix the expression matrix
#' @param abn.thresh a vector of abundance thresholds; must be the same length
#' as the number of columns of the expression matrix
#' @param add.thresh whether to add the noise threshold to all values in the expression matrix
#' (default), or set entries below the threshold to the threshold
#' @param average.thresh if TRUE (default), uses tthe average of the vector of thresholds across all samples;
#' if FALSE, uses the thresholds as supplied
#' @param remove.noisy.features logical, whether rows of the expression matrix that are
#' fully under the noise threshold should be removed (default TRUE)
#' @param export.csv option to write the matrix into a csv after the noise removal;
#' should be NULL or the name of the output file.
#' @return Returns a matrix of the same dims as the expression matrix,  with the noise removed.
#' This matrix has no entries remaining below the noise threshold.
#' @export
#' @examples remove_noise_matrix(
#'     expression.matrix = matrix(1:100, ncol=5),
#'     abn.thresh=c(5,30,45,62,83))
remove_noise_matrix <- function(expression.matrix, abn.thresh,
                                add.thresh=TRUE, average.thresh=TRUE,
                                remove.noisy.features=TRUE, export.csv=NULL){
  expression.matrix.noNoise <- expression.matrix
  thresh.mat <- base::matrix(base::rep(abn.thresh, base::nrow(expression.matrix)),
                             ncol=base::ncol(expression.matrix), byrow = TRUE)
  if(remove.noisy.features){
    above.noise.threshold <- base::as.vector(
      base::rowSums(expression.matrix >= thresh.mat) > 0)
    expression.matrix.noNoise <- expression.matrix.noNoise[above.noise.threshold,]
    thresh.mat <- thresh.mat[above.noise.threshold,]
  }
  if(average.thresh){
    abn.thresh.mean <- base::mean(abn.thresh)
    thresh.mat[] <- abn.thresh.mean
  }
  if(!add.thresh){
    expression.matrix.noNoise <- base::pmax(expression.matrix.noNoise, thresh.mat)
  }else{
    expression.matrix.noNoise <- expression.matrix.noNoise + thresh.mat

  }
  if(!base::is.null(export.csv)){
    utils::write.csv(expression.matrix.noNoise,
                     file=export.csv)
  }
  return(expression.matrix.noNoise)
}
