#' Function to remove the noisy reads from the expression matrix
#' @description This function is used to remove the noisy reads from the expression matrix.
#' It uses as input a vector of abundance thresholds;
#' all entries below the noise threshold are replaced with the noise threshold.
#' @param expression.matrix the expression matrix
#' @param abn.thr a vector of abundance thresholds; must be the same length
#' as the number of columns of the expression matrix
#' @param add.thr whether to add the noise threshold to all values in the expression matrix
#' (default), or set entries below the threshold to the threshold
#' @param average.thr if TRUE (default), uses tthe average of the vector of thresholds across all samples;
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
#'     abn.thr=c(5,30,45,62,83))
remove_noise_matrix <- function(expression.matrix, abn.thr,
                                add.thr=TRUE, average.thr=TRUE,
                                remove.noisy.features=TRUE, export.csv=NULL){
    if(base::length(abn.thr) == 1)
    abn.thr <- base::rep(abn.thr, base::ncol(expression.matrix))
  if(base::length(abn.thr) != base::ncol(expression.matrix)){
    base::message("abn.thr needs to be length 1 or ncol(expression.matrix)")
    break
  }
  expression.matrix.noNoise <- expression.matrix
  thr.mat <- base::matrix(base::rep(abn.thr, base::nrow(expression.matrix)),
                             ncol=base::ncol(expression.matrix), byrow = TRUE)
  if(remove.noisy.features){
    above.noise.threshold <- base::as.vector(
      base::rowSums(expression.matrix >= thr.mat) > 0)
    expression.matrix.noNoise <- expression.matrix.noNoise[above.noise.threshold,]
    thr.mat <- thr.mat[above.noise.threshold,]
  }
  if(average.thr){
    abn.thr.mean <- base::mean(abn.thr)
    thr.mat[] <- abn.thr.mean
  }
  if(!add.thr){
    expression.matrix.noNoise <- base::pmax(expression.matrix.noNoise, thr.mat)
  }else{
    expression.matrix.noNoise <- expression.matrix.noNoise + thr.mat

  }
  if(!base::is.null(export.csv)){
    utils::write.csv(expression.matrix.noNoise,
                     file=export.csv)
  }
  return(expression.matrix.noNoise)
}
