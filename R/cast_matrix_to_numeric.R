#' Cast a matrix of any type to numeric
#' @description Transforms values in the expression matrix to numeric,
#' to make it compatible with the rest of the functions.
#' @param expression.matrix The expression matrix (usually read from a file)
#' @return The expression matrix transformed to numeric, preserving row and column names.
#' Any values that are not coercible to numeric are replaced by 0.
#' @export
#' @examples
#' cast_matrix_to_numeric(matrix(
#'     c(1, "2", 3.0, 4),
#'     ncol=2,
#'     dimnames=list(paste0("X", 1:2),
#'                   paste0("Y", 1:2))))
cast_matrix_to_numeric <- function(expression.matrix){
  expression.matrix.numeric = base::suppressWarnings(
    base::apply(expression.matrix, 2, base::as.numeric))
  matrix.na <- base::is.na(expression.matrix.numeric)
  if(base::any(matrix.na)){
    base::warning("Count matrix contains non-numeric values: converted to 0.")
    expression.matrix.numeric[matrix.na] <- 0
  }
  base::rownames(expression.matrix.numeric) = base::rownames(expression.matrix)
  return(expression.matrix.numeric)
}
