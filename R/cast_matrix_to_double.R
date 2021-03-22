#' Cast a matrix of any type to double
#' @description Transforms values in the expression matrix to double,
#' to make it compatible with the rest of the functions.
#' @param expression.matrix The expression matrix (usually read from a file)
#' @return The expression matrix transformed to double, preserving row and column names.
#' Any values that are not coercible to numeric are replaced by 0.
#' @export
#' @examples
#' cast_matrix_to_double(matrix(
#'     c(1, "2", 3.0, 4),
#'     ncol=2,
#'     dimnames=list(paste0("X", 1:2),
#'                   paste0("Y", 1:2))))
cast_matrix_to_double <- function(expression.matrix){
  expr.mat.int = base::suppressWarnings(
    base::apply(expression.matrix, 2, base::as.numeric))
  mat.nonint <- base::is.na(expr.mat.int)
  if(base::any(mat.nonint)){
    base::warning("Count matrix contains non-numeric values: converted to 0.")
    expr.mat.int[mat.nonint] <- 0
  }
  base::rownames(expr.mat.int) = base::rownames(expression.matrix)
  return(expr.mat.int)
}
