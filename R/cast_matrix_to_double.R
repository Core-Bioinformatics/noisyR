#' Cast a matrix of any type to double
#' @description Transforms values in the expression matrix to double,
#' to make it compatible with the rest of the functions.
#' @param expression.matrix The expression matrix (usually read from a file)
#' @return The expression matrix transformed to double, preserving row and column names
#' @export
#' @examples
#' cast_matrix_to_double(matrix(
#'     c(1, "2", 3.0, 4),
#'     ncol=2,
#'     dimnames=list(paste0("X", 1:2),
#'                   paste0("Y", 1:2))))
cast_matrix_to_double <- function(expression.matrix){
  expr.mat.int = base::matrix(base::rep(0,
                                        base::ncol(expression.matrix)*
                                          base::nrow(expression.matrix)),
                              ncol=base::ncol(expression.matrix))
  for(i in base::seq_len(base::nrow(expression.matrix))){
    for(j in base::seq_len(base::ncol(expression.matrix))){
      expr.mat.int[i,j] = base::as.numeric(base::as.character(expression.matrix[i,j]))
    }
  }

  base::colnames(expr.mat.int) = base::colnames(expression.matrix)
  base::rownames(expr.mat.int) = base::rownames(expression.matrix)
  return(expr.mat.int)
}
