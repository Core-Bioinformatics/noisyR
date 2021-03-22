#' Calcualte the distance matrices using the count matrix
#' @description This function generates an average correlation/distance coefficient
#' for every sliding window, for each sample in the expression matrix.
#' That is done by comparing the distribution of genes in each window across samples.
#' @param expression.matrix expression matrix, can be normalized or not
#' @param method one of the correlation or distance metrics to be used,
#' defaults to pearson correlation; list of all methods in
#' get_methods_correlation_distance()
#' @param n.elements.per.window number of elements to have in a window,
#' default 10\% of the number of rows
#' @param nstep step size to slide across, default 1\% of n.elements.per.window
#' @param nstep.frac an alternative way to specify the step size, as a fraction of
#' the window length; default is 5\%
#' @return A list with three elements: the first element is the expression matrix,
#'         as supplied; the other two are the abundance and distance matrix;
#'         they have the same # of columns as the expression matrix,
#'         and n.elements.per.window*nstep rows.
#' @export
#' @examples
#' calculate_distance_matrices_counts(
#'     expression.matrix = matrix(1:100, ncol=5),
#'     method="correlation_pearson",
#'     n.elements.per.window=3)
calculate_distance_matrices_counts = function(expression.matrix,
                                              method = "correlation_pearson",
                                              n.elements.per.window=NULL,
                                              nstep=NULL,
                                              nstep.frac=0.05)
{
  base::message("The input matrix has ", base::nrow(expression.matrix), " rows and ",
            base::ncol(expression.matrix), " cols")
  base::message("    number of genes: ", base::nrow(expression.matrix))
  base::message("    number of samples: ", base::ncol(expression.matrix))

  if(base::is.null(n.elements.per.window))
  {
    base::message("Calculating the number of elements per window")
    n.elements.per.window = base::round(base::nrow(expression.matrix) / 10)
    nstep <- NULL
  }
  base::message("    the number of elements per window is ", n.elements.per.window)
  if(base::is.null(nstep))
  {
    nstep = base::max(base::floor(n.elements.per.window * nstep.frac), 1)
  }else{
    nstep.frac <- base::round(nstep / n.elements.per.window, 2)
  }
  base::message("    the step size is ", nstep)

  if(method %in% noisyr::get_methods_correlation_distance()){
    base::message("Method chosen: ", method)
  }else{
    base::message("Invalid method, reverting to Pearson correlation")
    method <- "correlation_pearson"
  }
  use.corr.dist <- base::strsplit(method, "_")[[1]][1]
  base.method <- base::sub(paste0(use.corr.dist,"_"), "", method)
  if(use.corr.dist=="correlation"){
    fun_corr_dist = function(cols.in){
      stats::cor(cols.in, method=base.method)[1,2]
    }
  }else if(use.corr.dist=="distance"){
    fun_corr_dist = function(cols.in){
      vecsums <- base::colSums(cols.in)
      if(base::any(is.na(vecsums) | vecsums==0)){
        return(NA)
      }else{
        return(suppressMessages(
          philentropy::distance(base::t(cols.in)/vecsums,
                                method=base.method, p=3, test.na=FALSE)))
      }
    }
  }

  nsequence = base::seq(1, (base::nrow(expression.matrix)-n.elements.per.window+1), nstep)
  ##the matrix contains the distance calculated for each sample, at each abundance level.
  ##the abundance level is maintained in abn.matrix variable
  dist.matrix = base::rep(0, base::ncol(expression.matrix)*base::length(nsequence))
  dist.matrix = base::matrix(dist.matrix, ncol=base::ncol(expression.matrix))
  abn.matrix = base::rep(0, base::ncol(expression.matrix)*base::length(nsequence))
  abn.matrix = base::matrix(abn.matrix, ncol=base::ncol(expression.matrix))

  for(j in base::seq_len(base::ncol(expression.matrix))){
    base::message("Working with sample ",j)
    sorted.matrix = expression.matrix[base::order(expression.matrix[,j]),]
    for(idx in base::seq_len(base::length(nsequence))){
      #focus on a sliding window, initialize the distance vector
      distance.vector = base::vector(mode="numeric", length=base::ncol(expression.matrix))
      for(k in base::seq_len(base::ncol(expression.matrix))){
        if(j != k){
          col.j = sorted.matrix[nsequence[idx]:(nsequence[idx]+n.elements.per.window-1), j];
          base::names(col.j)=base::rep("", base::length(col.j))
          col.k = sorted.matrix[nsequence[idx]:(nsequence[idx]+n.elements.per.window-1), k];
          base::names(col.k)=base::rep("", base::length(col.k))
          cols.in = base::cbind(col.j, col.k)
          distance.vector[k] = base::suppressWarnings(fun_corr_dist(cols.in))
        }#end if
      }#end for k
      dist.matrix[idx,j] = base::mean(distance.vector[distance.vector != 0])
      abn.matrix[idx,j] = base::mean(sorted.matrix[nsequence[idx]:
                                                     (nsequence[idx]+n.elements.per.window-1),j])
    }#end for on the sliding windows
  }

  if(noisyr::get_methods_correlation_distance(names=FALSE)[
    base::match(method, noisyr::get_methods_correlation_distance())] == "d"){
    base::message("Chosen method ", method, " is a dissimilarity measure, outputting inverse...")
    dist.matrix <- 1/dist.matrix
  }

  returnObject <- base::list("exp" = expression.matrix,
                             "abn" = abn.matrix,
                             "dist" = dist.matrix)
  return(returnObject)
}
