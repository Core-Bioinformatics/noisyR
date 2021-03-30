#' Calcualte the expression levels and expression levels similarity matrices using the count matrix
#' @description This function generates an average similarity (correlation/inverse distance) coefficient
#' for every sliding window, for each sample in the expression matrix.
#' That is done by comparing the distribution of genes in each window across samples.
#' @param expression.matrix the expression matrix, can be normalized or not
#' @param similarity.measure one of the correlation or distance metrics to be used,
#' defaults to pearson correlation; list of all methods in
#' get_methods_correlation_distance()
#' @param n.elements.per.window number of elements to have in a window,
#' default 10\% of the number of rows
#' @param n.step step size to slide across, default 1\% of n.elements.per.window
#' @param n.step.fraction an alternative way to specify the step size, as a fraction of
#' the window length; default is 5\%
#' @param ... arguments passed on to other methods
#' @return A list with three elements: the first element is the expression matrix,
#'         as supplied; the other two are the expression levels matrix and
#'         expression levels similarity matrix;
#'         they have the same # of columns as the expression matrix,
#'         and n.elements.per.window * n.step rows.
#' @export
#' @examples
#' calculate_expression_similarity_counts(
#'     expression.matrix = matrix(1:100, ncol = 5),
#'     similarity.measure = "correlation_pearson",
#'     n.elements.per.window = 3)
calculate_expression_similarity_counts = function(
  expression.matrix,
  similarity.measure = "correlation_pearson",
  n.elements.per.window=NULL,
  n.step=NULL,
  n.step.fraction=0.05,
  ...
){
  base::message("The input matrix has ", base::nrow(expression.matrix), " rows and ",
                base::ncol(expression.matrix), " cols")
  base::message("    number of genes: ", base::nrow(expression.matrix))
  base::message("    number of samples: ", base::ncol(expression.matrix))

  if(base::is.null(n.elements.per.window))
  {
    base::message("Calculating the number of elements per window")
    n.elements.per.window = base::round(base::nrow(expression.matrix) / 10)
    n.step <- NULL
  }
  base::message("    the number of elements per window is ", n.elements.per.window)
  if(base::is.null(n.step))
  {
    n.step = base::max(base::floor(n.elements.per.window * n.step.fraction), 1)
  }else{
    n.step.fraction <- base::round(n.step / n.elements.per.window, 2)
  }
  base::message("    the step size is ", n.step)

  if(similarity.measure %in% noisyr::get_methods_correlation_distance()){
    base::message("    the selected similarity metric is ", similarity.measure)
  }else{
    base::message("Invalid similarity metric, reverting to Pearson correlation")
    similarity.measure <- "correlation_pearson"
  }
  use.corr.dist <- base::strsplit(similarity.measure, "_")[[1]][1]
  base.method <- base::sub(paste0(use.corr.dist,"_"), "", similarity.measure)
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

  nsequence = base::seq(1, (base::nrow(expression.matrix)-n.elements.per.window+1), n.step)
  ## the average expression level per window is maintained in expression.levels variable
  ## the expression.levels.similarity matrix contains the similarity calculated for each sample,
  ##     at each expression level.
  expression.levels = base::rep(0, base::ncol(expression.matrix)*base::length(nsequence))
  expression.levels = base::matrix(expression.levels, ncol=base::ncol(expression.matrix))
  expression.levels.similarity = base::rep(0, base::ncol(expression.matrix)*base::length(nsequence))
  expression.levels.similarity = base::matrix(expression.levels.similarity, ncol=base::ncol(expression.matrix))

  for(j in base::seq_len(base::ncol(expression.matrix))){
    base::message("  Working with sample ",j)
    expression.matrix.sorted = expression.matrix[base::order(expression.matrix[,j]),]
    for(idx in base::seq_len(base::length(nsequence))){
      #focus on a sliding window, initialize the similarity vector
      similarity.vector = base::vector(mode="numeric", length=base::ncol(expression.matrix))
      for(k in base::seq_len(base::ncol(expression.matrix))){
        if(j != k){
          col.j = expression.matrix.sorted[nsequence[idx]:(nsequence[idx]+n.elements.per.window-1), j];
          base::names(col.j)=base::rep("", base::length(col.j))
          col.k = expression.matrix.sorted[nsequence[idx]:(nsequence[idx]+n.elements.per.window-1), k];
          base::names(col.k)=base::rep("", base::length(col.k))
          cols.in = base::cbind(col.j, col.k)
          similarity.vector[k] = base::suppressWarnings(fun_corr_dist(cols.in))
        }#end if
      }#end for k
      expression.levels.similarity[idx,j] = base::mean(similarity.vector[similarity.vector != 0])
      expression.levels[idx,j] = base::mean(expression.matrix.sorted[nsequence[idx]:
                                                                       (nsequence[idx]+n.elements.per.window-1),j])
    }#end for on the sliding windows
  }

  if(noisyr::get_methods_correlation_distance(names=FALSE)[
    base::match(similarity.measure, noisyr::get_methods_correlation_distance())] == "d"){
    base::message("Chosen similarity metric ", similarity.measure, " is a dissimilarity, outputting inverse...")
    expression.levels.similarity <- 1/expression.levels.similarity
  }

  expression.summary <- base::list("expression.matrix" = expression.matrix,
                                   "expression.levels" = expression.levels,
                                   "expression.levels.similarity" = expression.levels.similarity)
  return(expression.summary)
}
