#' Calcualate the expression levels and expression levels similarity matrices using the count matrix
#' @description This function generates an average similarity (correlation/inverse distance) coefficient
#' for every sliding window, for each sample in the expression matrix.
#' That is done by comparing the distribution of genes in each window across samples.
#' @param expression.matrix the expression matrix, can be normalized or not
#' @param similarity.measure one of the correlation or distance metrics to be used,
#' defaults to pearson correlation; list of all methods in
#' \code{\link{get_methods_correlation_distance}}
#' @param n.elements.per.window number of elements to have in a window,
#' default 10\% of the number of rows
#' @param n.step step size to slide across, default 1\% of n.elements.per.window
#' @param n.step.fraction an alternative way to specify the step size, as a fraction of
#' the window length; default is 5\%
#' @param n.cores number of cores used to process samples in parallel (forked
#' processes via \code{parallel::mclapply}, so only used on unix-alikes); default 1
#' @param ... arguments passed on to other methods
#' @return A list with three elements: the first element is the expression matrix,
#'         as supplied; the other two are the expression levels matrix and
#'         expression levels similarity matrix;
#'         they have the same # of columns as the expression matrix,
#'         and n.elements.per.window * n.step rows.
#' @seealso \code{\link{calculate_expression_similarity_transcript}}
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
  n.cores=1,
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

  ## each function compares the window of the sorted sample (x)
  ## against the same rows of every other sample (columns of Y) in one call,
  ## instead of one call per pair of samples
  if(use.corr.dist=="correlation"){
    fun_corr_dist = function(x, Y){
      base::as.vector(stats::cor(x, Y, method=base.method))
    }
  }else if(use.corr.dist == "distance"){
    fun_corr_dist = function(x, Y){
      out <- base::rep(NA_real_, base::ncol(Y))
      sum.x <- base::sum(x)
      sums.Y <- base::colSums(Y)
      valid <- !(base::is.na(sums.Y) | sums.Y == 0)
      if(base::is.na(sum.x) || sum.x == 0 || !base::any(valid)) return(out)
      out[valid] <- base::suppressMessages(philentropy::dist_one_many(
        x / sum.x, base::t(Y[, valid, drop=FALSE]) / sums.Y[valid],
        method=base.method, p=3, testNA=FALSE))
      out
    }
  }

  nsequence = base::seq(1, (base::nrow(expression.matrix)-n.elements.per.window+1), n.step)
  window.offsets <- base::seq_len(n.elements.per.window) - 1
  n.samples <- base::ncol(expression.matrix)
  expression.matrix.unnamed <- base::unname(expression.matrix)

  ## for one sample, returns the average expression level of each sliding window
  ## and the average similarity of the window to the other samples
  process_sample <- function(j){
    base::message("  Working with sample ",j)
    expression.matrix.sorted = expression.matrix.unnamed[
      base::order(expression.matrix.unnamed[,j]), , drop=FALSE]
    levels <- base::numeric(base::length(nsequence))
    similarity <- base::numeric(base::length(nsequence))
    for(idx in base::seq_along(nsequence)){
      window <- expression.matrix.sorted[nsequence[idx] + window.offsets, , drop=FALSE]
      similarity.vector <- base::suppressWarnings(
        fun_corr_dist(window[, j], window[, -j, drop=FALSE]))
      similarity[idx] <- base::mean(similarity.vector[similarity.vector != 0])
      levels[idx] <- base::mean(window[, j])
    }
    base::list(levels=levels, similarity=similarity)
  }

  if(n.cores > 1 && .Platform$OS.type == "unix"){
    by.sample <- parallel::mclapply(base::seq_len(n.samples), process_sample, mc.cores=n.cores)
  }else{
    by.sample <- base::lapply(base::seq_len(n.samples), process_sample)
  }
  ## the average expression level per window is maintained in expression.levels variable
  ## the expression.levels.similarity matrix contains the similarity calculated for each sample,
  ##     at each expression level.
  expression.levels <- base::matrix(
    base::unlist(base::lapply(by.sample, `[[`, "levels")), ncol=n.samples)
  expression.levels.similarity <- base::matrix(
    base::unlist(base::lapply(by.sample, `[[`, "similarity")), ncol=n.samples)

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
