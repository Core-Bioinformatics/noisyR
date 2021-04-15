#' Optimise the elements per window for the count matrix approach
#' @description This function optimises the number of elements per window
#' that is used in \code{\link{calculate_expression_similarity_counts}}, by requiring
#' the distribution of correlations/distances to stabilise to a
#' uniform distribution. The Jensen-Shannon divergence is used to assess
#' the stability.
#' @param expression.matrix expression matrix, can be normalized or not
#' @param similarity.measure one of the correlation or distance metrics to be used,
#' defaults to pearson correlation; list of all methods in
#' \code{\link{get_methods_correlation_distance}}
#' @param window.length.min,window.length.max,window.length.by definition of the parameter search space;
#' default is between 1\% and 33\% of the number of rows in the expression matrix,
#' incremented by 1\%
#' @param n.step.fraction step size to slide across, as a fraction of
#' the window length; default is 5\%
#' @param iteration.number number of iterations for the subsampling and calculation of JSE;
#' subsampling is needed because shorter windows have fewer points; default is 100
#' @param minimum.similar.windows number of windows that a window needs to be similar to (including itself)
#' in order to be accepted as optimal; default is 3, but can be reduced to 2 if no optimum is found
#' @param save.plot name of the pdf in which to print the output plot
#' showing the distribution of JSE by window; output to the console by default
#' @return A single value of the optimal number of elements per window.
#' If no optimal value was found, this function returns NULL.
#' @export
#' @examples
#' optimise_window_length(
#'   matrix(1:100+runif(100), ncol=5, byrow=TRUE),
#'   window.length.min=3, window.length.max=5, iteration.number=5
#' )
optimise_window_length = function(
  expression.matrix,
  similarity.measure = "correlation_pearson",
  window.length.min=NULL,
  window.length.max=NULL,
  window.length.by=NULL,
  n.step.fraction=0.05,
  iteration.number=50,
  minimum.similar.windows=3,
  save.plot=NULL
)
{
  nrows <- base::nrow(expression.matrix)

  if(base::is.null(window.length.min)) window.length.min = base::max(base::floor(nrows * 0.01), 1)

  if(base::is.null(window.length.max)) window.length.max = base::max(base::floor(nrows * 0.33), 1)

  if(base::is.null(window.length.by)) window.length.by = base::max(base::floor(nrows * 0.01), 1)

  window.lengths <- base::seq(from=window.length.min, to=window.length.max, by=window.length.by)

  base::message("Window length optimisation")
  base::message("    number of windows: ", length(window.lengths))
  base::message("    minimum window: ", window.length.min)
  base::message("    maximum window: ", window.length.max)
  base::message("    window step: ", window.length.by)
  base::message("    # of iterations: ", iteration.number)
  base::message("    minimum similar windows: ", minimum.similar.windows)
  base::message("Calculating expression summary and JSE for each window...")

  tbl.jse.all <- tibble::tibble()
  for(cnt in base::rev(base::seq_along(window.lengths))){

    window.length <- window.lengths[cnt]
    tbl.jse.temp <- tibble::tibble(window.length=base::rep(window.length, iteration.number),
                                   iteration=1:iteration.number,
                                   JSE=rep(NA_real_, iteration.number))
    base::suppressMessages(
      expression.summary <- noisyr::calculate_expression_similarity_counts(
        expression.matrix=expression.matrix,
        n.elements.per.window=window.length,
        n.step.fraction=n.step.fraction,
        similarity.measure=similarity.measure)
    )
    row.means <- base::rowMeans(expression.summary$expression.levels.similarity, na.rm=TRUE)
    x=base::abs(row.means[!base::is.na(row.means)])
    x.length <- base::length(x)
    if(cnt == base::length(window.lengths)) subsize <- x.length
    yy=base::seq(from=base::min(abs(row.means), na.rm=TRUE),
                 to=base::max(row.means, na.rm=TRUE),
                 length.out=subsize)

    for(iter in base::seq_len(iteration.number)){
      xx <- x[base::sort(base::sample(1:x.length, subsize))]
      tbl.jse.temp[iter, "JSE"] <- philentropy::jensen_shannon(xx, yy, testNA=FALSE, unit="log")
    }

    tbl.jse.all <- base::rbind(tbl.jse.temp, tbl.jse.all)
  }

  JSE=NULL
  p <- ggplot2::ggplot(tbl.jse.all) +
    ggplot2::theme_minimal() +
    ggplot2::geom_boxplot(ggplot2::aes(x=base::as.factor(window.length), y=JSE)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
    ggplot2::xlab("Window length")
  if(!base::is.null(save.plot)) grDevices::pdf(save.plot)
  print(p)
  if(!base::is.null(save.plot)) grDevices::dev.off()

  jse.by.window.list <- base::list()
  for(i in base::seq_along(window.lengths)){
    jse.by.window.list[[i]] = dplyr::filter(tbl.jse.all, window.length==window.lengths[i])$JSE
  }

  base::message("Performing t-tests...")
  ttests.matrix = base::matrix(0, nrow=length(window.lengths), ncol=length(window.lengths))
  for(i in base::seq_len(base::length(window.lengths))){
    ttest.positive.by.window <- base::vector()
    for(j in base::seq_len(base::length(window.lengths))){
      if(i!=j){
        ttest <- stats::t.test(jse.by.window.list[[i]], jse.by.window.list[[j]])
        ttest.positive.by.window <- base::c(ttest.positive.by.window, ttest$p.value >= 0.05)
      }else{
        ttest.positive.by.window <- base::c(ttest.positive.by.window, 1)
      }
    }
    ttests.matrix[i,] <- ttest.positive.by.window
  }

  base::message("The number of similar windows found for each window (including itself) were:")
  base::message("    ", base::paste(base::rowSums(ttests.matrix), collapse=" "))
  enough.similar.windows <- base::rowSums(ttests.matrix) >= minimum.similar.windows
  optimal.window.length.id <- match(TRUE, enough.similar.windows)

  if(!base::is.na(optimal.window.length.id)){
    optimal.window.length <- window.lengths[optimal.window.length.id]
    base::message("Optimal window found!")
  }else{
    optimal.window.length <- NULL
    base::message("Optimal window not found, consider less stringent parameters, e.g.")
    base::message("    > reducing minimum.smaller.windows")
    base::message("    > reducing iteration.number")
    base::message("    > reducing window.length.by")
  }
  return(optimal.window.length)
}
