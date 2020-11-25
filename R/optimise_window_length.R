#' Optimise the elements per window for the count matrix approach
#' @description This function optimises the number of elements per window
#' that is used in calculate_distance_matrices_counts(), by requiring
#' the distribution of correlations/distances to stabilise to a
#' uniform distribution. The Jensen-Shannon divergence is used to assess
#' the stability.
#' @param expression.matrix expression matrix, can be normalized or not
#' @param method one of the correlation or distance metrics to be used,
#' defaults to pearson correlation; list of all methods in
#' get_methods_correlation_distance()
#' @param winlen.min,winlen.max,winlen.by definition of the parameter search space
#' @param nstep.frac an alternative way to specify the step size, as a fraction of
#' the window length; default is 5\%
#' @param iternum number of iterations for the subsampling and calculation of JSE;
#' subsampling is needed because shorter windows have fewer points
#' @param save.plot name of the pdf in which to print the output plot
#' showing the distribution of JSE by window
#' @return A single value of the optimal number of elements per window
#' @export
#' @examples optimise_window_length(
#' matrix(1:100+runif(100), ncol=5, byrow=TRUE),
#'   winlen.min=3, winlen.max=5, iternum=5)
optimise_window_length = function(
  expression.matrix,
  method = "correlation_pearson",
  winlen.min=NULL,
  winlen.max=NULL,
  winlen.by=NULL,
  nstep.frac=0.05,
  iternum=1000,
  save.plot=NULL
)
{
  nrows <- base::nrow(expression.matrix)

  if(base::is.null(winlen.min)) winlen.min = base::max(base::floor(nrows * 0.01), 1)

  if(base::is.null(winlen.max)) winlen.max = base::max(base::floor(nrows * 0.33), 1)

  if(base::is.null(winlen.by)) winlen.by = base::max(base::floor(nrows * 0.01), 1)

  winlens <- base::seq(from=winlen.min, to=winlen.max, by=winlen.by)

  tbl <- tibble::tibble()
  for(cnt in base::rev(base::seq_along(winlens))){

    winlen <- winlens[cnt]
    diff <- tibble::tibble(winlen=base::rep(winlen, iternum),
                           iter=1:iternum,
                           JSE=rep(NA_real_, iternum))
    base::suppressMessages(
      obj <- noisyr::calculate_distance_matrices_counts(expression.matrix=expression.matrix,
                                                        n.elements.per.window=winlen,
                                                        nstep.frac=nstep.frac,
                                                        method=method)
    )
    rmeans <- base::rowMeans(obj$dist, na.rm=TRUE)
    x=base::abs(rmeans[!base::is.na(rmeans)])
    xlen <- base::length(x)
    if(cnt == length(winlens)) subsize <- xlen
    yy=seq(from=base::min(abs(rmeans), na.rm=TRUE),
           to=base::max(rmeans, na.rm=TRUE),
           length.out=subsize)

    for(iter in base::seq_len(iternum)){
      xx <- x[base::sort(base::sample(1:xlen, subsize))]
      diff[iter, "JSE"] <- philentropy::jensen_shannon(xx, yy, testNA=FALSE, unit="log")
    }

    tbl <- base::rbind(diff, tbl)
  }

  JSE=NULL
  p <- ggplot2::ggplot(tbl) +
    ggplot2::geom_boxplot(ggplot2::aes(x=base::as.factor(winlen), y=JSE)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
    ggplot2::xlab("Window length")
  if(!base::is.null(save.plot)) grDevices::pdf(save.plot)
  print(p)
  if(!base::is.null(save.plot)) grDevices::dev.off()

  dists <- base::list()
  for(i in base::seq_along(winlens)){
    dists[[i]] = dplyr::filter(tbl, winlen==winlens[i])$JSE
  }

  for(i in base::seq_len(base::length(winlens))){
    x <- base::vector()
    for(j in base::seq_len(base::length(winlens))){
      if(i!=j){
        x <- c(x, stats::t.test(dists[[i]], dists[[j]])$p.value < 0.05)
      }
    }
    if(base::sum(x)>1){
      break
    }
  }
  return(winlens[i])
}
