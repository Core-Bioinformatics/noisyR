#' Function to calculate the noise threshold for a given expression matrix
#' @description This function is used to calculate the noise threshold for a given expression matrix.
#' It uses as input a distance matrix and the corresponding abundance matrix.
#' A variety of methods are available to obtain an abundance threshold using an input distance threshold.
#' @param expression.matrix expression matrix, should be the one used in calculate_distance_matrices()
#' @param abn.matrix,dist.matrix the input distance and abundance matrices as calculated by
#' calculate_distance_matrices(); if either is not supplied, only a fixed threshold is
#' calculated based on the density
#' @param dist.thresh correlation threshold to be used to find corresponding abundance threshold.
#' The default, 0.25 is usually suitable for the Pearson correlation (the default method)
#' @param binsize size of each bin in the boxplot methods; defaults to 0.1 (on a log-scale)
#' @param min.pts.in.box minumum number of points allowed in each box in the boxplot;
#' if a box has fewer observations, it is merged with the one to its left; default is 20
#' @param dump.stats name of csv to export different thresholds calculated (optional)
#' @param method.chosen method to use to obtain a single vector of thresholds,
#' must be one of get_methods_calculate_noise_threshold();
#' if set, it skips all other methods; this is meant for speed, to be used internally or if rerunning
#' an analysis and is not recommended as a first approach.
#' @return Normal output is a tibble containing information on thresholds calculated using different
#' methods (returned silently and optionally written in a csv file).
#' If method.chosen is set to one of the methods in get_methods_calculate_noise_threshold(),
#' then the output is a vector of noise thresholds,
#' the same length as the number of columns in the expression matrix.
#' @export
#' @examples obj <- calculate_distance_matrices_counts(
#'     expression.matrix = matrix(1:100, ncol=5),
#'     method="correlation_pearson",
#'     n.elements.per.window=3)
#' calculate_threshold_noise(obj$exp, obj$abn, obj$dist,
#'     method.chosen="Boxplot-IQR")
calculate_threshold_noise <- function(
  expression.matrix, abn.matrix=NULL, dist.matrix=NULL,
  dist.thresh=0.25, binsize=0.1, min.pts.in.box=20,
  dump.stats=NULL, method.chosen=NULL)
{
  stats.df <- tibble::tibble()

  # Using the density to obtain a fixed threshold
  approach <- "Density_based_fixed_threshold"

  # No normalisation
  method <- "No_normalisation"
  if(base::is.null(method.chosen) ||
     method.chosen==base::paste(approach, method, sep="-")){
    mat <- expression.matrix
    thresholds.vec <- base::rep(noisyr::calculate_threshold_fixed_density(mat),
                                base::ncol(expression.matrix))
    stats.vec <- tibble::tibble(approach=approach,
                                method=method,
                                dist.thresh="N/A",
                                abn.thresh.min=base::min(thresholds.vec),
                                abn.thresh.mean=base::mean(thresholds.vec),
                                abn.thresh.coef.var=stats::sd(thresholds.vec)/
                                  base::mean(thresholds.vec),
                                abn.thresh.max=base::max(thresholds.vec),
                                abn.thresh.all=base::paste(thresholds.vec,
                                                           collapse=","))
    stats.df <- base::rbind(stats.df, stats.vec)
  }

  # RPM normalisation
  method <- "RPM_normalisation"
  if(base::is.null(method.chosen) ||
     method.chosen==base::paste(approach, method, sep="-")){
    mat <- expression.matrix/
      base::colSums(expression.matrix)*
      stats::median(base::colSums(expression.matrix))
    thresholds.vec <- base::rep(noisyr::calculate_threshold_fixed_density(mat),
                                base::ncol(expression.matrix))
    stats.vec <- tibble::tibble(approach=approach,
                                method=method,
                                dist.thresh="N/A",
                                abn.thresh.min=base::min(thresholds.vec),
                                abn.thresh.mean=base::mean(thresholds.vec),
                                abn.thresh.coef.var=stats::sd(thresholds.vec)/
                                  base::mean(thresholds.vec),
                                abn.thresh.max=base::max(thresholds.vec),
                                abn.thresh.all=base::paste(thresholds.vec,
                                                           collapse=","))
    stats.df <- base::rbind(stats.df, stats.vec)
  }

  # Quantile normalisation
  method <- "Quantile_normalisation"
  if(base::is.null(method.chosen) ||
     method.chosen==base::paste(approach, method, sep="-")){
    mat <- preprocessCore::normalize.quantiles(expression.matrix)
    rownames(mat) <- base::rownames(expression.matrix)
    colnames(mat) <- base::colnames(expression.matrix)
    thresholds.vec <- base::rep(noisyr::calculate_threshold_fixed_density(mat),
                                base::ncol(expression.matrix))
    stats.vec <- tibble::tibble(approach=approach,
                                method=method,
                                dist.thresh="N/A",
                                abn.thresh.min=base::min(thresholds.vec),
                                abn.thresh.mean=base::mean(thresholds.vec),
                                abn.thresh.coef.var=stats::sd(thresholds.vec)/
                                  base::mean(thresholds.vec),
                                abn.thresh.max=base::max(thresholds.vec),
                                abn.thresh.all=base::paste(thresholds.vec,
                                                           collapse=","))
    stats.df <- base::rbind(stats.df, stats.vec)
  }

  if(!(is.null(dist.matrix) | is.null(abn.matrix))){

    # Using the line plot
    approach <- "Line_plot"

    # No smoothing
    method <- "No_smoothing"
    if(base::is.null(method.chosen) ||
       method.chosen==base::paste(approach, method, sep="-")){
      dist.thr  = base::rep(0,base::ncol(expression.matrix))
      abn.thr = base::rep(0,base::ncol(expression.matrix))
      for(j in 1:base::ncol(expression.matrix)){
        dist.thr.raw = dist.matrix[,j] > dist.thresh
        dist.thr[j] = base::max(base::which(dist.thr.raw == FALSE))
        abn.thr[j] = abn.matrix[dist.thr[j],j]
      }
      thresholds.vec <- abn.thr
      stats.vec <- tibble::tibble(approach=approach,
                                  method=method,
                                  dist.thresh=dist.thresh,
                                  abn.thresh.min=base::min(thresholds.vec),
                                  abn.thresh.mean=base::mean(thresholds.vec),
                                  abn.thresh.coef.var=stats::sd(thresholds.vec)/
                                    base::mean(thresholds.vec),
                                  abn.thresh.max=base::max(thresholds.vec),
                                  abn.thresh.all=base::paste(thresholds.vec,
                                                             collapse=","))
      stats.df <- base::rbind(stats.df, stats.vec)
    }

    # loess10 smoothing
    method <- "loess10_smoothing"
    if(base::is.null(method.chosen) ||
       method.chosen==base::paste(approach, method, sep="-")){
      dist.thr  = base::rep(0,base::ncol(expression.matrix))
      abn.thr = base::rep(0,base::ncol(expression.matrix))
      for(j in 1:base::ncol(expression.matrix)){
        loessMod10 <- stats::loess(dist.matrix[,j] ~ log2(abn.matrix[,j]+1), span=0.10)
        smoothedx10 <- 2^loessMod10$x-1
        smoothedy10 <- stats::predict(loessMod10)
        dist.thr.10  = smoothedy10 > dist.thresh
        dist.thr[j] = base::max(base::which(dist.thr.10  == FALSE))
        abn.thr[j] = smoothedx10[dist.thr[j]]
      }
      thresholds.vec <- abn.thr
      stats.vec <- tibble::tibble(approach=approach,
                                  method=method,
                                  dist.thresh=dist.thresh,
                                  abn.thresh.min=base::min(thresholds.vec),
                                  abn.thresh.mean=base::mean(thresholds.vec),
                                  abn.thresh.coef.var=stats::sd(thresholds.vec)/
                                    base::mean(thresholds.vec),
                                  abn.thresh.max=base::max(thresholds.vec),
                                  abn.thresh.all=base::paste(thresholds.vec,
                                                             collapse=","))
      stats.df <- base::rbind(stats.df, stats.vec)
    }

    # loess25 smoothing
    method <- "loess25_smoothing"
    if(base::is.null(method.chosen) ||
       method.chosen==base::paste(approach, method, sep="-")){
      dist.thr  = base::rep(0,base::ncol(expression.matrix))
      abn.thr = base::rep(0,base::ncol(expression.matrix))
      for(j in 1:base::ncol(expression.matrix)){
        loessMod25 <- stats::loess(dist.matrix[,j] ~ log2(abn.matrix[,j]+1), span=0.25)
        smoothedx25 <- 2^loessMod25$x-1
        smoothedy25 <- stats::predict(loessMod25)
        dist.thr.25  = smoothedy25 > dist.thresh
        dist.thr[j] = base::max(base::which(dist.thr.25  == FALSE))
        abn.thr[j] = smoothedx25[dist.thr[j]]
      }
      thresholds.vec <- abn.thr
      stats.vec <- tibble::tibble(approach=approach,
                                  method=method,
                                  dist.thresh=dist.thresh,
                                  abn.thresh.min=base::min(thresholds.vec),
                                  abn.thresh.mean=base::mean(thresholds.vec),
                                  abn.thresh.coef.var=stats::sd(thresholds.vec)/
                                    base::mean(thresholds.vec),
                                  abn.thresh.max=base::max(thresholds.vec),
                                  abn.thresh.all=base::paste(thresholds.vec,
                                                             collapse=","))
      stats.df <- base::rbind(stats.df, stats.vec)
    }

    # loess50 smoothing
    method <- "loess50_smoothing"
    if(base::is.null(method.chosen) ||
       method.chosen==base::paste(approach, method, sep="-")){
      dist.thr  = base::rep(0,base::ncol(expression.matrix))
      abn.thr = base::rep(0,base::ncol(expression.matrix))
      for(j in 1:base::ncol(expression.matrix)){
        loessMod50 <- stats::loess(dist.matrix[,j] ~ log2(abn.matrix[,j]+1), span=0.50)
        smoothedx50 <- 2^loessMod50$x-1
        smoothedy50 <- stats::predict(loessMod50)
        dist.thr.50  = smoothedy50 > dist.thresh
        dist.thr[j] = base::max(base::which(dist.thr.50 == FALSE))
        abn.thr[j] = smoothedx50[dist.thr[j]]
      }
      thresholds.vec <- abn.thr
      stats.vec <- tibble::tibble(approach=approach,
                                  method=method,
                                  dist.thresh=dist.thresh,
                                  abn.thresh.min=base::min(thresholds.vec),
                                  abn.thresh.mean=base::mean(thresholds.vec),
                                  abn.thresh.coef.var=stats::sd(thresholds.vec)/
                                    base::mean(thresholds.vec),
                                  abn.thresh.max=base::max(thresholds.vec),
                                  abn.thresh.all=base::paste(thresholds.vec,
                                                             collapse=","))
      stats.df <- base::rbind(stats.df, stats.vec)
    }

    # Using the boxplot
    approach <- "Boxplot"
    if(base::is.null(method.chosen) ||
       base::substr(method.chosen,1,base::nchar(approach))==approach){
      abn.thr.med = base::rep(0,base::ncol(expression.matrix))
      abn.thr.iqr = base::rep(0,base::ncol(expression.matrix))
      for(j in base::seq_len(base::ncol(expression.matrix))){
        df <- tibble::tibble(x=base::log2(abn.matrix[,j]+1),
                             y=dist.matrix[,j])
        breaks <- base::seq(0, base::ceiling(base::max(df$x)), binsize)
        brk <- 2
        while(brk < base::length(breaks)){
          if(base::sum(df$x>=breaks[brk] & df$x<breaks[brk+1]) < min.pts.in.box){
            breaks <- breaks[base::seq_len(base::length(breaks)) != brk]
          }else{
            brk <- brk + 1
          }
        }
        x.binned <- base::cut(df$x, breaks, right=FALSE)
        df <-dplyr::mutate(df, x.binned=x.binned)
        medians <- base::rep(0, base::length(base::levels(df$x.binned)))
        iqrs25 <- base::rep(0, base::length(base::levels(df$x.binned)))
        for(l in base::seq_len(base::length(base::levels(df$x.binned)))){
          df.sub <- dplyr::filter(df, df$x.binned==base::levels(df$x.binned)[l])
          medians[l] = stats::quantile(df.sub$y, 1/2, na.rm=TRUE)
          iqrs25[l] = stats::quantile(df.sub$y, 1/4, na.rm=TRUE)
        }
        medians[base::is.na(medians)] <- 1
        iqrs25[base::is.na(iqrs25)] <- 1

        abn.thr.med[j] = 2^breaks[base::max(base::match(
          medians[medians<dist.thresh][base::sum(medians<dist.thresh)],medians) + 1, 1)]
        abn.thr.iqr[j] = 2^breaks[base::max(base::match(
          iqrs25[iqrs25<dist.thresh][base::sum(iqrs25<dist.thresh)],iqrs25) + 1, 1)]
      }

      # Using the median
      method <- "Median"
      if(base::is.null(method.chosen) ||
         method.chosen==base::paste(approach, method, sep="-")){
        thresholds.vec <- abn.thr.med
        stats.vec <- tibble::tibble(approach=approach,
                                    method=method,
                                    dist.thresh=dist.thresh,
                                    abn.thresh.min=base::min(thresholds.vec),
                                    abn.thresh.mean=base::mean(thresholds.vec),
                                    abn.thresh.coef.var=stats::sd(thresholds.vec)/
                                      base::mean(thresholds.vec),
                                    abn.thresh.max=base::max(thresholds.vec),
                                    abn.thresh.all=base::paste(thresholds.vec,
                                                               collapse=","))
        stats.df <- base::rbind(stats.df, stats.vec)
      }

      # Using the inter-quartile range (IQR)
      method <- "IQR"
      if(base::is.null(method.chosen) ||
         method.chosen==base::paste(approach, method, sep="-")){
        thresholds.vec <- abn.thr.iqr
        stats.vec <- tibble::tibble(approach=approach,
                                    method=method,
                                    dist.thresh=dist.thresh,
                                    abn.thresh.min=base::min(thresholds.vec),
                                    abn.thresh.mean=base::mean(thresholds.vec),
                                    abn.thresh.coef.var=stats::sd(thresholds.vec)/
                                      base::mean(thresholds.vec),
                                    abn.thresh.max=base::max(thresholds.vec),
                                    abn.thresh.all=base::paste(thresholds.vec,
                                                               collapse=","))
        stats.df <- base::rbind(stats.df, stats.vec)
      }
    }
  }

  if(base::is.null(method.chosen)){
    if(!base::is.null(dump.stats)){
      utils::write.table(stats.df, file=dump.stats, sep=",", append=FALSE,
                         row.names=FALSE, col.names = TRUE)
    }
    return(stats.df)
  }else{
    return(thresholds.vec)
  }
}
