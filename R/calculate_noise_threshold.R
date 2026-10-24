#' Function to calculate the noise threshold for a given expression matrix and parameters
#' @description This function is used to calculate the noise threshold for a given expression matrix.
#' It uses as input an expression profile, or just an expression matrix for a simple calculation based on density.
#' A variety of methods are available to obtain a noise threshold using an input similarity threshold.
#' @param expression either an expression summary (as calculated by
#' \code{\link{calculate_expression_similarity_counts}} or \code{\link{calculate_expression_similarity_transcript}}),
#' which should be a list with 3 slots: expression.matrix, expression.levels, expression.levels.similarity;
#' alternatively, just an expression matrix; only density based methods are available for the latter case
#' @param similarity.threshold similarity (correlation or inverse distance) threshold to be used
#' to find corresponding noise threshold; the default, 0.25 is usually suitable for the
#' Pearson correlation (the default similarity measure)
#' @param method.chosen method to use to obtain a vector of noise thresholds,
#' must be one of \code{\link{get_methods_calculate_noise_threshold}}; defaults to Boxplot-IQR
#' @param binsize size of each bin in the boxplot methods; defaults to 0.1 (on a log-scale)
#' @param minimum.observations.per.bin minumum number of observations allowed in each bin of the boxplot;
#' if a bin has fewer observations, it is merged with the one to its left; default is calculated as:
#' ceiling(number of observations / number of bins / 10)
#' @param ... arguments passed on to other methods
#' @return The output is a vector of noise thresholds, the same length as the number of columns in
#' the expression matrix, or a single value in the case of density based methods.
#' @seealso \code{\link{calculate_noise_threshold_method_statistics}}
#' @export
#' @examples
#' expression.summary <- calculate_expression_similarity_counts(
#'     expression.matrix = matrix(1:100, ncol=5),
#'     method = "correlation_pearson",
#'     n.elements.per.window = 3)
#' calculate_noise_threshold(expression.summary)
calculate_noise_threshold <- function(
  expression,
  similarity.threshold=0.25,
  method.chosen="Boxplot-IQR",
  binsize=0.1,
  minimum.observations.per.bin=NULL,
  ...
){
  curves <- noisyr::calculate_noise_threshold_curves(
    expression = expression,
    method.chosen = method.chosen,
    binsize = binsize,
    minimum.observations.per.bin = minimum.observations.per.bin
  )
  base::message("Calculating noise thresholds for ", curves$n.samples, " samples...")
  base::message("    similarity.threshold = ", similarity.threshold)
  base::message("    method.chosen = ", method.chosen)
  return(noisyr::noise_thresholds_from_curves(curves, similarity.threshold))
}

#' Precompute the similarity-independent part of a noise threshold method
#' @description The noise threshold calculation has two parts: an expensive one that only
#' depends on the expression summary and the method (density estimation, loess fits,
#' binning for the boxplots), and a cheap one that reads off the threshold for a given
#' similarity threshold. This function performs the first part, so that many similarity
#' thresholds can be evaluated quickly with \code{\link{noise_thresholds_from_curves}}
#' (e.g. in \code{\link{calculate_noise_threshold_method_statistics}} or interactively).
#' @inheritParams calculate_noise_threshold
#' @return A list of class \code{noisyr_threshold_curves}, holding for each sample the
#' expression levels (x) and similarity summary (y) from which thresholds are read off.
#' @seealso \code{\link{calculate_noise_threshold}}
#' @export
#' @examples
#' expression.summary <- calculate_expression_similarity_counts(
#'     expression.matrix = matrix(1:100, ncol=5),
#'     method = "correlation_pearson",
#'     n.elements.per.window = 3)
#' curves <- calculate_noise_threshold_curves(expression.summary)
#' noise_thresholds_from_curves(curves, 0.25)
calculate_noise_threshold_curves <- function(
  expression,
  method.chosen="Boxplot-IQR",
  binsize=0.1,
  minimum.observations.per.bin=NULL,
  ...
){
  if(base::is.matrix(expression)){
    expression.matrix <- expression
    expression.levels <- NULL
    expression.levels.similarity <- NULL
  }else if(base::is.list(expression) &
           identical(names(expression),
                     c("expression.matrix",
                       "expression.levels",
                       "expression.levels.similarity"))){
    expression.matrix <- expression$expression.matrix
    expression.levels <- expression$expression.levels
    expression.levels.similarity <- expression$expression.levels.similarity
  }else{
    stop("Please provide an expression.matrix or an expression.summary list")
  }

  if(!(method.chosen %in% noisyr::get_methods_calculate_noise_threshold())){
    stop("Please provide a method from get_methods_calculate_noise_threshold()")
  }
  approach <- base::strsplit(method.chosen, "-")[[1]][1]
  method <- base::strsplit(method.chosen, "-")[[1]][2]
  if(approach != "Density_based" & base::is.null(expression.levels)){
    stop("Only density based methods are available for a simple matrix input")
  }

  curves <- base::list(method.chosen = method.chosen,
                       approach = approach,
                       n.samples = base::ncol(expression.matrix),
                       all.na = FALSE)
  base::class(curves) <- "noisyr_threshold_curves"

  if(!base::is.null(expression.levels.similarity) &&
     base::all(base::is.na(base::rowSums(expression.levels.similarity)))){
    curves$all.na <- TRUE
    return(curves)
  }

  # Using the density to obtain a fixed threshold
  if(approach == "Density_based"){
    if(method == "RPM_normalisation"){
      expression.matrix <- expression.matrix/
        base::colSums(expression.matrix)*
        stats::median(base::colSums(expression.matrix))
    }else if(method == "Quantile_normalisation"){
      expression.matrix.normalised <- preprocessCore::normalize.quantiles(expression.matrix)
      rownames(expression.matrix.normalised) <- base::rownames(expression.matrix)
      colnames(expression.matrix.normalised) <- base::colnames(expression.matrix)
      expression.matrix <- expression.matrix.normalised
    }
    curves$fixed <- noisyr::calculate_first_minimum_density(expression.matrix)
  }

  # Using the line plot, optionally loess-smoothed on a log-scale
  if(approach == "Line_plot"){
    span <- base::switch(method,
                         "No_smoothing" = NULL,
                         "loess10_smoothing" = 0.10,
                         "loess25_smoothing" = 0.25,
                         "loess50_smoothing" = 0.50)
    curves$by.sample <- base::lapply(base::seq_len(curves$n.samples), function(j){
      if(base::is.null(span)){
        return(base::list(x = expression.levels[,j], y = expression.levels.similarity[,j]))
      }
      loessMod <- base::suppressWarnings(
        stats::loess(expression.levels.similarity[,j] ~ log2(expression.levels[,j]+1), span=span)
      )
      base::list(x = 2^base::as.vector(loessMod$x)-1, y = stats::predict(loessMod))
    })
  }

  # Using the boxplot; the summary statistic of each bin is the median, 25% or 5% quantile
  if(approach == "Boxplot"){
    if(base::is.null(minimum.observations.per.bin)){
      minimum.observations.per.bin <- base::ceiling(
        base::nrow(expression.levels) /
          (base::log2(base::max(expression.levels+1)) / binsize) / 10
      )
    }
    prob <- base::switch(method, "Median" = 1/2, "IQR" = 1/4, "Quant5" = 1/20)
    curves$by.sample <- base::lapply(base::seq_len(curves$n.samples), function(j){
      x <- base::log2(expression.levels[,j]+1)
      y <- expression.levels.similarity[,j]
      breaks <- base::seq(0, base::max(base::ceiling(base::max(x)), 1), binsize)
      brk <- 2
      while(brk < base::length(breaks)){
        if(base::sum(x>=breaks[brk] & x<breaks[brk+1]) < minimum.observations.per.bin){
          breaks <- breaks[base::seq_len(base::length(breaks)) != brk]
        }else{
          brk <- brk + 1
        }
      }
      x.binned <- base::cut(x, breaks, right=FALSE)
      quants <- base::vapply(base::split(y, x.binned), function(v){
        stats::quantile(v, prob, na.rm=TRUE, names=FALSE)
      }, base::numeric(1))
      quants[base::is.na(quants)] <- 1
      base::list(breaks = breaks, y = base::unname(quants))
    })
  }

  return(curves)
}

#' Read off noise thresholds from precomputed curves
#' @description Fast second part of the noise threshold calculation;
#' see \code{\link{calculate_noise_threshold_curves}}.
#' @param curves the output of \code{\link{calculate_noise_threshold_curves}}
#' @param similarity.threshold similarity (correlation or inverse distance) threshold
#' @return A vector of noise thresholds, one per sample, or a single value
#' for the density based methods, identical to \code{\link{calculate_noise_threshold}}.
#' @export
noise_thresholds_from_curves <- function(curves, similarity.threshold=0.25){
  if(curves$all.na){
    base::message("Similarity calculation produced too many NAs, returning zero...")
    return(base::rep(0, curves$n.samples))
  }
  if(curves$approach == "Density_based"){
    return(curves$fixed)
  }
  if(curves$approach == "Line_plot"){
    # last point whose similarity is not above the threshold
    return(base::vapply(curves$by.sample, function(curve){
      above <- curve$y > similarity.threshold
      if(base::any(base::is.na(above)) | base::all(above)){
        id <- 1
      }else{
        id <- base::max(base::which(!above))
      }
      curve$x[id]
    }, base::numeric(1)))
  }
  # Boxplot: start of the bin after the last bin below the threshold
  return(base::vapply(curves$by.sample, function(curve){
    v <- curve$y
    2^curve$breaks[base::max(base::match(
      v[v<similarity.threshold][base::sum(v<similarity.threshold)], v) + 1, 1)]
  }, base::numeric(1)))
}
