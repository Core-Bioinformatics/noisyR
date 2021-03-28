#' Function to calculate the noise threshold for a given expression matrix and parameters
#' @description This function is used to calculate the noise threshold for a given expression matrix.
#' It uses as input an expression profile (as calculated by calculate_expression_similarity_*()),
#' or just an expression matrix for a simple calculation based on density.
#' A variety of methods are available to obtain a noise threshold using an input similarity threshold.
#' @param expression either an expression summary (as calculated by calculate_expression_similarity_*()),
#' which should be a list with 3 slots: expression.matrix, expression.levels, expression.levels.similarity;
#' alternatively, just an expression matrix; only density based methods are available for the latter case
#' @param similarity.threshold similarity (correlation or inverse distance) threshold to be used
#' to find corresponding noise threshold; the default, 0.25 is usually suitable for the
#' Pearson correlation (the default similarity measure)
#' @param method.chosen method to use to obtain a vector of noise thresholds,
#' must be one of get_methods_calculate_noise_threshold(); defaults to Boxplot-IQR
#' @param binsize size of each bin in the boxplot methods; defaults to 0.1 (on a log-scale)
#' @param minimum.observations.per.bin minumum number of observations allowed in each bin of the boxplot;
#' if a bin has fewer observations, it is merged with the one to its left; default is calculated as:
#' ceiling(number of observations / number of bins / 10)
#' @return The output is a vector of noise thresholds, the same length as the number of columns in
#' the expression matrix, or a single value in the case of density based methods.
#' @export
#' @examples
#' expression.summary <- calculate_expression_similarity_counts(
#'     expression.matrix = matrix(1:100, ncol=5),
#'     method = "correlation_pearson",
#'     n.elements.per.window = 3)
#' calculate_noise_threshold_base(expression.summary)
calculate_noise_threshold_base <- function(
  expression,
  similarity.threshold=0.25,
  method.chosen="Boxplot-IQR",
  binsize=0.1,
  minimum.observations.per.bin=NULL
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
  }else if(base::strsplit(method.chosen, "-")[[1]][1] != "Density_based" &
           base::is.null(expression.levels)){
    stop("Only density based methods are available for a simple matrix input")
  }

  base::message("Calculating noise thresholds for ", base::ncol(expression.matrix), " samples...")
  base::message("    similarity.threshold = ", similarity.threshold)
  base::message("    method.chosen = ", method.chosen)

  # Using the density to obtain a fixed threshold
  approach <- "Density_based"

  # No normalisation
  method <- "No_normalisation"
  if(method.chosen == base::paste(approach, method, sep="-")){
    noise.thresholds <- noisyr::calculate_first_minimum_density(expression.matrix)
  }

  # RPM normalisation
  method <- "RPM_normalisation"
  if(method.chosen==base::paste(approach, method, sep="-")){
    expression.matrix.normalised <- expression.matrix/
      base::colSums(expression.matrix)*
      stats::median(base::colSums(expression.matrix))
    noise.thresholds <- noisyr::calculate_first_minimum_density(expression.matrix.normalised)
  }

  # Quantile normalisation
  method <- "Quantile_normalisation"
  if(base::is.null(method.chosen) ||
     method.chosen==base::paste(approach, method, sep="-")){
    expression.matrix.normalised <- preprocessCore::normalize.quantiles(expression.matrix)
    rownames(expression.matrix.normalised) <- base::rownames(expression.matrix)
    colnames(expression.matrix.normalised) <- base::colnames(expression.matrix)
    noise.thresholds <- noisyr::calculate_first_minimum_density(expression.matrix.normalised)
  }

  # Using the line plot
  approach <- "Line_plot"

  # No smoothing
  method <- "No_smoothing"
  if(method.chosen==base::paste(approach, method, sep="-")){
    similarity.vector  = base::rep(0,base::ncol(expression.matrix))
    noise.thresholds = base::rep(0,base::ncol(expression.matrix))
    for(j in 1:base::ncol(expression.matrix)){
      similarity.threshold.raw = expression.levels.similarity[,j] > similarity.threshold
      similarity.vector[j] = base::max(base::which(!similarity.threshold.raw))
      noise.thresholds[j] = expression.levels[similarity.vector[j],j]
    }
  }

  # loess10 smoothing
  method <- "loess10_smoothing"
  if(method.chosen==base::paste(approach, method, sep="-")){
    similarity.vector  = base::rep(0,base::ncol(expression.matrix))
    noise.thresholds = base::rep(0,base::ncol(expression.matrix))
    for(j in 1:base::ncol(expression.matrix)){
      loessMod10 <- stats::loess(expression.levels.similarity[,j] ~ log2(expression.levels[,j]+1), span=0.10)
      smoothedx10 <- 2^loessMod10$x-1
      smoothedy10 <- stats::predict(loessMod10)
      similarity.threshold.10  = smoothedy10 > similarity.threshold
      similarity.vector[j] = base::max(base::which(similarity.threshold.10  == FALSE))
      noise.thresholds[j] = smoothedx10[similarity.vector[j]]
    }
  }

  # loess25 smoothing
  method <- "loess25_smoothing"
  if(method.chosen==base::paste(approach, method, sep="-")){
    similarity.vector  = base::rep(0,base::ncol(expression.matrix))
    noise.thresholds = base::rep(0,base::ncol(expression.matrix))
    for(j in 1:base::ncol(expression.matrix)){
      loessMod25 <- stats::loess(expression.levels.similarity[,j] ~ log2(expression.levels[,j]+1), span=0.25)
      smoothedx25 <- 2^loessMod25$x-1
      smoothedy25 <- stats::predict(loessMod25)
      similarity.threshold.25  = smoothedy25 > similarity.threshold
      similarity.vector[j] = base::max(base::which(similarity.threshold.25  == FALSE))
      noise.thresholds[j] = smoothedx25[similarity.vector[j]]
    }
  }

  # loess50 smoothing
  method <- "loess50_smoothing"
  if(method.chosen==base::paste(approach, method, sep="-")){
    similarity.vector  = base::rep(0,base::ncol(expression.matrix))
    noise.thresholds = base::rep(0,base::ncol(expression.matrix))
    for(j in 1:base::ncol(expression.matrix)){
      loessMod50 <- stats::loess(expression.levels.similarity[,j] ~ log2(expression.levels[,j]+1), span=0.50)
      smoothedx50 <- 2^loessMod50$x-1
      smoothedy50 <- stats::predict(loessMod50)
      similarity.threshold.50  = smoothedy50 > similarity.threshold
      similarity.vector[j] = base::max(base::which(similarity.threshold.50 == FALSE))
      noise.thresholds[j] = smoothedx50[similarity.vector[j]]
    }
  }

  # Using the boxplot
  approach <- "Boxplot"
  if(base::substr(method.chosen, 1, base::nchar(approach))==approach){
    if(base::is.null(minimum.observations.per.bin)){
      minimum.observations.per.bin <- base::ceiling(
        base::nrow(expression.levels) /
          (base::log2(base::max(expression.levels+1)) / binsize) / 10
      )
    }
    noise.thresholds.median = base::rep(0,base::ncol(expression.matrix))
    noise.thresholds.iqrs25 = base::rep(0,base::ncol(expression.matrix))
    noise.thresholds.quant5 = base::rep(0,base::ncol(expression.matrix))
    for(j in base::seq_len(base::ncol(expression.matrix))){
      df <- tibble::tibble(x=base::log2(expression.levels[,j]+1),
                           y=expression.levels.similarity[,j])
      breaks <- base::seq(0, base::ceiling(base::max(df$x)), binsize)
      brk <- 2
      while(brk < base::length(breaks)){
        if(base::sum(df$x>=breaks[brk] & df$x<breaks[brk+1]) < minimum.observations.per.bin){
          breaks <- breaks[base::seq_len(base::length(breaks)) != brk]
        }else{
          brk <- brk + 1
        }
      }
      x.binned <- base::cut(df$x, breaks, right=FALSE)
      df <- dplyr::mutate(df, x.binned=x.binned)
      medians <- base::rep(0, base::length(base::levels(df$x.binned)))
      iqrs25 <- base::rep(0, base::length(base::levels(df$x.binned)))
      quant5 <- base::rep(0, base::length(base::levels(df$x.binned)))
      for(l in base::seq_len(base::length(base::levels(df$x.binned)))){
        df.sub <- dplyr::filter(df, df$x.binned==base::levels(df$x.binned)[l])
        medians[l] = stats::quantile(df.sub$y, 1/2, na.rm=TRUE)
        iqrs25[l] = stats::quantile(df.sub$y, 1/4, na.rm=TRUE)
        quant5[l] = stats::quantile(df.sub$y, 1/20, na.rm=TRUE)
      }
      medians[base::is.na(medians)] <- 1
      iqrs25[base::is.na(iqrs25)] <- 1
      quant5[base::is.na(quant5)] <- 1

      noise.thresholds.median[j] = 2^breaks[base::max(base::match(
        medians[medians<similarity.threshold][base::sum(medians<similarity.threshold)],medians) + 1, 1)]
      noise.thresholds.iqrs25[j] = 2^breaks[base::max(base::match(
        iqrs25[iqrs25<similarity.threshold][base::sum(iqrs25<similarity.threshold)],iqrs25) + 1, 1)]
      noise.thresholds.quant5[j] = 2^breaks[base::max(base::match(
        quant5[quant5<similarity.threshold][base::sum(quant5<similarity.threshold)],quant5) + 1, 1)]
    }

    # Using the median
    method <- "Median"
    if(method.chosen==base::paste(approach, method, sep="-")){
      noise.thresholds <- noise.thresholds.median
    }

    # Using the inter-quartile range (IQR)
    method <- "IQR"
    if(method.chosen==base::paste(approach, method, sep="-")){
      noise.thresholds <- noise.thresholds.iqrs25
    }

    # Using the 5& - 95% quantile range
    method <- "Quant5"
    if(method.chosen==base::paste(approach, method, sep="-")){
      noise.thresholds <- noise.thresholds.quant5
    }
  }

  return(noise.thresholds)
}
