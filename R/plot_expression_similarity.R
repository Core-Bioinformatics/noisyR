#' Plot the similarity against expression levels
#' @description Creates the expression-similarity line and box plots for each sample.
#' @param expression.summary list containing expression_levels and expression_levels_similarity
#' matrices, as calculated by \code{\link{calculate_expression_similarity_counts}} or
#' \code{\link{calculate_expression_similarity_transcript}}
#' @param sample.names names for the plots, defaults to the column names of the expression matrix
#' @param similarity.name similarity metric used (for the y-axis title)
#' @param log.transform should the count matrix be log-transformed? If not, boxplot is skipped
#' @param min.y,max.y limits for the y axis. If unset default to symmetric including all values
#'                    in expression.levels.similarity; min is set to 0 if there are no negative values
#' @param smooth.span span to be used for smoothing in the line plot; defaults to 0.1
#' @param only.boxplot option to skip the line plot (usually a good idea if there are too many points
#'                     and lines are too erratic); sets log.transform to TRUE
#' @param binsize size of each bin in the boxplot; defaults to 0.5
#' @param last.together groups observations so the highest abundance bin has at least this many
#' @param show.counts whether to show how many observations are in each bin
#' @param add.threshold adds a horizontal line at this value
#' @param file.name name of pdf to output the plots; if not provided (default), no printing is done
#' @return A list of all the plots (returned silently)
#' @export
#' @examples
#' plots <- plot_expression_similarity(
#'   expression.summary=list(
#'     "expression.levels" = matrix(2^(10*seq(0,1,length.out=100))),
#'     "expression.levels.similarity" = matrix(seq(0,1,length.out=100)+(runif(100)/5))))
#' plots[[1]]
#' plots[[2]]
plot_expression_similarity <- function(expression.summary,
                                       sample.names=NULL,
                                       similarity.name="Pearson correlation",
                                       log.transform=TRUE, min.y=NULL, max.y=NULL, smooth.span=0.1,
                                       only.boxplot=FALSE, binsize=1, last.together=0,
                                       show.counts=TRUE, add.threshold=NULL, file.name=NULL)
{
  expression.matrix <- expression.summary$expression.matrix
  expression.levels <- expression.summary$expression.levels
  expression.levels.similarity <- expression.summary$expression.levels.similarity

  if(is.null(sample.names)){
    sample.names <- base::colnames(expression.matrix)
  }

  x=NULL; y=NULL; median=NULL
  if(only.boxplot){log.transform <- TRUE}
  if(log.transform){expression.levels <- base::log2(expression.levels+1)}
  if(base::is.null(min.y) & base::is.null(max.y)){
    max.y <- base::max(base::abs(expression.levels.similarity), na.rm=TRUE)
  }
  if(base::is.null(min.y)){
    if(base::min(expression.levels.similarity, na.rm=TRUE)<0){
      min.y <- -max.y
    }else{
      min.y <- 0
    }
  }
  if(base::is.null(max.y)){
    if(min.y<0){
      max.y <- -min.y
    }else{
      max.y <- base::max(base::abs(expression.levels.similarity), na.rm=TRUE)
    }
  }
  if(!base::is.null(file.name)){grDevices::pdf(file.name, width=7, height=7)}
  all.plots <- list()
  plot.id <- 1
  for(j in base::seq_len(base::ncol(expression.levels)))
  {
    if(!base::is.null(file.name)){
      base::message("Plotting ", sample.names[j])
    }
    df <- tibble::tibble(x=expression.levels[,j], y=expression.levels.similarity[,j])
    if(!only.boxplot){
      plot <- ggplot2::ggplot(df) +
        ggplot2::theme_minimal() +
        ggplot2::geom_line(ggplot2::aes(x,y)) +
        ggplot2::geom_smooth(ggplot2::aes(x,y), method="loess",
                             formula= y ~ x, span=smooth.span) +
        ggplot2::ylim(base::c(min.y, max.y)) +
        ggplot2::ggtitle(sample.names[j]) +
        ggplot2::xlab(base::ifelse(log.transform, "expression (log2)", "expression (linear)")) +
        ggplot2::ylab(similarity.name)
      if(!base::is.null(add.threshold)){
        plot = plot + ggplot2::geom_hline(yintercept=add.threshold, color="green")
      }
      if(!base::is.null(file.name)){base::print(plot)}
      all.plots[[plot.id]] <- plot
      plot.id = plot.id + 1
    }
    if(log.transform){
      n.intervals <- base::ceiling(base::max(df$x))/binsize-1
      while(base::sum(df$x>(binsize*n.intervals))<last.together){
        n.intervals = n.intervals - 1}
      n.intervals = n.intervals + 1
      breaks <- base::unique(base::c(0,
                                     (1:n.intervals)*binsize,
                                     base::ceiling(base::max(df$x))))
      x.binned <- base::cut(df$x, breaks, right=FALSE)
      df <- df %>% dplyr::mutate(x.binned=x.binned)
      plot <- ggplot2::ggplot(df) +
        ggplot2::theme_minimal() +
        ggplot2::geom_boxplot(ggplot2::aes(x.binned, y)) +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=1)) +
        ggplot2::ylim(c(min.y, max.y*1.1)) +
        ggplot2::ggtitle(sample.names[j]) +
        ggplot2::xlab("expression (log2)") +
        ggplot2::ylab(similarity.name)
      if(show.counts){
        plot <- plot +
          ggplot2::stat_summary(ggplot2::aes(x.binned, y),
                                fun.data = function(x) return(base::c(y = max.y*1.05,
                                                                      label = base::length(x))),
                                geom = "text", fun = median)
      }
      if(!base::is.null(add.threshold)){
        plot = plot + ggplot2::geom_hline(yintercept=add.threshold, color="green")
      }
      if(!base::is.null(file.name)){base::print(plot)}
      all.plots[[plot.id]] <- plot
      plot.id = plot.id + 1
    }
  }
  if(!base::is.null(file.name)){grDevices::dev.off()}
  return(base::invisible(all.plots))
}
