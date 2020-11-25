#' Plot the distance
#' @description Creates the abundance-distance line and box plots for each sample.
#' @param abn.matrix,dist.matrix abundance and distance matrix,
#' as calculated by calculate_distance_matrices()
#' @param sample.names names for the plots, enumerates the samples by default
#' @param distance.name distance metric used (for the y-axis title)
#' @param log.transform should the count matrix be log-transformed? If not, boxplot is skipped
#' @param min.y,max.y limits for the y axis. If unset default to symmetric including all values
#'                    in dist.matrix; min is set to 0 if there are no negative values
#' @param smooth.span span to be used for smoothing in the line plot; defaults to 0.1
#' @param only.boxplot option to skip the line plot (usually a good idea if there are too many points
#'                     and lines are too erratic); sets log.transform to TRUE
#' @param binsize size of each bin in the boxplot; defaults to 0.5
#' @param last.together groups observations so the highest abundance bin has at least this many
#' @param show.counts whether to show how many observations are in each bin
#' @param add.threshold adds a horizontal line at this value
#' @param file.name name of pdf to output the plots (console by default)
#' @return A list of all the plots (returned silently), which are also plotted to the console,
#'         or specified pdf file
#' @export
#' @examples plot_distance_abundance(
#'   abn.matrix=matrix(2^(10*seq(0,1,length.out=100))),
#'   dist.matrix=matrix(seq(0,1,length.out=100)+(runif(100)/5)))
plot_distance_abundance <- function(abn.matrix, dist.matrix,
                               sample.names=paste("Sample", 1:ncol(abn.matrix)),
                               distance.name="Pearson correlation",
                               log.transform=TRUE, min.y=NULL, max.y=NULL, smooth.span=0.1,
                               only.boxplot=FALSE, binsize=0.5, last.together=30,
                               show.counts=TRUE, add.threshold=NULL, file.name=NULL)
{
  x=NULL; y=NULL; median=NULL
  if(only.boxplot){log.transform <- TRUE}
  if(log.transform){abn.matrix <- base::log2(abn.matrix+1)}
  if(base::is.null(min.y) & base::is.null(max.y)){
    max.y <- base::max(base::abs(dist.matrix), na.rm=TRUE)
  }
  if(base::is.null(min.y)){
    if(base::min(dist.matrix, na.rm=TRUE)<0){
      min.y <- -max.y
    }else{
      min.y <- 0
    }
  }
  if(base::is.null(max.y)){
    if(min.y<0){
      max.y <- -min.y
    }else{
      max.y <- base::max(base::abs(dist.matrix), na.rm=TRUE)
    }
  }
  if(!base::is.null(file.name)){grDevices::pdf(file.name, width=7, height=7)}
  all.plots <- list()
  plot.id <- 1
  for(j in 1:base::ncol(abn.matrix))
  {
    base::message("Plotting ", sample.names[j])
    df <- tibble::tibble(x=abn.matrix[,j], y=dist.matrix[,j])
    if(!only.boxplot){
      plot <- ggplot2::ggplot(df) +
        ggplot2::geom_line(ggplot2::aes(x,y)) +
        ggplot2::geom_smooth(ggplot2::aes(x,y), method="loess",
                             formula= y ~ x, span=smooth.span) +
        ggplot2::ylim(base::c(min.y, max.y)) +
        ggplot2::ggtitle(sample.names[j]) +
        ggplot2::xlab(base::ifelse(log.transform, "Abn(log2)", "Abn(linear)")) +
        ggplot2::ylab(distance.name)
      if(!base::is.null(add.threshold)){
        plot = plot + ggplot2::geom_hline(yintercept=add.threshold, color="green")
      }
      base::print(plot)
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
      df %>% dplyr::mutate(x.binned=x.binned)
      plot <- ggplot2::ggplot(df) +
        ggplot2::geom_boxplot(ggplot2::aes(x.binned, y)) +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=1)) +
        ggplot2::ylim(c(min.y, max.y*1.1)) +
        ggplot2::ggtitle(sample.names[j]) +
        ggplot2::xlab(base::ifelse(log.transform, "Abn(log2)", "Abn(linear)")) +
        ggplot2::ylab(distance.name)
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
      base::print(plot)
      all.plots[[plot.id]] <- plot
      plot.id = plot.id + 1
    }
  }
  if(!base::is.null(file.name)){grDevices::dev.off()}
  return(base::invisible(all.plots))
}
