#' Function to find the first local minimum of the density of a vector
#' @description This function is used to estimate the first local minimum of the
#' density of a vector.
#' It is meant to be used on the distribution of expression of genes in a sample;
#' since the distribution tails off, finding the global minimum is not appropriate.
#' The plot option can be used to visualise the process.
#' @param mat matrix whose columns will be used; usually an expression matrix; it
#' can also be a vector
#' @param log.transform whether to log-transform the data before the density estimation;
#' default is TRUE
#' @param adjust adjust factor for the smoothing, passed to density(); default is 2
#' @param makeplots a logical value of whether a plot with a vertical line on the minimum found
#' should be printed for each column of the matrix.
#' @return The function outputs a single value corresponding to the median of the minima calculated
#' for each column of the matrix. floor() is taken as a conservative estimate
#' @export
#' @examples calculate_first_minimum_density(
#'     matrix(c(rep(0,100),rep(3,30),rep(10,50),12,13,15,20),ncol=1),
#'     log.transform=FALSE, makeplots=TRUE)
calculate_first_minimum_density = function(mat, log.transform=TRUE, adjust=2, makeplots=FALSE){
  something.is.wrong <- FALSE
  if(base::is.vector(mat)){mat <- base::as.matrix(mat, ncol=1)}
  firstmin.vector <- base::vector(length=base::ncol(mat))
  for(j in base::seq_len(base::ncol(mat))){
    vec <- mat[,j]
    if(log.transform){vec = base::log2(vec+1)}
    dens <- stats::density(base::sort(vec), adjust=2)
    i <- 2
    while(dens$y[i]==base::max(dens$y[1:i])){
      i = i + 1
      if(i > base::length(dens$y)){
        something.is.wrong <- TRUE
        break
      }
    }
    i = i - 1
    firstmax <- dens$x[i]
    idmax <- base::match(firstmax, dens$x)
    i <- idmax
    while(dens$y[i]==base::min(dens$y[idmax:i])){
      i = i + 1
      if(i > base::length(dens$y)){
        something.is.wrong <- TRUE
        break
      }
    }
    i = i - 1
    if(!something.is.wrong){
      xmain <- dens$x[i:length(dens$x)]
      ymain <- dens$y[i:length(dens$x)]
      secondmax <- xmain[ymain==base::max(ymain)]
      firstmin <- xmain[ymain==base::min(ymain[xmain<secondmax])]
      if(makeplots){
        x=NULL; y=NULL
        print(ggplot2::ggplot(data=base::data.frame("x"=dens$x, "y"=dens$y)) +
                ggplot2::theme_minimal() +
                ggplot2::geom_line(mapping=ggplot2::aes(x=x, y=y)) +
                ggplot2::geom_vline(xintercept=firstmin, color="blue"))
      }
      if(log.transform){
        firstmax = 2^firstmax
        secondmax = 2^secondmax
        firstmin = 2^firstmin
      }
      firstmin.vector[j] <- firstmin
    }
  }
  if(base::sum(firstmin.vector == 0) > base::length(firstmin.vector) / 2){
    base::message("Distribution does not appear to have two peaks, returning zero")
    firstmin.vector.median <- 0
  }else{
    firstmin.vector <- firstmin.vector[firstmin.vector > 0]
    firstmin.vector.median <- base::floor(stats::median(firstmin.vector))
  }

  return(firstmin.vector.median)
}
