#' Show the methods for calculating correlation or distance
#' @description This function outputs the methods available for the calculation of
#' the correlation or distance. The standard correlation methods use stats::cor and
#' a wide variety of distance methods are available using the philentropy package.
#' To be used as input in calculate_distance_matrices_*().
#' @param names whether to output names (default) or characterisation
#' as similarity or dissimilarity (used internally to invert dissimilarity measures)
#' @return A character vector of options for the method arguement of calculate_distance_matrices_*();
#' if names=FALSE, a vector of types (similarity/dissimilarity measure) of the same length
#' @export
#' @examples get_methods_correlation_distance()

get_methods_correlation_distance = function(names=TRUE){
  if(names){
    methods.cor <- base::c("pearson", "kendall", "spearman")
    methods.dist <- philentropy::getDistMethods()
    methods.all <- base::c(base::paste0("correlation_", methods.cor),
                           base::paste0("distance_", methods.dist))
    return(methods.all)
  }else{
    is.sim <- base::c(1, 1, 1,                 # correlations
                      0, 0, 0, 0,              # Lp Minkowski Family
                      0, 0, 0, 0, 0, 0,        # L1 Family
                      1, 0, 0, 0, 0, 0, 0, 1,  # Intersection Family
                      1, 1, 1, 1, 0, 0,        # Inner Product Family
                      1, 0, 0, 0, 0,           # Squared-chord Family
                      0, 0, 0, 0, 0, 0, 0, 0,  # Squared L2 family (X2 squared family)
                      0, 0, 0, 0, 0, 0,        # Shannon’s Entropy Family
                      0, 0, 0                  # Combinations
    )
    types.all <- base::rep("d", length(is.sim))
    types.all[is.sim==1] <- "s"
    return(types.all)
  }
}
