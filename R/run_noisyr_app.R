#' Launch the interactive noisyR app
#' @description Opens a Shiny app that walks through the count matrix approach:
#' loading a count matrix, computing the window similarities, choosing a noise
#' threshold method and similarity threshold interactively (with all methods compared
#' side by side), and denoising, with downloads of the denoised matrix and of an
#' R script that reproduces the chosen parameters.
#' @param ... arguments passed on to \code{shiny::runApp}, e.g. \code{port} or \code{launch.browser}
#' @return Called for its side effect of running the app.
#' @export
#' @examples
#' if(interactive()) run_noisyr_app()
run_noisyr_app <- function(...){
  needed <- c("shiny", "bslib", "plotly")
  missing <- needed[!base::vapply(needed, base::requireNamespace, base::logical(1), quietly=TRUE)]
  if(base::length(missing) > 0){
    base::stop("The app needs the following packages: ", base::paste(missing, collapse=", "),
               "; install them with install.packages()")
  }
  shiny::runApp(base::system.file("shiny", "noisyR_app", package="noisyr", mustWork=TRUE), ...)
}
