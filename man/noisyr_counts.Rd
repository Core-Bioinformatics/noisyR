% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/noisyr_counts.R
\name{noisyr_counts}
\alias{noisyr_counts}
\title{Run the noisyR pipeline for the count matrix approach}
\usage{
noisyr_counts(
  expression.matrix,
  n.elements.per.window = NULL,
  optimise.window.length.logical = FALSE,
  similarity.threshold = 0.25,
  method.chosen = "Boxplot-IQR",
  ...
)
}
\arguments{
\item{expression.matrix}{the expression matrix used as input for the similarity calculation;
this argument is required}

\item{n.elements.per.window}{number of elements to have in a window passed to
calculate_expression_similarity_counts(); default 10\% of the number of rows}

\item{optimise.window.length.logical}{whether to call optimise_window_length to try and
optimise the value of n.elements.per.window}

\item{similarity.threshold, method.chosen}{parameters passed on to
\code{\link{calculate_noise_threshold}}; they can be single values or vectors;
if they are vectors optimal values are computed by calling
\code{\link{calculate_noise_threshold_method_statistics}} and
minimising the coefficient of variation across samples; all possible values for
method.chosen can be viewed by \code{\link{get_methods_calculate_noise_threshold}}}

\item{...}{arguments to be passed on to individual pipeline steps}
}
\value{
The denoised expression matrix.
}
\description{
Calls the functions to run each of the three steps of the pipeline
(similarity calculation, noise quantification, noise removal), with the specified parameters.
See the individual function documentation for more details and required arguments.
Required steps: \code{\link{calculate_expression_similarity_counts}},
\code{\link{calculate_noise_threshold}}. \code{\link{remove_noise_from_matrix}}.
Optional steps: \code{\link{optimise_window_length}},
\code{\link{calculate_noise_threshold_method_statistics}}
}
\examples{
noisyr_counts(
    expression.matrix = matrix(1:100, ncol = 5),
    similarity.measure = "correlation_pearson",
    n.elements.per.window = 3)
}
\seealso{
\code{\link{noisyr}}, \code{\link{noisyr_transcript}}
}
