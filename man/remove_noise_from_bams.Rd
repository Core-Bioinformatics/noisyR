% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_noise_from_bams.R
\name{remove_noise_from_bams}
\alias{remove_noise_from_bams}
\title{Function to remove the noisy reads from the BAM files}
\usage{
remove_noise_from_bams(
  bams,
  genes,
  expression,
  noise.thresholds,
  destination.files = base::paste0(base::basename(bams), ".noisefiltered.bam"),
  filter.by = c("gene", "exon"),
  make.index = FALSE,
  unique.only = TRUE,
  mapq.unique = 255,
  ...
)
}
\arguments{
\item{bams}{a character vector of the BAM file names}

\item{genes}{a tibble of the exons extracted from the gtf file;
(usually the the output of \code{\link{cast_gtf_to_genes}})}

\item{expression}{the expression matrix or expression summary list,
as calculated by \code{\link{calculate_expression_similarity_transcript}}}

\item{noise.thresholds}{a vector of expression thresholds by sample;
must be the same length as the number of BAM files,
or a singular value to be used as a fixed noise threshold}

\item{destination.files}{names for the output denoised BAM files; by default the same as
the original files, appended with ".noisefiltered.bam", but created in the working directory}

\item{filter.by}{Either "gene" (default) or "exon"; if filter.by="gene", a gene is removed from all BAM files
if and only if all of its exons are below the corresponding noise thresholds;
if filter.by="exon", then each exon is individually removed (from all samples)
if it is below the corresponding noise thresholds.}

\item{make.index}{whether a BAM index should be generated; if this is FALSE (the default)
and no index exists, the function will exit with an error; the index needs to have
the same name as each BAM file, but ending with .bam.bai}

\item{unique.only}{whether only uniquely mapped reads should contribute to the expression
of a gene/exon; default is TRUE}

\item{mapq.unique}{The values of the mapping quality field in the BAM file that corresponds
to uniquely mapped reads; by default, values of 255 are used as these correspond to
the most popular aligners, but an adjustment might be needed;
the mapq scores should be as follows: 255 for STAR, 60 for hisat2,
255 for bowtie in -k mode, 40 for bowtie2 default, 50 for tophat}

\item{...}{arguments passed on to other methods}
}
\value{
Returns a matrix of the same dims as the expression matrix, with the noise removed.
This matrix has no entries remaining below the noise threshold.
}
\description{
This function is used to remove the noisy reads from the BAM files.
It uses as input the BAM file names, a gene table (usually containing individual exons,
made using \code{\link{cast_gtf_to_genes}}), an expression matrix for each of these genes and
a vector of abundance thresholds.
}
\examples{
bams <- rep(system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE), 2)
genes <- data.frame("id" = 1:2,
                    "gene_id" = c("gene1", "gene2"),
                    "seqid" = c("seq1", "seq2"),
                    "start" = 1,
                    "end" = 1600)
noise.thresholds <- c(0, 1)
expression.summary = calculate_expression_similarity_transcript(
  bams = bams,
  genes = genes,
  mapq.unique = 99
)
remove_noise_from_bams(
    bams = bams,
    genes = genes,
    expression = expression.summary,
    noise.thresholds = noise.thresholds,
    destination.files = paste0(tempdir(), "/", basename(bams), ".noisefiltered.bam"),
    mapq.unique = 99
)

}
\seealso{
\code{\link{remove_noise_from_matrix}}
}
