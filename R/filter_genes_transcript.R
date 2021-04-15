#' Function to filter the gene table for the transcript approach
#' @description This function is used to filter the gene table (usually created with
#' \code{\link{cast_gtf_to_genes}}), only keeping genes above the noise thresholds.
#' It uses as input the gene table (usually containing individual exons),
#' an expression matrix for each of these and a vector of abundance thresholds.
#' This function is used internally by \code{\link{remove_noise_from_bams}} to determine
#' which genes to retain.
#' @param genes a tibble of the exons extracted from the gtf file;
#' (usually the the output of \code{\link{cast_gtf_to_genes}})
#' @param expression.matrix the expression matrix, usually
#' calculated by \code{\link{calculate_expression_similarity_transcript}}
#' @param noise.thresholds a vector of expression thresholds by sample
#' @param filter.by Either "gene" (default) or "exon"; if filter.by="gene", a gene
#' (as determined by its ENSEMBL id) is removed
#' if and only if all of its exons are below the corresponding noise thresholds;
#' if filter.by="exon", then each exon is individually removed
#' if it is below the corresponding noise thresholds.
#' @param ... arguments passed on to other methods
#' @return Returns a filtered tibble of exons, with the noise removed.
#' @export
#' @examples
#' bams <- rep(system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE), 2)
#' genes <- data.frame("id" = 1:2,
#'                     "gene_id" = c("gene1", "gene2"),
#'                     "seqid" = c("seq1", "seq2"),
#'                     "start" = 1,
#'                     "end" = 1600)
#' noise.thresholds <- c(0, 1)
#' expression.summary = calculate_expression_similarity_transcript(
#'   bams = bams,
#'   genes = genes,
#'   mapq.unique = 99
#' )
#' filter_genes_transcript(
#'     genes = genes,
#'     expression.matrix = expression.summary$expression.matrix,
#'     noise.thresholds = noise.thresholds,
#' )
#'
filter_genes_transcript = function(
  genes,
  expression.matrix,
  noise.thresholds,
  filter.by = c("gene", "exon"),
  ...
){
  base::message("  filtering genes using the noise thresholds")
  ngenes <- base::nrow(genes)
  threshold.matrix <- base::matrix(base::rep(noise.thresholds, base::nrow(expression.matrix)),
                                   ncol=base::ncol(expression.matrix), byrow = TRUE)
  above.noise.threshold <- base::as.vector(
    base::rowSums(expression.matrix >= threshold.matrix) > 0)
  if(filter.by[1] == "gene"){
    base::message("    doing filtering by gene")
    all.unique.ensids <- base::unique(genes$gene_id)
    all.unique.ensids.above.noise <-
      sapply(X = all.unique.ensids, USE.NAMES = FALSE, FUN = function(gene){
        base::any(above.noise.threshold[genes$gene_id == gene])
      })
    genes <- genes[genes$gene_id %in% all.unique.ensids[all.unique.ensids.above.noise], ]
  }else if(filter.by[1] == "exon"){
    base::message("    doing filtering by exon")
    genes <- genes[above.noise.threshold, ]
  }else stop("filter.by must be either gene or exon")
  base::message("  kept ", base::nrow(genes), " entries out of ", ngenes)
  return(genes)
}
