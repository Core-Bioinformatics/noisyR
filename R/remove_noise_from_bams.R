#' Function to remove the noisy reads from the BAM files
#' @description This function is used to remove the noisy reads from the BAM files.
#' It uses as input the BAM file names, a gene table (usually containing individual exons,
#' made using cast_gtf_to_genes()), an expression matrix for each of these genes and
#' a vector of abundance thresholds.
#' @param bams a character vector of the BAM file names
#' @param genes a tibble of the exons extracted from the gtf file;
#' (usually the the output of cast_gtf_to_genes())
#' @param expression the expression matrix or expression summary list,
#' as calculated by calculate_expression_similarity_transcript()
#' @param noise.thresholds a vector of expression thresholds by sample;
#' must be the same length as the number of BAM files,
#' or a singular value to be used as a fixed noise threshold
#' @param destination.files names for the output denoised BAM files; by default the same as
#' the original files, appended with ".noisefiltered.bam", but created in the working directory
#' @param filter.by Either "gene" (default) or "exon"; if filter.by="gene", a gene is removed from all BAM files
#' if and only if all of its exons are below the corresponding noise thresholds;
#' if filter.by="exon", then each exon is individually removed (from all samples)
#' if it is below the corresponding noise thresholds.
#' @param make.index whether a BAM index should be generated; if this is FALSE (the default)
#' and no index exists, the function will exit with an error; the index needs to have
#' the same name as each BAM file, but ending with .bam.bai
#' @param unique.only whether only uniquely mapped reads should contribute to the expression
#' of a gene/exon; default is TRUE
#' @param mapq.unique The values of the mapping quality field in the BAM file that corresponds
#' to uniquely mapped reads; by default, values of 255 are used as these correspond to
#' the most popular aligners, but an adjustment might be needed;
#' the mapq scores should be as follows: 255 for STAR, 60 for hisat2,
#' 255 for bowtie in -k mode, 40 for bowtie2 default, 50 for tophat
#' @param ... arguments passed on to other methods
#' @return Returns a matrix of the same dims as the expression matrix, with the noise removed.
#' This matrix has no entries remaining below the noise threshold.
#' @export
#' @examples
#' expression.matrix <- matrix(1:100, ncol=5)
#' noise.thresholds <- c(5,30,45,62,83)
#' expression.matrix.denoised <- remove_noise_from_matrix(
#'     expression.matrix = expression.matrix,
#'     noise.thresholds = noise.thresholds
#' )
#'
remove_noise_from_bams = function(
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
){
  if(base::is.matrix(expression)){
    expression.matrix <- expression
  }else if(base::is.list(expression) &
           identical(names(expression)[1], "expression.matrix")){
    expression.matrix <- expression$expression.matrix
  }else{
    stop("Please provide an expression.matrix or an expression.summary list")
  }

  if(base::length(noise.thresholds) == 1){
    base::message("noise.thresholds only has 1 value, using a fixed threshold...")
    noise.thresholds <- base::rep(noise.thresholds, base::ncol(expression.matrix))
  }else if(base::length(noise.thresholds) != base::ncol(expression.matrix)){
    base::stop("noise.thresholds needs to be length 1 or ncol(expression.matrix)")
  }

  base::message("Denoising ", base::length(bams), " BAM files...")

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
  base::message("Kept ", base::nrow(genes), " entries out of ", ngenes)

  bam.indexes <- base::paste(bams, ".bai", sep="")
  for(j in base::seq_len(base::length(bam.indexes))){
    bam.index <- bam.indexes[j]
    if(!base::file.exists(bam.index)){
      if(make.index){
        base::message("Creating BAM index for ", bams[j])
        bams[j] <- Rsamtools::sortBam(bams[j], destination=base::paste0(bams[j], ".sorted"))
        bam.index <- Rsamtools::indexBam(bams[j])
      }else{
        stop("BAM index not found. It can be generated by setting make.index = TRUE")
      }
    }
  }
  gr <- GenomicRanges::GRanges(seqnames = genes$seqid,
                               ranges = IRanges::IRanges(start=genes$start,
                                                               end=genes$end))
  params <- Rsamtools::ScanBamParam(which = gr,
                                    what = Rsamtools::scanBamWhat(),
                                    mapqFilter = ifelse(unique.only, mapq.unique, 0))
  for(j in base::seq_len(base::length(bams))){
    base::message("  denoising BAM file ", j, " of ", base::length(bams))
    Rsamtools::filterBam(file = bams[j],
                         destination = destination.files[j],
                         indexDestination = make.index,
                         param = params)
  }
  return(base::invisible(NULL))
}
