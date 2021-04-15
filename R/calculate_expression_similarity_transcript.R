#' Calcualte the distance matrices using the BAM files
#' @description This function generates an average correlation/distance coefficient
#' for every exon present in the BAM files. This is done by calculating
#' the point-to-point correlation/distance of the distribution of reads
#' across the transcript of each exon and comparing it across samples.
#' The reason why exons are used instead of full length genes is that long
#' intronic regions artificially increase the correlation since there is
#' consistently no expression there, across samples. The user has the
#' option to use genes instead, by running \code{\link{cast_gtf_to_genes}} separately,
#' with non default parameters.
#' @param bams,path.bams either a path to the directory where the BAM files are
#' or a vector of paths to each individual file; if a path is specified,
#' it extracts all files that end in .bam; looks in the working directory by default
#' @param genes a tibble of the exons extracted from the gtf file;
#' this is meant for speed if the output of \code{\link{cast_gtf_to_genes}} is already generated,
#' or if the user wants to only calculate similarity for a subset of exons
#' @param path.gtf the path to the gtf/gff annotation file (only used if genes is not
#' provided); if unspecified, looks for one in the working directory
#' @param expression.matrix expression matrix; not necessary but is used to filter the
#' gtf to fewer entries and for subsampling if subsample.genes=TRUE;
#' if not provided, raw read counts per exon are extracted from the BAM files
#' @param subsample.genes logical, whether to subsample low abundance genes to decrease
#' computational time; the first minimum of the distribution of abundances is calculated,
#' and genes lower than it are subsampled to match the number of genes higher than it;
#' the expression matrix needs to be provided for this calculation;
#' a plot is generated to show that minimum
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
#' @param slack slack needs to be >=readLength, adjust for efficiency; the default is 200,
#' as it is higher than most modern sequencing experiments
#' @param similarity.measure one of the similarity metrics to be used, defaults to pearson correlation;
#' currently, only correlation is supported
#' @param save.image.every.1000 whether to save a workspace image after every 1000 exons
#' are processed; default is FALSE
#' @param ncores Number of cores for parallel computation; defaults to sequential computation,
#' but parallelisation is highly encouraged; it is set to detectCores() if higher
#' @param ... arguments passed on to other methods
#' @return A list with three elements: the first element is the expression matrix,
#'         as supplied or calculated; the other two are the expression levels matrix and
#'         expression levels similarity matrix;
#'         they have the same # of columns as the expression matrix,
#'         and as many rows as exons processed.
#' @seealso \code{\link{calculate_expression_similarity_counts}}
#' @export
#' @examples
#' bams <- rep(system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE), 2)
#' genes <- data.frame("id" = 1:2,
#'                     "gene_id" = c("gene1", "gene2"),
#'                     "seqid" = c("seq1", "seq2"),
#'                     "start" = 1,
#'                     "end" = 1600)
#' expression.summary <- calculate_expression_similarity_transcript(
#'   bams = bams,
#'   genes = genes,
#'   mapq.unique = 99
#' )

calculate_expression_similarity_transcript <- function(
  bams=NULL,
  path.bams=".",
  genes=NULL,
  path.gtf=list.files(".", pattern="\\.g[tf]f$"),
  expression.matrix=NULL,
  subsample.genes=FALSE,
  make.index=FALSE,
  unique.only=TRUE,
  mapq.unique=255,
  slack=200,
  similarity.measure = "correlation_pearson",
  save.image.every.1000=FALSE,
  ncores=1,
  ...
){
  if(base::is.null(bams)){
    bams <- base::list.files(path.bams, pattern="\\.bam$", full.names=TRUE)
  }
  if(base::length(bams)<2) base::stop("Please provide at least 2 BAM files")
  for(j in base::seq_len(base::length(bams))){
    bam.index <- base::paste0(bams[j], ".bai")
    if(!base::file.exists(bam.index)){
      if(make.index){
        base::message("Creating BAM index for ", bams[j])
        bams[j] <- Rsamtools::sortBam(bams[j], destination=base::paste0(bams[j], ".sorted"))
        bam.index <- Rsamtools::indexBam(bams[j])
      }else{
        stop("BAM index not found. It can be generated by setting make.index=TRUE")
      }
    }
  }
  if(base::is.null(genes)){
    base::message("Creating gene table from gtf file...")
    genes <- noisyr::cast_gtf_to_genes(path.gtf, ...)
  }
  if(base::nrow(genes)<2) base::stop("Please provide at least 2 genes")
  if(!base::is.null(expression.matrix)){
    genes.subset <- genes[genes$gene_id %in% base::rownames(expression.matrix), ]
    if(subsample.genes){
      rSum <- base::rowSums(expression.matrix)
      densSumMin <- noisyr::calculate_first_minimum_density(rSum)
      subsampled.gene.ids <-
        c(base::names(rSum)[rSum>=densSumMin],
          base::names(rSum[base::sample(base::names(rSum)[rSum<densSumMin],
                                        size=base::sum(rSum>=densSumMin))]))
      genes.subset <- genes.subset[genes.subset$gene_id %in% subsampled.gene.ids,]
    }
  }else{
    genes.subset <- genes
  }

  ngenes <- base::nrow(genes.subset)

  use.corr.dist <- base::strsplit(similarity.measure, "_")[[1]][1]
  if(use.corr.dist!="correlation")
    stop(paste("Distance measures are currently not supported for transcript approach.",
               "Please use a correlation measure instead."))
  base.method <- base::sub(paste0(use.corr.dist,"_"), "", similarity.measure)

  expression.levels <- base::matrix(nrow=ngenes, ncol=base::length(bams))
  base::rownames(expression.levels) <- genes.subset$gene_id
  expression.levels.similarity <- base::matrix(nrow=ngenes,
                                               ncol=base::length(bams))

  base::message("Calculating expression similarity for ", ngenes, " genes...")
  base::message("    this process may take a long time...")
  start_time <- base::Sys.time()
  if(ncores>1){
    ncores <- base::min(ncores, parallel::detectCores())
    base::message("    ncores=", ncores)

    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl = cl,
                                   cores = ncores)
    concatenated.matrices <- foreach::foreach(
      n=1:ngenes,
      .combine=rbind,
      .inorder=TRUE) %dopar% {
        profile_expression_list <-
          noisyr::calculate_expression_profile(gene=genes.subset[n,],
                                               bams=bams,
                                               unique.only=unique.only,
                                               mapq.unique=mapq.unique,
                                               slack=slack)
        profile <- profile_expression_list$profile
        cors <- base::suppressWarnings(
          stats::cor(base::unname(profile), method=base.method))
        cors[cors==1] = NA
        avCors <- base::colMeans(cors, na.rm=TRUE)
        avCors[base::is.nan(avCors)] = NA
        expression.levels.similarity[n,] <- avCors
        if(n%%1000==0){
          if(save.image.every.1000){
            base::save.image()
          }
        }
        c(avCors, profile_expression_list$expression)
      }
    parallel::stopCluster(cl)

    expression.levels[] <- concatenated.matrices[, (base::length(bams)+1):(2*base::length(bams))]
    expression.levels.similarity[] <- concatenated.matrices[, 1:base::length(bams)]

  }else{
    base::message("    ncores=1, running sequentially...")
    if(ngenes > 5000){
      base::message("    consider increasing ncores for more than a few thousand exons")
      base::message("    this process may take a long time...")
    }
    for(n in base::seq_len(ngenes)){
      profile_expression_list <-
        noisyr::calculate_expression_profile(gene=genes.subset[n,],
                                             bams=bams,
                                             unique.only=unique.only,
                                             mapq.unique=mapq.unique,
                                             slack=slack)
      profile <- profile_expression_list$profile

      cors <- stats::cor(base::unname(profile), method=base.method)
      cors[cors==1] = NA
      avCors <- base::colMeans(cors, na.rm=TRUE)
      avCors[base::is.nan(avCors)] = NA

      expression.levels[n,] <- profile_expression_list$expression
      expression.levels.similarity[n,] = avCors
      if(n%%1000==0){
        if(save.image.every.1000){
          base::save.image()
        }
        part_time <- base::Sys.time()
        base::message("Done ", n, " genes out of ", ngenes)
        time_elapsed <- part_time - start_time
        base::message(" Time elapsed: ", base::round(time_elapsed, 2),
                      " ", base::units(time_elapsed))
      }
    }
  }
  if(save.image.every.1000){
    base::save.image()
  }
  end_time <- base::Sys.time()
  time_elapsed <- end_time - start_time
  base::message("Finished! Time elapsed: ", base::round(time_elapsed, 2),
                " ", base::units(time_elapsed))

  if(base::is.null(expression.matrix)){
    expression.matrix <- expression.levels
  }

  expression.levels <- base::unname(expression.levels)
  expression.levels.sorted <- expression.levels
  expression.levels.similarity.sorted <- expression.levels.similarity
  for(j in base::seq_len(base::ncol(expression.levels))){
    ordering <- base::order(expression.levels[,j])
    expression.levels.sorted[,j] <- expression.levels[ordering,j]
    expression.levels.similarity.sorted[,j] <- expression.levels.similarity[ordering,j]
  }

  if(noisyr::get_methods_correlation_distance(names=FALSE)[
    base::match(similarity.measure, noisyr::get_methods_correlation_distance())] == "d"){
    base::message("Chosen similarity metric ", similarity.measure, " is a dissimilarity, outputting inverse...")
    expression.levels.similarity <- 1/expression.levels.similarity
  }

  expression.summary <- base::list("expression.matrix" = expression.matrix,
                                   "expression.levels" = expression.levels.sorted,
                                   "expression.levels.similarity" = expression.levels.similarity.sorted)
  return(expression.summary)
}
