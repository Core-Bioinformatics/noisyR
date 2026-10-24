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
#' @param genes.per.batch number of exons whose reads are fetched from each BAM file in
#' a single query; larger batches mean fewer queries but more memory; default 500
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
  genes.per.batch=500,
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

  ## exons are processed in batches: the reads of a whole batch are fetched from each BAM
  ## file in one query, instead of one query per exon and BAM file
  batches <- base::split(base::seq_len(ngenes), base::ceiling(base::seq_len(ngenes) / genes.per.batch))
  process_batch <- function(rows){
    .expression_similarity_batch(genes.subset[rows, ], bams=bams, unique.only=unique.only,
                                 mapq.unique=mapq.unique, slack=slack, base.method=base.method)
  }

  base::message("Calculating expression similarity for ", ngenes, " genes...")
  start_time <- base::Sys.time()
  if(ncores>1){
    ncores <- base::min(ncores, parallel::detectCores())
    base::message("    ncores=", ncores)
    if(.Platform$OS.type == "unix"){
      by.batch <- parallel::mclapply(batches, process_batch, mc.cores=ncores)
    }else{
      cl <- parallel::makePSOCKcluster(ncores)
      on.exit(parallel::stopCluster(cl), add=TRUE)
      by.batch <- parallel::parLapply(cl, batches, process_batch)
    }
  }else{
    base::message("    ncores=1, running sequentially...")
    if(ngenes > 5000){
      base::message("    consider increasing ncores for more than a few thousand exons")
    }
    by.batch <- base::vector("list", base::length(batches))
    for(b in base::seq_along(batches)){
      by.batch[[b]] <- process_batch(batches[[b]])
      n <- base::max(batches[[b]])
      time_elapsed <- base::Sys.time() - start_time
      base::message("Done ", n, " genes out of ", ngenes, "; time elapsed: ",
                    base::round(time_elapsed, 2), " ", base::units(time_elapsed))
    }
  }
  failed <- base::vapply(by.batch, base::inherits, base::logical(1), what="try-error")
  if(base::any(failed)) base::stop("Similarity calculation failed: ", by.batch[[base::which(failed)[1]]])

  expression.levels <- base::do.call(base::rbind, base::lapply(by.batch, `[[`, "expression"))
  base::rownames(expression.levels) <- genes.subset$gene_id
  expression.levels.similarity <- base::do.call(base::rbind, base::lapply(by.batch, `[[`, "similarity"))
  if(save.image.every.1000){
    base::save.image()
  }
  time_elapsed <- base::Sys.time() - start_time
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

  expression.summary <- base::list("expression.matrix" = expression.matrix,
                                   "expression.levels" = expression.levels.sorted,
                                   "expression.levels.similarity" = expression.levels.similarity.sorted)
  return(expression.summary)
}

# Raw expression and average profile correlation with the other samples,
# for a batch of exons (rows of genes); one scanBam query per BAM file for the whole batch
.expression_similarity_batch <- function(genes, bams, unique.only, mapq.unique, slack, base.method){
  profiles <- .expression_profiles_batch(genes, bams, unique.only=unique.only,
                                         mapq.unique=mapq.unique, slack=slack)
  similarity <- base::t(base::vapply(profiles$profiles, function(profile){
    cors <- base::suppressWarnings(stats::cor(profile, method=base.method))
    cors[cors==1] = NA
    avCors <- base::colMeans(cors, na.rm=TRUE)
    avCors[base::is.nan(avCors)] = NA
    avCors
  }, base::numeric(base::length(bams))))
  base::list(expression=profiles$expression, similarity=similarity)
}

# Expression profiles (reads overlapping each position) and read counts of a batch of exons;
# same output as calling calculate_expression_profile() on each exon
.expression_profiles_batch <- function(genes, bams, unique.only=TRUE, mapq.unique=c(50,255), slack=200){
  lengths <- genes$end - genes$start + 1
  gr <- GenomicRanges::GRanges(seqnames = genes$seqid,
                               ranges = IRanges::IRanges(start=genes$start-slack,
                                                         end=genes$end))
  params <- Rsamtools::ScanBamParam(which=gr, what=c("pos", "qwidth", "mapq"))
  # scanBam returns one element per range, grouped by seqname (in the order of the
  # seqlevels) and in input order within each seqname
  result.order <- base::order(base::as.integer(GenomicRanges::seqnames(gr)), base::seq_along(gr))
  coord <- function(x) base::format(x, scientific=FALSE, trim=TRUE)
  expected.names <- base::paste0(genes$seqid, ":", coord(genes$start-slack), "-", coord(genes$end))[result.order]
  profiles <- base::lapply(lengths, function(len) base::matrix(0, nrow=len, ncol=base::length(bams)))
  expression <- base::matrix(0, nrow=base::nrow(genes), ncol=base::length(bams))
  for(j in base::seq_along(bams)){
    bamFile <- Rsamtools::BamFile(bams[j], base::paste0(bams[j], ".bai"))
    aln <- Rsamtools::scanBam(bamFile, param = params)
    if(!identical(base::names(aln), expected.names)){
      base::stop("Unexpected order of the reads returned by Rsamtools::scanBam")
    }
    for(k in base::seq_along(aln)){
      n <- result.order[k]
      reads <- aln[[k]]
      keep <- if(unique.only) reads$mapq %in% mapq.unique else base::rep(TRUE, base::length(reads$pos))
      expression[n, j] <- base::sum(keep)
      profiles[[n]][, j] <- .coverage_from_reads(reads$pos[keep] - genes$start[n] + 1,
                                                 reads$qwidth[keep], lengths[n])
    }
  }
  base::list(profiles=profiles, expression=expression)
}
