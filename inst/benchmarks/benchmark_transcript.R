# Compares the original (v1.1.0) and optimised transcript (BAM) approach on simulated data.
# Usage: Rscript benchmark_transcript.R <original R/ dir> <new R/ dir> <output dir for simulated BAMs>
args <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "make_test_bams.R"))
load_src <- function(dir){
  env <- new.env(parent = globalenv()); env$`%>%` <- magrittr::`%>%`; env$`%dopar%` <- foreach::`%dopar%`  # v1.1.0 used %dopar% without importing it
  for(f in list.files(dir, "\\.R$", full.names = TRUE)){
    eval(parse(text = gsub("noisyr::", "", readLines(f))), envir = env)
  }
  env
}
old <- load_src(args[1]); new <- load_src(args[2])
q <- function(e) suppressMessages(suppressWarnings(e))
time_it <- function(e) unname(system.time(q(e))["elapsed"])
results <- data.frame()
add <- function(step, t.old, t.new) {
  results <<- rbind(results, data.frame(step = step, original.s = t.old, optimised.s = t.new,
                                        speedup = round(t.old / t.new, 1)))
  print(tail(results, 1))
}

bams <- make_test_data(args[3], n.exons = 2000, n.samples = 4, reads.per.sample = 4e5)
genes <- q(new$cast_gtf_to_genes(file.path(args[3], "sim.gtf")))
add("similarity, 2000 exons x 4 BAMs",
    time_it(old$calculate_expression_similarity_transcript(bams, genes = genes)),
    time_it(new$calculate_expression_similarity_transcript(bams, genes = genes)))
# the parallel path of v1.1.0 needs the package installed (its PSOCK workers call
# noisyr::calculate_expression_profile), so compare installed versions for ncores > 1

# GTF parsing: a GTF with 200,000 exon lines
big.gtf <- tempfile(fileext = ".gtf")
lines <- readLines(file.path(args[3], "sim.gtf"))[-1]
writeLines(rep(lines, length.out = 2e5), big.gtf)
add("cast_gtf_to_genes, 200k lines",
    time_it(old$cast_gtf_to_genes(big.gtf)), time_it(new$cast_gtf_to_genes(big.gtf)))

# gene filtering: 20,000 genes with 10 exons each
n <- 2e5
g <- data.frame(gene_id = paste0("G", rep(1:2e4, each = 10)), seqid = "chr1", start = 1, end = 100)
m <- matrix(rexp(n * 4, 1/50), ncol = 4)
add("filter_genes_transcript, 20k genes / 200k exons",
    time_it(old$filter_genes_transcript(g, m, rep(30, 4))), time_it(new$filter_genes_transcript(g, m, rep(30, 4))))
print(results, row.names = FALSE)
