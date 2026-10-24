# Compares the original (v1.1.0, master) and optimised counts pipeline.
# Usage: Rscript benchmark_counts.R <dir with original R/ sources> <dir with new R/ sources>
args <- commandArgs(trailingOnly = TRUE)
load_src <- function(dir){
  env <- new.env(parent = globalenv())
  for(f in list.files(dir, "\\.R$", full.names = TRUE)){
    eval(parse(text = gsub("noisyr::", "", readLines(f))), envir = env)
  }
  env
}
sim_counts <- function(n.genes, n.samples, seed = 1){
  set.seed(seed)
  mu <- rlnorm(n.genes, meanlog = 3, sdlog = 2)
  sapply(seq_len(n.samples), function(j) rnbinom(n.genes, mu = mu * runif(1, .7, 1.3), size = 2))
}
old <- load_src(args[1]); new <- load_src(args[2])
q <- function(e) suppressMessages(suppressWarnings(e))
time_it <- function(e) unname(system.time(q(e))["elapsed"])

m <- sim_counts(20000, 12)
thrs <- seq(0.1, 0.5, 0.05)
results <- data.frame()
add <- function(step, t.old, t.new) {
  results <<- rbind(results, data.frame(step = step, original.s = t.old, optimised.s = t.new,
                                        speedup = round(t.old / t.new, 1)))
  print(tail(results, 1))
}
for(meas in c("correlation_pearson", "correlation_spearman", "distance_euclidean")){
  add(paste("similarity", meas),
      time_it(old$calculate_expression_similarity_counts(m, meas)),
      time_it(new$calculate_expression_similarity_counts(m, meas)))
}
add("similarity correlation_pearson, 4 cores",
    time_it(old$calculate_expression_similarity_counts(m)),
    time_it(new$calculate_expression_similarity_counts(m, n.cores = 4)))
es <- q(new$calculate_expression_similarity_counts(m))
add("method statistics (9 thresholds x 10 methods)",
    time_it(old$calculate_noise_threshold_method_statistics(es, thrs)),
    time_it(new$calculate_noise_threshold_method_statistics(es, thrs)))
m.small <- sim_counts(5000, 6)
add("optimise_window_length (5k genes, 33 windows)",
    time_it(old$optimise_window_length(m.small, iteration.number = 20)),
    time_it(new$optimise_window_length(m.small, iteration.number = 20)))
pdf(NULL)
print(results, row.names = FALSE)
