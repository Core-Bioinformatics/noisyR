## Runtime / memory as the number of genes grows.
## Run from the repository root:
##   Rscript analysis/grn_noise/run_scaling.R [G values, default 500,2000,5000]
## Writes analysis/grn_noise/results/scaling.csv
## Memory is the peak of the R heap of the main process during each stage
## (gc "max used"); forked workers (mclapply) are not included.

source("analysis/grn_noise/R/load_noisyr.R")
source("analysis/grn_noise/R/adjacency_noise.R")
source("analysis/grn_noise/R/simulate.R")
ns <- load_noisyr(".")

args <- commandArgs(trailingOnly = TRUE)
G.values <- if(length(args)) as.integer(strsplit(args[1], ",")[[1]]) else c(500, 2000, 5000)
n.cells <- 1000
n.rep <- 5
n.cores <- 5

rows <- list()
stage <- function(G, inference, name, expr, n.edges = NA){
  invisible(gc(reset = TRUE))
  t <- system.time(val <- force(expr))[["elapsed"]]
  g <- gc()
  mem <- sum(g[, ncol(g)])   # "max used" in Mb
  rows[[length(rows) + 1]] <<- data.frame(G = G, inference = inference, stage = name,
                                          seconds = round(t, 2), peak.Mb = round(mem),
                                          n.edges = n.edges)
  message(sprintf("  %-6s %-40s %8.2f s %7.0f Mb", inference, name, t, mem))
  val
}

for(G in G.values){
  message("=== G = ", G)
  net <- simulate_network(G, seed = G)
  X <- simulate_expression(net, n.cells, seed = G)
  for(inf in c("corr", "ridge")){
    symmetric <- inf == "corr"
    f <- if(symmetric) infer_corr else function(x) infer_ridge(x, net$n.tf)
    tl <- truth_index(net, symmetric)
    reps <- stage(G, inf, "inference: 5 split replicates (5 cores)",
                  replicate_inferences(X, f, n.rep, "split", seed = G, n.cores = n.cores))
    n <- length(edge_index(reps[[1]], symmetric))
    edges <- stage(G, inf, "prepare_edges, top.m = 10 G", prepare_edges(reps, symmetric, 10 * G), n)
    res <- stage(G, inf, "concordance + threshold (top.m = 10 G)",
                 noisyr_adjacency(reps, ns, symmetric, top.m = 10 * G,
                                  similarity.measure = "concordance",
                                  similarity.threshold = 0.3), n)
    rows[[length(rows)]]$F1.direct <- eval_selection(res$selected, tl)$F1
    stage(G, inf, "concordance + threshold (no pre-filter)",
          noisyr_adjacency(reps, ns, symmetric, similarity.measure = "concordance",
                           similarity.threshold = 0.3), n)
    stage(G, inf, "native noisyR (top.m = 5 G, 5 cores)",
          noisyr_adjacency(reps, ns, symmetric, top.m = 5 * G, window.frac = 0.1,
                           n.cores = n.cores), n)
    stage(G, inf, "native noisyR (no pre-filter, 5 cores)",
          noisyr_adjacency(reps, ns, symmetric, window.frac = 0.1, n.cores = n.cores), n)
    curves <- res$curves
    stage(G, inf, "re-threshold 19 similarity thresholds",
          lapply(seq(0.05, 0.95, 0.05), function(st) edges_from_curves(curves, res$edges, st, ns)), n)
    stage(G, inf, "stability selection", stability_selection(edges, ev = G / 20), n)
    stage(G, inf, "permutation null (1 inference + FDR cut)", {
      full <- f(X)
      idx <- edge_index(full, symmetric)
      permutation_threshold(full[idx], f(permute_expression(X, seed = G))[idx])
    }, n)
    rm(reps, edges, res); invisible(gc())
  }
}

res <- do.call(rbind, lapply(rows, function(d){ if(is.null(d$F1.direct)) d$F1.direct <- NA; d }))
dir.create("analysis/grn_noise/results", showWarnings = FALSE)
utils::write.csv(res, "analysis/grn_noise/results/scaling.csv", row.names = FALSE)
options(width = 200)
print(res, row.names = FALSE)
