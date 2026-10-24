## Accuracy benchmark: how well do noisyR-style edge-weight thresholds recover the true
## network, compared with a permutation null, stability selection and the oracle cut?
## Run from the repository root:
##   Rscript analysis/grn_noise/run_benchmark.R [G values, default 500,2000]
## Writes analysis/grn_noise/results/{benchmark,threshold_sweep,pr_curves}.csv

source("analysis/grn_noise/R/load_noisyr.R")
source("analysis/grn_noise/R/adjacency_noise.R")
source("analysis/grn_noise/R/simulate.R")
ns <- load_noisyr(".")

args <- commandArgs(trailingOnly = TRUE)
G.values <- if(length(args)) as.integer(strsplit(args[1], ",")[[1]]) else c(500, 2000)
n.cells <- 1000
n.rep <- 5
n.cores <- 5
out.dir <- "analysis/grn_noise/results"

## the two noisyR variants compared throughout
SETTINGS <- list(
  native = list(similarity.measure = "correlation_pearson", top.m.per.gene = 5,
                window.frac = 0.1, similarity.threshold = 0.25),
  concordance = list(similarity.measure = "concordance", top.m.per.gene = 10,
                     window.frac = 0.02, similarity.threshold = 0.3)
)
SWEEP <- seq(0.05, 0.95, 0.05)

rows <- list(); sweep <- list(); prs <- list()

## evaluate a selection against every truth definition of the scenario
record <- function(base, method, sel, truths, extra = list()){
  for(tn in names(truths)){
    row <- cbind(base, data.frame(method = method, truth = tn),
                 eval_selection(sel, truths[[tn]]))
    if(length(extra)) row <- cbind(row, as.data.frame(extra))
    rows[[length(rows) + 1]] <<- row
  }
}

## full.rank: candidate edges ordered by the full-data inference, if there is one;
## the noisyR edge count K is then also applied to it (top K of the full ranking)
run_replicates <- function(base, reps, symmetric, truths, G, full.rank = NULL){
  for(nm in names(SETTINGS)){
    s <- SETTINGS[[nm]]
    res <- noisyr_adjacency(reps, ns, symmetric, top.m = s$top.m.per.gene * G,
                            window.frac = s$window.frac,
                            similarity.measure = s$similarity.measure,
                            similarity.threshold = s$similarity.threshold, n.cores = n.cores)
    record(base, paste0("noisyR-", nm), res$selected, truths,
           list(weight.threshold = mean(res$weight.threshold),
                n.ge.cv = stats::sd(res$n.ge) / mean(res$n.ge),
                top.m = res$edges$top.m))
    if(!is.null(full.rank)){
      record(base, paste0("noisyR-", nm, "+full-rank"), full.rank[seq_len(res$n.edges)], truths)
    }
    edges <- res$edges
    for(st in SWEEP){
      r <- suppressWarnings(edges_from_curves(res$curves, edges, st, ns))
      for(tn in names(truths)){
        sweep[[length(sweep) + 1]] <<- cbind(base, data.frame(
          method = paste0("noisyR-", nm), truth = tn, similarity.threshold = st),
          eval_selection(r$selected, truths[[tn]]))
      }
    }
    if(nm == "concordance"){
      ## consensus ranking quality (within the pre-filtered edges)
      score <- rowMeans(log2(edges$X + 1))
      for(tn in names(truths)){
        pr <- pr_curve(score, edges$idx[edges$kept] %in% truths[[tn]],
                       n.true = length(truths[[tn]]))
        record(base, "oracle-consensus", edges$idx[edges$kept[
          order(score, decreasing = TRUE)[seq_len(pr$best.k)]]], truths[tn],
          list(aupr = pr$aupr))
        cv <- pr$curve[unique(round(exp(seq(0, log(nrow(pr$curve)), length.out = 300)))), ]
        prs[[length(prs) + 1]] <<- cbind(base, data.frame(ranking = "consensus", truth = tn), cv)
      }
      ss <- stability_selection(edges, ev = G / 20, pi.thr = 0.8)
      record(base, "stability-selection", ss$selected, truths)
    }
  }
}

for(G in G.values){
  message("=== G = ", G)
  net <- simulate_network(G, seed = G)
  X <- simulate_expression(net, n.cells, seed = G)
  P <- population_correlation(net)
  X.perm <- permute_expression(X, seed = G)

  ## 1. correlation network (symmetric) and 2. GENIE3-like importance (directed)
  for(inf in c("corr", "ridge")){
    symmetric <- inf == "corr"
    f <- if(symmetric) infer_corr else function(x) infer_ridge(x, net$n.tf)
    full <- f(X)
    idx <- edge_index(full, symmetric)
    truths <- list(direct = truth_index(net, symmetric))
    if(symmetric){
      ## true associations: |population correlation| >= 0.1 (includes indirect edges)
      truths$association <- idx[abs(P[idx]) >= 0.1]
    }
    base <- data.frame(G = G, inference = inf, scheme = NA)
    w <- full[idx]
    for(tn in names(truths)){
      pr <- pr_curve(w, idx %in% truths[[tn]])
      o <- order(w, decreasing = TRUE)
      record(base, "oracle-full", idx[o[seq_len(pr$best.k)]], truths[tn], list(aupr = pr$aupr))
      cv <- pr$curve[unique(round(exp(seq(0, log(nrow(pr$curve)), length.out = 300)))), ]
      prs[[length(prs) + 1]] <- cbind(base, data.frame(ranking = "full", truth = tn), cv)
    }
    perm <- permutation_threshold(w, f(X.perm)[idx], fdr = 0.05)
    record(base, "permutation-FDR5%", idx[w >= perm$weight.threshold], truths,
           list(weight.threshold = perm$weight.threshold))

    for(scheme in c("split", "bootstrap")){
      message("  ", inf, " / ", scheme)
      reps <- replicate_inferences(X, f, n.rep, scheme, seed = G, n.cores = n.cores)
      run_replicates(data.frame(G = G, inference = inf, scheme = scheme),
                     reps, symmetric, truths, G, full.rank = idx[order(w, decreasing = TRUE)])
    }
  }

  ## 3. direct noisy-weight simulation: independent heavy-tailed noise per replicate
  message("  direct weights")
  reps <- simulate_direct_weights(net, n.rep, seed = G)
  truths <- list(direct = truth_index(net, TRUE))
  idx <- edge_index(reps[[1]], TRUE)
  w <- Reduce(`+`, reps)[idx] / n.rep
  null <- Reduce(`+`, simulate_direct_weights(net, n.rep, seed = G + 1, signal = 0))[idx] / n.rep
  base <- data.frame(G = G, inference = "direct", scheme = NA)
  pr <- pr_curve(w, idx %in% truths$direct)
  record(base, "oracle-full", idx[order(w, decreasing = TRUE)[seq_len(pr$best.k)]], truths,
         list(aupr = pr$aupr))
  perm <- permutation_threshold(w, null, fdr = 0.05)
  record(base, "permutation-FDR5%", idx[w >= perm$weight.threshold], truths)
  run_replicates(data.frame(G = G, inference = "direct", scheme = "independent"),
                 reps, TRUE, truths, G)
}

bind <- function(l){
  cols <- unique(unlist(lapply(l, names)))
  do.call(rbind, lapply(l, function(d){ d[setdiff(cols, names(d))] <- NA; d[cols] }))
}
dir.create(out.dir, showWarnings = FALSE)
res <- bind(rows)
utils::write.csv(res, file.path(out.dir, "benchmark.csv"), row.names = FALSE)
utils::write.csv(bind(sweep), file.path(out.dir, "threshold_sweep.csv"), row.names = FALSE)
utils::write.csv(bind(prs), file.path(out.dir, "pr_curves.csv"), row.names = FALSE)
num <- vapply(res, is.numeric, logical(1))
res[num] <- lapply(res[num], function(v) signif(v, 3))
options(width = 220)
print(res, row.names = FALSE)
