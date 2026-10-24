## Prototype: noisyR-style noise thresholds for weighted adjacency matrices.
##
## Replicate inferences (disjoint cell subsets, tools, ...) play the role of samples,
## edges play the role of genes. Each replicate's weights are put on a tool-agnostic
## "tail score" scale, x = n / n.ge - 1, where n.ge is the number of edges in that
## replicate at least as strong (so log2(x + 1) = -log2(tail fraction), 0 for the
## weakest edge, log2(n) for the strongest). The noisyR counts machinery then finds
## the tail score above which edges are reproducible across replicates.
##
## Requires `ns`, the environment returned by load_noisyr().

## Linear indices of the candidate edges of an adjacency matrix:
## symmetric -> upper triangle; directed -> all entries except self-loops
## (matched by dimnames for rectangular regulator x target matrices, else the diagonal).
edge_index <- function(adj, symmetric){
  if(symmetric){
    return(which(upper.tri(adj)))
  }
  rn <- rownames(adj); cn <- colnames(adj)
  if(!is.null(rn) && !is.null(cn)){
    self <- cbind(match(intersect(rn, cn), rn), match(intersect(rn, cn), cn))
  }else{
    d <- seq_len(min(dim(adj)))
    self <- cbind(d, d)
  }
  self.lin <- (self[, 2] - 1) * nrow(adj) + self[, 1]
  setdiff(seq_len(length(adj)), self.lin)
}

is_symmetric_adj <- function(adj){
  nrow(adj) == ncol(adj) && isTRUE(all.equal(adj[1:min(50, nrow(adj)), ],
                                             t(adj)[1:min(50, nrow(adj)), ]))
}

## Turn a list of replicate adjacency matrices into an edge x replicate matrix of
## tail scores, keeping only the union of each replicate's top.m edges.
## Each replicate is sorted once; only one full weight vector is held at a time.
prepare_edges <- function(adj.list, symmetric = NULL, top.m = NULL){
  if(is.null(symmetric)) symmetric <- is_symmetric_adj(adj.list[[1]])
  idx <- edge_index(adj.list[[1]], symmetric)
  n <- length(idx)
  n.rep <- length(adj.list)
  top.m <- if(is.null(top.m)) n else min(top.m, n)

  ## pass 1: partial sort to find each replicate's top.m cut
  keep <- logical(n)
  for(r in seq_len(n.rep)){
    w <- adj.list[[r]][idx]
    cut <- sort.int(w, partial = n - top.m + 1)[n - top.m + 1]
    keep <- keep | (w >= cut)
  }
  kept <- which(keep)

  ## pass 2: exact tail scores of the kept edges, one full sort per replicate
  X <- matrix(0, length(kept), n.rep)
  top.weights <- matrix(NA_real_, top.m, n.rep)
  for(r in seq_len(n.rep)){
    w <- adj.list[[r]][idx]
    s <- sort.int(w, method = "radix")
    n.ge <- n - findInterval(w[kept], s, left.open = TRUE)
    X[, r] <- n / n.ge - 1
    top.weights[, r] <- s[n:(n - top.m + 1)]
  }
  list(X = X, kept = kept, idx = idx, n = n, n.rep = n.rep, top.m = top.m,
       top.weights = top.weights, symmetric = symmetric, dims = dim(adj.list[[1]]))
}

## Concordance similarity, an alternative to the within-window correlation: for each
## window in replicate j's ordering, the fraction of (edge, other replicate) pairs in
## which the edge's log2 tail score is within `delta` of (or above) its score in j.
## Noise edges at tail fraction f reproduce with probability ~ 2^delta * f, so the curve
## is ~0 in the noise region and approaches 1 for reproducible edges.
## Output has the same layout as calculate_expression_similarity_counts().
## O(n.edges) per replicate after one sort (cumulative sums instead of windows).
calculate_edge_concordance <- function(X, n.elements.per.window, n.step.fraction = 0.05,
                                       delta = 1){
  L <- log2(X + 1)
  n <- nrow(L)
  w <- n.elements.per.window
  starts <- seq(1, n - w + 1, max(1, floor(w * n.step.fraction)))
  levels <- similarity <- matrix(NA_real_, length(starts), ncol(L))
  for(j in seq_len(ncol(L))){
    o <- order(L[, j])
    hit <- rowMeans(L[o, -j, drop = FALSE] >= L[o, j] - delta)
    cs <- c(0, cumsum(hit))
    cx <- c(0, cumsum(X[o, j]))
    similarity[, j] <- (cs[starts + w] - cs[starts]) / w
    levels[, j] <- (cx[starts + w] - cx[starts]) / w
  }
  list(expression.matrix = X, expression.levels = levels,
       expression.levels.similarity = similarity)
}

## noisyR on an edge x replicate tail-score matrix.
## Returns the per-replicate threshold (as edge count and as raw weight), the consensus
## number of edges above noise and the selected edges (linear indices into the matrix).
## similarity.measure is any noisyR correlation/distance or "concordance".
noisyr_adjacency <- function(
  adj.list,
  ns,
  symmetric = NULL,
  top.m = NULL,
  window.frac = 0.02,
  n.step.fraction = 0.05,
  similarity.measure = "correlation_pearson",
  method.chosen = "Boxplot-IQR",
  similarity.threshold = 0.25,
  binsize = 0.1,
  delta = 1,
  n.cores = 1
){
  repeat{
    edges <- prepare_edges(adj.list, symmetric, top.m)
    X <- edges$X
    n.window <- max(20, round(window.frac * nrow(X)))
    if(similarity.measure == "concordance"){
      summ <- calculate_edge_concordance(X, n.window, n.step.fraction, delta)
    }else{
      summ <- suppressMessages(ns$calculate_expression_similarity_counts(
        X, similarity.measure = similarity.measure, n.elements.per.window = n.window,
        n.step.fraction = n.step.fraction, n.cores = n.cores))
    }
    curves <- suppressMessages(ns$calculate_noise_threshold_curves(
      summ, method.chosen = method.chosen, binsize = binsize))
    out <- suppressWarnings(edges_from_curves(curves, edges, similarity.threshold, ns))
    ## the pre-filter must leave a stretch of noise below the threshold:
    ## if not, widen it and start again
    if(edges$top.m >= edges$n || max(out$n.ge) <= edges$top.m / 2) break
    top.m <- 4 * edges$top.m
    message("    threshold near the top.m boundary, retrying with top.m = ", top.m)
  }
  c(out, list(summary = summ, curves = curves, edges = edges))
}

## Read the edge set off precomputed curves (cheap: call it for many thresholds)
edges_from_curves <- function(curves, edges, similarity.threshold, ns){
  X <- edges$X
  thr <- suppressMessages(ns$noise_thresholds_from_curves(curves, similarity.threshold))
  n.ge <- pmax(1, round(edges$n / (thr + 1)))
  if(edges$top.m < edges$n && any(n.ge > 0.9 * edges$top.m)){
    warning("noise threshold close to the top.m pre-filter boundary; increase top.m")
  }
  weight.thr <- vapply(seq_len(edges$n.rep), function(r){
    edges$top.weights[min(n.ge[r], edges$top.m), r]
  }, numeric(1))

  ## consensus: mean log tail score; keep the median number of edges above noise
  k <- min(round(stats::median(n.ge)), nrow(X))
  score <- rowMeans(log2(X + 1))   # (could be cached across thresholds)
  sel <- edges$kept[order(score, decreasing = TRUE)[seq_len(k)]]
  list(n.edges = k, n.ge = n.ge, tail.threshold = thr, weight.threshold = weight.thr,
       selected = edges$idx[sel])
}

## Permutation-null baseline: smallest weight cut with estimated FDR <= fdr, where
## FDR(t) = #null edges >= t / #observed edges >= t (monotonised, as for q-values).
permutation_threshold <- function(w.obs, w.null, fdr = 0.05){
  o <- sort.int(w.obs, decreasing = TRUE, method = "radix")
  s.null <- sort.int(w.null, method = "radix")
  n.null.ge <- length(s.null) - findInterval(o, s.null, left.open = TRUE)
  est <- n.null.ge / seq_along(o) * length(w.obs) / length(w.null)
  q <- rev(cummin(rev(est)))
  k <- sum(q <= fdr)
  list(n.edges = k, weight.threshold = if(k > 0) o[k] else Inf)
}

## Stability selection (Meinshausen & Buhlmann 2010) on the replicates: each replicate
## selects its top q edges, edges selected in >= pi.thr of replicates are kept.
## q is set from the bound E[V] <= q^2 / ((2 pi.thr - 1) n).
stability_selection <- function(edges, ev = 50, pi.thr = 0.8){
  q <- floor(sqrt(ev * (2 * pi.thr - 1) * edges$n))
  q <- min(q, edges$top.m)
  freq <- rowMeans(edges$X >= edges$n / q - 1)
  sel <- edges$kept[freq >= pi.thr]
  list(n.edges = length(sel), q = q, selected = edges$idx[sel])
}

## Precision/recall/F1 of a selected edge set
eval_selection <- function(selected, truth.lin){
  tp <- sum(selected %in% truth.lin)
  k <- length(selected)
  p <- if(k > 0) tp / k else NA
  rcl <- tp / length(truth.lin)
  data.frame(n.edges = k, TP = tp, precision = p, recall = rcl,
             F1 = if(k > 0 && tp > 0) 2 * p * rcl / (p + rcl) else 0)
}

## Precision-recall curve of an edge ranking; n.true is the total number of true
## edges, if some are missing from `score` (e.g. outside the top.m pre-filter)
pr_curve <- function(score, is.true, n.true = sum(is.true)){
  o <- order(score, decreasing = TRUE)
  tp <- cumsum(is.true[o])
  k <- seq_along(o)
  prec <- tp / k
  rec <- tp / n.true
  f1 <- ifelse(tp > 0, 2 * prec * rec / (prec + rec), 0)
  list(aupr = sum(prec[is.true[o]]) / n.true,   # average precision
       best.k = which.max(f1), best.f1 = max(f1),
       curve = data.frame(k = k, precision = prec, recall = rec))
}
