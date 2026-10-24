## Simulated benchmark: a sparse true GRN, expression from a linear structural model,
## and two fast inference proxies (a correlation network and a GENIE3-like
## regulator -> target importance matrix).

## Sparse true network: the first n.tf genes are regulators (TFs). Every gene gets
## 1 + Poisson(mean.indegree - 1) regulators among the TFs that precede it
## (TFs are ordered, so the network is acyclic). Effect sizes are signed U(0.1, 1).
simulate_network <- function(G, frac.tf = 0.1, mean.indegree = 2, seed = 1){
  set.seed(seed)
  n.tf <- max(5, round(frac.tf * G))
  edges <- do.call(rbind, lapply(2:G, function(j){
    cand <- seq_len(min(j - 1, n.tf))
    k <- min(length(cand), 1 + stats::rpois(1, mean.indegree - 1))
    from <- cand[sample.int(length(cand), k)]
    data.frame(from = from, to = rep(j, k))
  }))
  edges$beta <- sample(c(-1, 1), nrow(edges), replace = TRUE) *
    stats::runif(nrow(edges), 0.1, 1)
  list(G = G, n.tf = n.tf, edges = edges)
}

## Expression (cells x genes) from x_j = sum_i beta_ij x_i + e_j, generated in
## topological order and standardised gene by gene.
simulate_expression <- function(net, n.cells = 1000, seed = 1){
  set.seed(seed)
  X <- matrix(stats::rnorm(n.cells * net$G), n.cells, net$G)
  by.target <- split(net$edges, net$edges$to)
  for(j in as.integer(names(by.target))){
    e <- by.target[[as.character(j)]]
    X[, j] <- X[, j] + X[, e$from, drop = FALSE] %*% e$beta
    X[, j] <- (X[, j] - mean(X[, j])) / stats::sd(X[, j])
  }
  colnames(X) <- paste0("g", seq_len(net$G))
  X
}

## Population correlation matrix implied by simulate_expression() (exact recursion
## over the topological order), used to define "true associations" for correlation
## networks, which include indirect (co-regulation, cascade) associations.
population_correlation <- function(net){
  C <- diag(net$G)
  by.target <- split(net$edges, net$edges$to)
  for(j in as.integer(names(by.target))){
    e <- by.target[[as.character(j)]]
    b <- e$beta
    s <- sqrt(drop(t(b) %*% C[e$from, e$from, drop = FALSE] %*% b) + 1)
    cj <- drop(b %*% C[e$from, seq_len(j - 1), drop = FALSE]) / s
    C[j, seq_len(j - 1)] <- cj
    C[seq_len(j - 1), j] <- cj
  }
  C
}

## |Pearson correlation| network, symmetric G x G
infer_corr <- function(X){
  abs(stats::cor(X))
}

## GENIE3-like importance proxy, directed n.tf x G: per target, ridge regression on
## all TFs except itself; importance = |standardised coefficient|, normalised to sum
## to 1 per target (as GENIE3 importances are). TF targets use the precision-matrix
## identity beta_{-j} = -Omega_{-j,j} / Omega_jj, so one T x T inverse serves all targets.
infer_ridge <- function(X, n.tf, lambda = 0.2){
  Z <- scale(X)
  S.TA <- crossprod(Z[, seq_len(n.tf)], Z) / (nrow(Z) - 1)
  Omega <- solve(S.TA[, seq_len(n.tf)] + diag(lambda, n.tf))
  B <- Omega %*% S.TA
  B[, seq_len(n.tf)] <- -sweep(Omega, 2, diag(Omega), "/")
  diag(B[, seq_len(n.tf)]) <- 0
  W <- abs(B)
  W <- sweep(W, 2, pmax(colSums(W), .Machine$double.eps), "/")
  dimnames(W) <- list(colnames(X)[seq_len(n.tf)], colnames(X))
  W
}

## Replicate inferences: disjoint cell subsets (independent noise) or bootstraps
## (resamples of the same cells, so spurious edges are shared between replicates)
replicate_inferences <- function(X, infer, n.rep = 5, scheme = c("split", "bootstrap"),
                                 seed = 1, n.cores = 1){
  scheme <- match.arg(scheme)
  set.seed(seed)
  cells <- if(scheme == "split"){
    split(sample.int(nrow(X)), rep_len(seq_len(n.rep), nrow(X)))
  }else{
    lapply(seq_len(n.rep), function(r) sample.int(nrow(X), replace = TRUE))
  }
  parallel::mclapply(cells, function(cc) infer(X[cc, , drop = FALSE]), mc.cores = n.cores)
}

## Permutation null: every gene permuted independently across cells
permute_expression <- function(X, seed = 1){
  set.seed(seed)
  apply(X, 2, sample)
}

## Linear indices of the true edges in the inferred matrix
truth_index <- function(net, symmetric){
  e <- net$edges
  if(symmetric){
    a <- pmin(e$from, e$to); b <- pmax(e$from, e$to)
    unique((b - 1) * net$G + a)
  }else{
    (e$to - 1) * net$n.tf + e$from
  }
}

## Direct weight simulation (no expression): weights of true edges are reproducible
## (true strength x lognormal multiplicative noise), all edges also carry independent
## heavy-tailed lognormal background noise; symmetric G x G.
simulate_direct_weights <- function(net, n.rep = 5, seed = 1, signal = 0.2,
                                    noise.meanlog = -4, noise.sdlog = 1.2){
  set.seed(seed)
  G <- net$G
  tl <- truth_index(net, symmetric = TRUE)
  strength <- abs(net$edges$beta)[!duplicated(
    (pmax(net$edges$from, net$edges$to) - 1) * G + pmin(net$edges$from, net$edges$to))]
  lapply(seq_len(n.rep), function(r){
    W <- matrix(stats::rlnorm(G * G, noise.meanlog, noise.sdlog), G, G)
    W[tl] <- W[tl] + signal * strength * stats::rlnorm(length(tl), 0, 0.5)
    W[lower.tri(W)] <- t(W)[lower.tri(W)]
    W
  })
}
