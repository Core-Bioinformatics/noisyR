simulate_counts <- function(n.genes = 2000, n.samples = 5, seed = 1){
  set.seed(seed)
  mu <- stats::rlnorm(n.genes, meanlog = 3, sdlog = 2)
  sapply(seq_len(n.samples), function(j) stats::rnbinom(n.genes, mu = mu, size = 2))
}

# the pairwise loop of noisyr <= 1.1.0, kept as the reference implementation
reference_similarity <- function(expression.matrix, n.elements.per.window, n.step, fun){
  nsequence <- seq(1, nrow(expression.matrix) - n.elements.per.window + 1, n.step)
  out <- matrix(0, length(nsequence), ncol(expression.matrix))
  for(j in seq_len(ncol(expression.matrix))){
    sorted <- expression.matrix[order(expression.matrix[, j]), ]
    for(idx in seq_along(nsequence)){
      rows <- nsequence[idx]:(nsequence[idx] + n.elements.per.window - 1)
      sim <- numeric(ncol(expression.matrix))
      for(k in seq_len(ncol(expression.matrix))) if(j != k){
        sim[k] <- suppressWarnings(fun(cbind(sorted[rows, j], sorted[rows, k])))
      }
      out[idx, j] <- mean(sim[sim != 0])
    }
  }
  out
}

test_that("vectorised correlations match the pairwise loop", {
  m <- simulate_counts()
  for(method in c("pearson", "spearman")){
    es <- suppressMessages(calculate_expression_similarity_counts(
      m, similarity.measure = paste0("correlation_", method), n.elements.per.window = 200))
    ref <- reference_similarity(m, 200, 10, function(x) stats::cor(x, method = method)[1, 2])
    expect_equal(es$expression.levels.similarity, ref)
  }
})

test_that("vectorised distances match the pairwise loop", {
  m <- simulate_counts()
  es <- suppressMessages(calculate_expression_similarity_counts(
    m, similarity.measure = "distance_euclidean", n.elements.per.window = 200))
  ref <- reference_similarity(m, 200, 10, function(x){
    s <- colSums(x)
    if(any(is.na(s) | s == 0)) return(NA)
    suppressMessages(philentropy::distance(t(x) / s, method = "euclidean", test.na = FALSE))
  })
  expect_equal(es$expression.levels.similarity, 1 / ref)
})

test_that("parallel and serial similarity agree", {
  skip_on_os("windows")
  m <- simulate_counts()
  serial <- suppressMessages(calculate_expression_similarity_counts(m, n.elements.per.window = 200))
  forked <- suppressMessages(calculate_expression_similarity_counts(m, n.elements.per.window = 200, n.cores = 2))
  expect_identical(serial, forked)
})

test_that("curves give the same thresholds as calculate_noise_threshold for every method", {
  m <- simulate_counts()
  es <- suppressMessages(calculate_expression_similarity_counts(m, n.elements.per.window = 200))
  for(method in get_methods_calculate_noise_threshold()){
    curves <- suppressMessages(calculate_noise_threshold_curves(es, method.chosen = method))
    for(thr in c(0.1, 0.25, 0.5)){
      expect_equal(suppressMessages(noise_thresholds_from_curves(curves, thr)),
                   suppressMessages(calculate_noise_threshold(es, thr, method)))
    }
  }
})

test_that("method statistics keep the threshold-major row order", {
  m <- simulate_counts()
  es <- suppressMessages(calculate_expression_similarity_counts(m, n.elements.per.window = 200))
  methods <- c("Boxplot-IQR", "Line_plot-No_smoothing")
  st <- suppressMessages(calculate_noise_threshold_method_statistics(es, c(0.2, 0.3), methods))
  expect_equal(st$similarity.threshold, c(0.2, 0.2, 0.3, 0.3))
  expect_equal(paste(st$approach, st$method, sep = "-"), rep(methods, 2))
})

test_that("coverage from reads matches incrementing each read", {
  set.seed(2)
  len <- 300
  pos <- sample(-100:350, 500, TRUE)
  qwidth <- sample(0:100, 500, TRUE)
  ref <- numeric(len)
  for(r in seq_along(pos)){
    add <- seq(from = pos[r], length.out = qwidth[r])
    add <- add[add > 0 & add <= len]
    ref[add] <- ref[add] + 1
  }
  expect_equal(noisyr:::.coverage_from_reads(pos, qwidth, len), ref)
})

test_that("casting keeps dimnames and works for a single row", {
  x <- matrix(c("1", "2", "a"), 1, dimnames = list("g1", c("A", "B", "C")))
  out <- suppressWarnings(cast_matrix_to_numeric(x))
  expect_equal(dim(out), c(1, 3))
  expect_equal(dimnames(out), dimnames(x))
  expect_equal(as.vector(out), c(1, 2, 0))
})
