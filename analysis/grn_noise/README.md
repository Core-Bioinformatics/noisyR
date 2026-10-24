# noisyR-style noise thresholds for GRN adjacency matrices

Exploratory analysis and prototype. The question: can noisyR's idea carry over to
weighted adjacency matrices from GRN inference? The idea is that a signal is real
where replicate measurements agree, and noise sets in at the level below which
they stop agreeing. GRN examples are GENIE3/GRNBoost2 importances, ARACNe/CLR
mutual information, WGCNA adjacency/TOM, PIDC and correlation networks. The goal
is the edge-weight level below which edges cannot be told apart from noise.

Nothing in `R/` was changed. The prototype sources the package functions from this
checkout (`R/load_noisyr.R`).

```
analysis/grn_noise/
  R/load_noisyr.R       source the package R/ files without installing
  R/adjacency_noise.R   prototype: edge extraction, tail-score transform, noisyr_adjacency(),
                        concordance similarity, permutation / stability baselines, evaluation
  R/simulate.R          sparse true GRN, linear SEM expression, correlation + GENIE3-like
                        inference proxies, replicate schemes, direct noisy-weight simulation
  run_benchmark.R       accuracy (writes results/benchmark.csv, threshold_sweep.csv, pr_curves.csv)
  run_scaling.R         runtime / memory by stage (results/scaling.csv)
  plot_results.R        results/*.png
```

Run from the repository root with an R library that has the noisyR dependencies
(philentropy, preprocessCore, magrittr, ggplot2). The full benchmark at 500/2k/5k
genes takes a few minutes on 5 cores.

```bash
Rscript analysis/grn_noise/run_benchmark.R 500,2000,5000
```

---

## 1. Mapping noisyR onto adjacency matrices

| noisyR (counts) | adjacency version |
|---|---|
| sample (column) | one **replicate inference**: the same tool run on a disjoint subset of cells, or another tool on the same data |
| gene (row) | one **candidate edge**: the upper triangle for symmetric matrices; all off-diagonal regulator→target pairs for directed ones |
| expression level | the edge weight, on a common **tail-score** scale (below) |
| windowed similarity across samples | agreement of the same edges across replicates, in windows of edges sorted by weight in one replicate |
| noise threshold (expression) | edge-weight threshold, reported as a weight per replicate and as a number of edges |

**What counts as a replicate matters more than anything else.** noisyR assumes that
noise is independent between samples. For networks that only holds when the
replicates come from **disjoint sets of cells** (split-half or K-fold partitions,
biological replicates, batches). Bootstrap resamples, or different random seeds of
GENIE3/GRNBoost2, share the same cells. A spurious edge caused by finite-sample
correlation then shows up in every replicate, so it looks reproducible. Seeds
capture only the algorithm's randomness. In the benchmark, bootstrap replicates
pushed the concordance threshold to the bottom of the weight range. On the
correlation network that selected 62k–6.2M edges at a precision of 0.002–0.015. The
GENIE3-like proxy at 500 genes selected 12.5k edges (F1 0.14). The same data split
into 5 disjoint parts gave sensible thresholds.

Different tools on the same cells (GENIE3 vs CLR vs PIDC) are "replicates" of the
method, not of the data. Their agreement measures consensus between tools, not
noise. That is still useful (an ensemble threshold), but it is a different question.

**Rows: edges vs per-regulator/per-target rows.** Edges are the natural unit, and
the prototype uses them. The alternative is to treat each regulator's (or target's)
weight vector as a "gene", with the similarity being the correlation of that row
across replicates. That answers "which regulators have a reproducible target
profile", which is a regulator-level QC, not an edge threshold. It also costs
O(G) rows of length G instead of O(G²) rows of length R. It is worth adding later
as a diagnostic, not as the main method.

## 2. What carries over and what has to change

**Scale: use rank-based tail scores, not raw weights.** Edge weights are heavy
tailed and tool specific. GENIE3 importances sum to 1 per target, MI is in nats,
TOM is in [0,1], and |r| is bounded with most mass near 0. noisyR's boxplot and
loess methods bin on `log2(x + 1)`, which squashes weights of 0.001–0.05 into a
single bin. The prototype converts each replicate to a tail score
`x = n / n.ge − 1`, where `n.ge` is the number of edges at least as strong in that
replicate. Then `log2(x + 1) = −log2(tail fraction)`, which runs from 0 for the
weakest edge to log2(n) for the strongest, with one bin of 0.1 per 7% change in
tail rank. This is quantile normalisation onto a fixed, tool-agnostic reference.
It makes GENIE3, CLR and |r| directly comparable, and it handles ties (for
example the many zero importances in GRNBoost2) without special cases.
Thresholds are mapped back to a raw weight per replicate (`weight.threshold`).
In the direct-weight simulation, the mean of tail scores ranked edges much better
than the mean of raw heavy-tailed weights: AUPR 0.70 vs 0.26 at 500 genes.

**Similarity: noisyR's within-window correlation is fragile on edges.** The real
edges are a tiny top tail, about 0.1–1% of the O(G²) candidates. noisyR windows are
linear in rank, so almost every window sits in the noise. The top few windows mix
the whole signal region, and the correlation inside a narrow rank window is weak
because the range is restricted. Pearson, Spearman and philentropy distances all
behave the same way, because the problem is the windowing, not the metric. With
the best settings from a small grid (top 5·G edges, window 10%, Boxplot-IQR,
threshold 0.25), the native version did well at 500 genes. On the GENIE3-like
proxy it got F1 0.85 (oracle 0.90). At 2k and 5k genes it collapsed to a
high-precision core: recall 0.19 and 0.12, F1 0.31 and 0.21.

The prototype adds a **concordance** similarity that plugs into the same
`calculate_noise_threshold_curves()` / `noise_thresholds_from_curves()` machinery.
For each window of edges in replicate j's ordering, it is the fraction of (edge,
other replicate) pairs in which the edge's log tail score is within `delta = 1`
(one doubling of tail rank) of its score in j, or above it. A noise edge at tail
fraction f reproduces with probability about 2f. That is close to 0 everywhere
that matters, while reproducible edges approach 1. The threshold is
therefore readable as "at least 30% of replicate comparisons reproduce the edge".
It costs O(n) per replicate after one sort (cumulative sums instead of a
`cor()` per window).

**Threshold methods.** All 10 noisyR methods run unchanged on the curves. The
boxplot methods (default Boxplot-IQR) were the most stable, and loess-25 was
similar. The density-based methods do not apply, because the "first minimum of
the density" has no meaning on rank scores. In concordance space, the similarity
threshold behaves like a minimum local reproducibility. Thresholds of 0.2–0.4 gave
F1 within 0.05 of the best for the directed proxy (`results/threshold_sweep.png`).
0.3 is the prototype default. The split into curves and read-off lets all 19
thresholds be re-evaluated in 0.4 s at 5k genes (directed) or 1.5 s
(12.5M-edge symmetric).

**Symmetric vs asymmetric.** Symmetric matrices (correlation, MI, CLR, WGCNA, TOM)
use the upper triangle only, which halves the work and avoids counting each edge
twice. Directed matrices (GENIE3, GRNBoost2: regulators × targets) use all entries
except self-loops, matched by dimnames for rectangular TF × gene matrices.
If you symmetrise a GENIE3 matrix first (max or mean of the two directions), it
becomes a symmetric input.

**Sparsity.** Tools that prune (GRNBoost2 output, ARACNe after DPI, kNN graphs)
give edge lists with implicit zeros. Zeros all get tail score 0 (ties at the
bottom), so only the edge universe size `n` is needed, not dense storage. The
union of the replicates' edge lists is the row set.

**Direct vs indirect edges: a limit of any reproducibility criterion.** Indirect
associations (co-regulation, cascades) are real features of the data. They
reproduce across disjoint cell sets just as well as direct edges do. On the
correlation network, neither noisyR nor permutation can separate them. Measured
against direct edges, the best possible F1 of the |r| ranking is 0.59–0.61. The
reproducibility threshold separates "statistically real association" from noise.
Whether an association is direct has to be settled by the inference method
(partial correlation, CLR/DPI, regression-based importances), before thresholding.

## 3. Computational optimisations

| Optimisation | What the prototype does | Effect (G = 5000) |
|---|---|---|
| Upper triangle only | `edge_index()` for symmetric inputs | 12.5M instead of 25M rows |
| Sort each replicate once | one partial sort to find the top-M cut, then one radix sort for exact tail ranks; `findInterval` for lookups | 5 s for 5 × 12.5M edges |
| **Top-M pre-filter** | keep the union of each replicate's top M = 10·G edges, with exact tail scores computed against all n edges. **Adaptive**: if the threshold lands in the top half of the kept set, M is multiplied by 4 and the step is redone | correlation network: concordance 16 s instead of 149 s, peak heap 2.3 GB instead of 4.9 GB; native noisyR 5 s instead of 168 s |
| Window statistic via cumulative sums | concordance curves in O(n) per replicate | negligible next to the sort |
| Curve / read-off split | `edges_from_curves()` reuses the curves for any similarity threshold | 19 thresholds in 1.5 s |
| Parallelism over replicates | inference with `mclapply`; native similarity with the new `n.cores` | inference of 5 replicates in 3.4 s (corr) and 1.4 s (ridge) |
| One replicate in memory at a time | `prepare_edges` pulls one weight vector at a time | peak about one dense matrix plus the edge × replicate matrix of kept edges |

Not yet in the prototype, but needed for 20k × 20k (200M symmetric edges,
1.6 GB per replicate as doubles):

- **Streaming / chunking.** Compute or read the adjacency in column blocks. Keep a
  running top-M heap per replicate (exact), plus a 10⁴-point quantile sketch of
  everything else, which gives approximate tail scores for edges that are only in
  another replicate's top M. This needs O(M·R) memory instead of O(G²·R). The
  error is only in tail scores far below the threshold, which do not affect it.
- **Sparse/edge-list input.** Read arboreto/ARACNe outputs directly as `(i, j, w)`
  and never densify.
- **Float storage** (`float` package, or integer tail ranks) halves memory again.
- Rank the full-data network once and only use the replicates to set K (see
  section 5). The replicates then never need to be stored densely.

Measured runtimes are in `results/scaling.csv` and `results/scaling.png`. The
setup was 1000 cells, 5 replicates and 5 cores on an Apple-silicon Mac with
reference BLAS. Peak memory is the main R process's heap; forked workers are not
included.

| G | edges (corr / ridge) | inference ×5 | concordance threshold, top-M | without pre-filter | permutation null |
|---|---|---|---|---|---|
| 500 | 125k / 25k | 0.05 s | 0.96 s / 0.04 s | 0.23 s / 0.05 s | 0.33 s / 0.16 s |
| 2000 | 2.0M / 0.4M | 0.57 s / 0.20 s | 2.6 s / 0.28 s | 8.9 s / 0.87 s | 5.0 s / 1.6 s |
| 5000 | 12.5M / 2.5M | 3.4 s / 1.4 s | 16 s / 2.0 s | 149 s / 13 s | 31 s / 10 s |

(At 500 genes the correlation network is not sparse at the threshold, so the
adaptive pre-filter widened M three times. That is why it is slower than no
pre-filter there.)

## 4. Benchmark

**Simulation** (`R/simulate.R`): a sparse acyclic GRN. 10% of genes are TFs, each
gene has 1 + Poisson(1) TF regulators, and effects are ±U(0.1, 1). This gives about
2·G true edges. Expression for 1000 cells comes from a linear structural model,
standardised gene by gene. Three inference settings:

- **ridge**: a GENIE3-like directed TF × gene importance. Per target, a ridge
  regression on all TFs except itself, |standardised β| normalised to sum to 1
  per target. Truth is the direct edges.
- **corr**: the symmetric |Pearson r| network. Truth is either the direct edges
  or "associations", meaning |population r| ≥ 0.1 computed exactly from the model.
- **direct**: no expression. Each replicate's weights are signal on true edges
  (with log-normal multiplicative noise) plus independent heavy-tailed log-normal
  background noise on every edge.

Replicates: 5 disjoint splits of 200 cells each, or 5 bootstraps. Baselines:
permutation null (each gene permuted across cells, same inference, threshold at
estimated FDR ≤ 5% on the full-data network); stability selection (top-q per
replicate, kept if selected in ≥ 80% of replicates, q from the E[V] ≤ G/20 bound);
and the oracle best-F1 cut of the full-data and consensus rankings.

**Results.** Each cell is precision / recall / F1 (number of edges), for disjoint
splits. "+full-rank" means the noisyR edge count K applied to the full-data
ranking. Full table: `results/benchmark.csv`. Figures:
`results/pr_operating_points.png` and `results/threshold_sweep.png`.

GENIE3-like directed network, direct-edge truth:

| method | G = 500 | G = 2000 | G = 5000 |
|---|---|---|---|
| oracle (full-data ranking) | 0.90/0.90/**0.90** (974) | 0.89/0.87/**0.88** (3875) | 0.89/0.85/**0.87** (9567) |
| **noisyR-concordance + full-rank** | 0.92/0.87/**0.89** (924) | 0.93/0.73/**0.82** (3099) | 0.92/0.78/**0.84** (8471) |
| noisyR-concordance (consensus) | 0.91/0.86/0.88 | 0.92/0.72/0.80 | 0.84/0.72/0.78 |
| noisyR-native (Pearson in window) | 0.79/0.92/0.85 | 0.95/0.19/0.31 | 0.94/0.12/0.21 |
| permutation, FDR 5% | 0.92/0.86/**0.89** (917) | 0.83/0.91/**0.87** (4312) | 0.77/0.91/**0.83** (11732) |
| stability selection | 0.95/0.53/0.68 | 0.92/0.67/0.78 | 0.87/0.68/0.77 |

Correlation network:

| method | association truth, G = 2000 | direct truth, G = 2000 | association, G = 5000 |
|---|---|---|---|
| oracle (full-data ranking) | 0.82/0.75/0.79 (80k) | 0.73/0.50/0.59 (2.7k) | 0.81/0.72/0.76 (214k) |
| noisyR-concordance | 1.00/0.27/0.42 (23k) | 0.14/0.81/0.23 | 1.00/0.18/0.30 (42k) |
| permutation, FDR 5% | 0.61/0.88/0.72 (127k) | 0.03/0.96/0.06 | 0.61/0.84/0.71 (329k) |

What this shows:

1. **Directed, importance-type networks (the GENIE3/GRNBoost2 use case) are where
   it works.** Concordance-noisyR on disjoint splits, with K applied to the
   full-data ranking, reached F1 0.89 / 0.82 / 0.84 (oracle 0.90 / 0.88 / 0.87).
   Its precision held at 0.92–0.93 at every size. The permutation FDR cut had
   slightly better F1 at 2k, but its realised precision fell from 0.92 to 0.77 as
   G grew, against a 5% FDR target. Indirect TF–target signal in the data is not
   in a permutation null, so permutation FDRs are optimistic for GRNs. The noisyR
   threshold is conservative in the direction of precision. It also needs no
   extra inference on permuted data, which costs more than a whole replicate set
   for the expensive tools: GENIE3 on 20k genes runs for hours.
2. **The native noisyR within-window correlation does not scale to O(G²) edges.**
   It was fine at 500 genes but missed 80–90% of true edges at 2k–5k. It should
   not be shipped as-is for adjacency matrices.
3. **Correlation and other association networks.** A reproducibility threshold
   (noisyR or permutation) finds associations, not direct edges. Permutation was
   the better association detector here. Concordance-noisyR selected a
   precision-1.00 core of associations (recall 0.18–0.27). The 500-gene
   correlation case was unstable: the curve never dropped below 0.3, so all edges
   above the weakest third were kept. Thresholding a raw correlation network is
   not the right target. Threshold after CLR, DPI or partial correlation.
4. **Ranking with all cells beats ranking a consensus of splits.** For the
   directed proxy, AUPR was 0.90 vs 0.84 at 2k and 0.89 vs 0.80 at 5k. So:
   replicates set the threshold K, and the full-data network is what gets
   thresholded.
5. **The direct heavy-tailed weight simulation is hard by design.** The oracle F1
   was 0.71 / 0.42 / 0.18 for the consensus ranking. noisyR returned a
   high-precision core (precision 0.98 / 0.93 / 0.52 at 500 / 2k / 5k), with low
   recall. At 5k the signal is buried, and no method found much.

## 5. Alternatives worth comparing against

- **Permutation / null networks.** These are the default baseline and are well
  calibrated for "is there any association". They are optimistic for direct
  edges, as shown above. Each null costs one full inference. Column-wise
  permutation also breaks the cell-level structure (batches, cell types), unless
  it is done within strata.
- **Stability selection** (Meinshausen–Bühlmann; used by e.g. the Inferelator) has
  a finite-sample bound on E[V]. In this setting it was precise but conservative,
  with recall 0.53–0.68. It needs subsamples of about n/2 cells, not the n/5 used
  here.
- **IDR (irreproducible discovery rate, Li et al. 2011)** is the closest
  statistical relative of the concordance curve. It is a copula mixture model on
  replicate ranks, and it gives a local reproducibility rate per edge. It is a
  natural "model-based" version of the noisyR concordance curve and worth adding
  as a comparator. The `idr` CRAN package is O(n) but assumes 2 replicates.
- **Empirical-Bayes / mixture models on weights.** For example, a two-group local
  fdr (`locfdr`, `fdrtool`) on Fisher-z correlations or log importances. This needs
  a parametric null shape, which is easy for correlations but not for GENIE3
  importances. A two- or three-component mixture on the per-edge consensus tail
  score is cheap and could be added as a single-network fallback, when no
  replicates exist.
- **Top-k per target** (e.g. GENIE3 users keeping the top 100k edges, or the top 3
  regulators per target) is the de-facto practice. noisyR's K gives a
  data-driven replacement for it.

## 6. Summary and recommended next steps

- The noisyR principle (threshold where replicate agreement collapses) **transfers
  well to directed importance networks, provided replicates are disjoint cell
  subsets**. Use a rank-based tail-score scale and a concordance similarity. The
  existing `calculate_noise_threshold_curves()` / `noise_thresholds_from_curves()`
  work unchanged on the concordance curves. With K applied to the full-data
  network, F1 was 0.82–0.89, within 0.01–0.06 of the oracle. Precision held at
  about 0.92 up to 5k genes, where the permutation FDR cut drifted down to 0.77.
- The **native within-window correlation should not be reused** for edges, and
  **bootstrap or seed replicates must be rejected or warned about**.
- For association networks, document that the threshold separates associations
  from noise, not direct from indirect edges.
- **Recommendation: add `noisyr_adjacency()`**, in a reduced form of this
  prototype:
  - inputs are a list of matrices or edge lists, plus `symmetric`, `top.m`
    (adaptive) and `method.chosen`;
  - similarity is `"concordance"`, with the noisyR metrics kept as options;
  - it returns K, per-replicate weight thresholds and the curves (for plotting with
    the existing `plot_expression_similarity`), plus a helper that splits cells
    into K folds and calls a user-supplied inference function;
  - streaming/edge-list input is needed before it is advertised for 20k genes.
- **App tab: yes, but after the function.** The curve/read-off split already makes
  a similarity-threshold slider instant (1.5 s for 19 thresholds on 12.5M edges).
  The tab would show the concordance curve per replicate, K and the weight
  threshold, plus a download of the thresholded edge list.
- Before publishing: validate on a real dataset with a reference network (e.g.
  BEELINE scRNA-seq benchmarks with ChIP/STRING ground truth) using real
  GRNBoost2/GENIE3 runs on 5-fold cell splits. Add IDR and a locfdr mixture as
  comparators. Test a non-linear or count-based simulation (dropouts, NB noise),
  since the linear Gaussian SEM here favours the ridge proxy.
