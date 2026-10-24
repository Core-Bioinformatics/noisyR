# Optimisation opportunities in noisyR

Status: **done** = implemented on this branch with output verified identical to v1.1.0;
**open** = identified, not implemented (reason given).

## Count matrix approach

| # | Where | Opportunity | Status |
|---|---|---|---|
| C1 | `calculate_expression_similarity_counts` | One `cor()`/`distance()` call per sample pair per window → one call per window comparing the sorted sample with all others (`cor(x, Y)`, `philentropy::dist_one_many`) | done |
| C2 | same | Samples are independent → `n.cores` (forked `mclapply`) | done |
| C3 | same | Pearson over all windows at once via prefix sums or banded matrix products (~100x) | open: one-pass sums lose precision on non-integer data; a C++ (Rcpp) two-pass kernel would be exact but adds compiled code |
| C4 | same | Spearman re-ranks every window (dominates its runtime) | open: ranks change per window; only a C++ kernel helps |
| C5 | `calculate_noise_threshold` | Loess fits, boxplot bins and densities do not depend on the similarity threshold → `calculate_noise_threshold_curves()` + `noise_thresholds_from_curves()` | done |
| C6 | `calculate_noise_threshold_method_statistics` | Recomputed every method for every threshold → curves once per method (56x) | done |
| C7 | same | Growing a tibble with `rbind` in a loop → `bind_rows` once | done |
| C8 | `optimise_window_length` | Symmetric t-test matrix computed twice → upper triangle only | done |
| C9 | same | Each window length recomputes the full similarity → can run window lengths in parallel | open: already benefits from C1/C2; nested parallelism not worth the memory |
| C10 | `cast_matrix_to_numeric` | `apply()` column by column → one `as.numeric()` on the whole matrix (also fixes single-row input) | done |
| C11 | `calculate_first_minimum_density` | `adjust` argument ignored; failure flag never reset across columns | done (bug fixes) |

## Transcript (point-by-point / BAM) approach

| # | Where | Opportunity | Status |
|---|---|---|---|
| T1 | `calculate_expression_profile` | Reads every BAM field (`scanBamWhat()`: sequences, qualities, tags) but uses only `pos`, `qwidth`, `mapq` | done |
| T2 | same | Coverage built with an R loop over every read → difference array + `cumsum` | done |
| T3 | `calculate_expression_similarity_transcript` | One `scanBam` query per exon per BAM file (open, seek, parse) → one query per batch of exons (`genes.per.batch`, default 500) per BAM file | done |
| T4 | same | Parallel `foreach(.combine = rbind)` over exons copies the growing result matrix at every step (quadratic) → parallelise over batches and bind once | done |
| T4b | same | `ncores > 1` crashed in v1.1.0 (`could not find function "%dopar%"`): `%dopar%` was never imported in NAMESPACE, so it only worked after `library(foreach)` | done (bug fix: no foreach needed) |
| T5 | same | PSOCK cluster on all platforms (workers re-load packages, data serialised) → forked `mclapply` on unix, PSOCK only on Windows | done |
| T6 | same | `save.image()` every 1000 exons saves the global environment, which never contains the partial results | open: kept as a final save for backwards compatibility; suggest deprecating |
| T7 | same | `cors[cors == 1] <- NA` is meant to drop the diagonal but also drops genuine correlations of exactly 1 between different samples (common for low-coverage exons with identical profiles) | open: would change results; recommend `diag(cors) <- NA` in a future version |
| T8 | same | Dead code: dissimilarity inversion after `stop()` for distances (and applied to the unsorted matrix) | done (removed) |
| T9 | `cast_gtf_to_genes` | Per-row `strsplit` of attributes via `lapply` → split all rows at once and index the first match | done |
| T10 | same | `read.table` type guessing on a large GTF → explicit `colClasses` | done |
| T11 | same | `grep(att_of_interest, ...)` matches substrings (e.g. `gene_id` inside `havana_gene_id`), and when the attribute is the last one on the line its value keeps the trailing `;` (`B;`) | open: kept for identical results; an anchored regex would fix both |
| T12 | `filter_genes_transcript` | `sapply` over unique genes comparing against all exons (O(genes x exons)) → one `%in%` | done |
| T13 | `remove_noise_from_bams` | Overlapping exon ranges queried repeatedly → `GenomicRanges::reduce` first (output unchanged; `filterBam` does not duplicate reads) | done |
| T14 | same | BAM files filtered one after another → could run in parallel (`mclapply` over files) | open: I/O bound; benefit depends on storage |
| T15 | `remove_noise_from_bams` vs profiles | `mapqFilter = mapq.unique` keeps mapq >= value, while the profiles keep mapq %in% values → inconsistent definitions of a unique read | open: behaviour change; worth harmonising |
| T16 | whole approach | Exons sharing a region are scanned once per exon | open: could scan reduced regions once and slice coverage per exon |

## Benchmarks

See `inst/benchmarks/benchmark_counts.R` and `inst/benchmarks/benchmark_transcript.R`
(the latter simulates a GTF and indexed BAMs with `make_test_bams.R`).
