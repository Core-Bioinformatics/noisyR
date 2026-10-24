# noisyr 1.2.0

## Performance (same results as 1.1.0)

* `calculate_expression_similarity_counts()` compares each window with all other samples in one
  call (`stats::cor(x, Y)` / `philentropy::dist_one_many()`) instead of one call per pair of samples.
  The new `n.cores` argument processes samples in parallel (unix-alikes).
* `calculate_noise_threshold()` is split into `calculate_noise_threshold_curves()` (the costly part,
  independent of the similarity threshold) and `noise_thresholds_from_curves()` (cheap). The
  duplicated loess code is shared, and boxplot bin quantiles are computed with `split()`.
* `calculate_noise_threshold_method_statistics()` computes the curves once per method, not once
  per (method, threshold) pair.
* `optimise_window_length()` gains `n.cores` and runs only half of the symmetric t-tests.
* `calculate_expression_profile()` builds read coverage with a difference array instead of a loop over reads.
* `cast_matrix_to_numeric()` converts the whole matrix at once and now also works for single-row matrices.
* Transcript approach: `calculate_expression_similarity_transcript()` fetches the reads of a batch of
  exons (new `genes.per.batch`, default 500) from each BAM file in one query, reads only the
  position, width and mapping quality fields, and parallelises over batches (forked on unix)
  instead of combining per-exon results with `rbind`.
* `cast_gtf_to_genes()` parses the attributes of all rows at once; `filter_genes_transcript()` no
  longer compares every gene with every exon; `remove_noise_from_bams()` merges overlapping exons
  before filtering.

## Bug fixes

* `calculate_first_minimum_density()` ignored its `adjust` argument (always used 2), and once one
  column had no clear minimum, every later column was skipped as well.

* `library(noisyr)` failed when utils was not attached (`.onAttach` called `citation()` unqualified).
* doParallel and foreach are no longer needed; `slides/` is excluded from the built package (was 10 MB).
* `calculate_expression_similarity_transcript(ncores > 1)` failed with
  `could not find function "%dopar%"` unless foreach was attached.

## New features

* `run_noisyr_app()`: interactive Shiny app for the count matrix approach.
* testthat tests that check the optimised code against the original implementations (counts and transcript approaches).
* `OPTIMISATION.md` lists the optimisation opportunities found, implemented or still open.
