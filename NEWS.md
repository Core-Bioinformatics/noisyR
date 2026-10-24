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

## Bug fixes

* `calculate_first_minimum_density()` ignored its `adjust` argument (always used 2), and once one
  column had no clear minimum, every later column was skipped as well.

## New features

* `run_noisyr_app()`: interactive Shiny app for the count matrix approach.
* testthat tests that check the optimised code against the original pairwise implementation.
