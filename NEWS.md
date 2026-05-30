# aapa 0.1.0

## New Features

* Initial release of the AAPA (Anchor-Assisted Pedigree Assignment) package.
* Core scoring engine with Mendelian conflict rate computation
  (`mendelian_conflict()`).
* Anchor-based IBS kinship scoring (`anchor_kinship()`).
* Composite score combining conflict and kinship (`composite_score()`).
* Full assignment pipeline with top-k pruning and rejection filtering
  (`aapa_assign()`).
* Data I/O helpers: `read_genotype()`, `read_parents()`, `read_anchors()`.
* Simulation function for testing: `simulate_aapa_data()`.
* Quality control filtering: `qc_filter()`.
* Visualization: `plot_score_distribution()`, `plot_topk()`,
  `plot_rejection_diagnostics()`.

## Infrastructure

* Configured `testthat` edition 3 for unit testing.
* Added GitHub Actions CI workflows for R CMD check, test coverage,
  linting, and pkgdown site deployment.
* Added `lintr` configuration for tidyverse-style code checking.
* Added `pkgdown` site configuration.
* C++ placeholders for Phase 2 Rcpp optimization.
