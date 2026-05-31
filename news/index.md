# Changelog

## aapa 0.1.0

### New Features

- Initial release of the AAPA (Anchor-Assisted Pedigree Assignment)
  package.
- Core scoring engine with Mendelian conflict rate computation
  ([`mendelian_conflict()`](https://luansheng.github.io/AAPA/reference/mendelian_conflict.md)).
- Anchor-based IBS kinship scoring
  ([`anchor_kinship()`](https://luansheng.github.io/AAPA/reference/anchor_kinship.md)).
- Composite score combining conflict and kinship
  ([`composite_score()`](https://luansheng.github.io/AAPA/reference/composite_score.md)).
- Full assignment pipeline with top-k pruning and rejection filtering
  ([`aapa_assign()`](https://luansheng.github.io/AAPA/reference/aapa_assign.md)).
- Data I/O helpers:
  [`read_genotype()`](https://luansheng.github.io/AAPA/reference/read_genotype.md),
  [`read_parents()`](https://luansheng.github.io/AAPA/reference/read_parents.md),
  [`read_anchors()`](https://luansheng.github.io/AAPA/reference/read_anchors.md).
- Simulation function for testing:
  [`simulate_aapa_data()`](https://luansheng.github.io/AAPA/reference/simulate_aapa_data.md).
- Quality control filtering:
  [`qc_filter()`](https://luansheng.github.io/AAPA/reference/qc_filter.md).
- Visualization:
  [`plot_score_distribution()`](https://luansheng.github.io/AAPA/reference/plot_score_distribution.md),
  [`plot_topk()`](https://luansheng.github.io/AAPA/reference/plot_topk.md),
  [`plot_rejection_diagnostics()`](https://luansheng.github.io/AAPA/reference/plot_rejection_diagnostics.md).

### Infrastructure

- Configured `testthat` edition 3 for unit testing.
- Added GitHub Actions CI workflows for R CMD check, test coverage,
  linting, and pkgdown site deployment.
- Added `lintr` configuration for tidyverse-style code checking.
- Added `pkgdown` site configuration.
- C++ placeholders for Phase 2 Rcpp optimization.
