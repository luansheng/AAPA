# Package index

## Data Input/Output

Functions for reading genotype, parent, and anchor data

- [`read_genotype()`](https://luansheng.github.io/AAPA/reference/read_genotype.md)
  : Read genotype dosage matrix
- [`read_parents()`](https://luansheng.github.io/AAPA/reference/read_parents.md)
  : Read parents table
- [`read_anchors()`](https://luansheng.github.io/AAPA/reference/read_anchors.md)
  : Read anchors table
- [`simulate_aapa_data()`](https://luansheng.github.io/AAPA/reference/simulate_aapa_data.md)
  : Simulate AAPA example data

## Core Scoring

Mendelian conflict, kinship scoring, and composite scoring

- [`mendelian_conflict()`](https://luansheng.github.io/AAPA/reference/mendelian_conflict.md)
  : Compute Mendelian conflict rate
- [`anchor_kinship()`](https://luansheng.github.io/AAPA/reference/anchor_kinship.md)
  : Compute anchor kinship scores
- [`composite_score()`](https://luansheng.github.io/AAPA/reference/composite_score.md)
  : Compute composite assignment score

## Assignment Pipeline

Main assignment, rejection filtering, and result methods

- [`aapa_assign()`](https://luansheng.github.io/AAPA/reference/aapa_assign.md)
  : AAPA family assignment
- [`print(`*`<aapa_result>`*`)`](https://luansheng.github.io/AAPA/reference/print.aapa_result.md)
  : Print method for aapa_result
- [`summary(`*`<aapa_result>`*`)`](https://luansheng.github.io/AAPA/reference/summary.aapa_result.md)
  : Summary method for aapa_result

## Quality Control

Genotype data quality filtering

- [`qc_filter()`](https://luansheng.github.io/AAPA/reference/qc_filter.md)
  : Quality control filter for genotype data

## Visualization

Plotting functions for results exploration

- [`plot_score_distribution()`](https://luansheng.github.io/AAPA/reference/plot_score_distribution.md)
  : Plot score distribution
- [`plot_topk()`](https://luansheng.github.io/AAPA/reference/plot_topk.md)
  : Plot top-k candidates for an individual
- [`plot_rejection_diagnostics()`](https://luansheng.github.io/AAPA/reference/plot_rejection_diagnostics.md)
  : Plot rejection diagnostics
