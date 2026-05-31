# Compute composite assignment score

Combine Mendelian conflict rate and anchor kinship score into a single
composite score for family assignment.

## Usage

``` r
composite_score(conflict_mat, kinship_mat, alpha = 1, beta = 1)
```

## Arguments

- conflict_mat:

  Numeric matrix of Mendelian conflict rates (from
  [`mendelian_conflict()`](https://luansheng.github.io/AAPA/reference/mendelian_conflict.md)).

- kinship_mat:

  Numeric matrix of anchor kinship scores (from
  [`anchor_kinship()`](https://luansheng.github.io/AAPA/reference/anchor_kinship.md)).

- alpha:

  Weight for conflict penalty (default: 1.0).

- beta:

  Weight for kinship reward (default: 1.0).

## Value

A numeric matrix (n_test x n_families) of composite scores. Higher
scores indicate better match.

## Details

\$\$S\_{i,f} = -\alpha \cdot C\_{i,f} + \beta \cdot K\_{i,f}\$\$

## See also

Other scoring:
[`anchor_kinship()`](https://luansheng.github.io/AAPA/reference/anchor_kinship.md),
[`mendelian_conflict()`](https://luansheng.github.io/AAPA/reference/mendelian_conflict.md)
