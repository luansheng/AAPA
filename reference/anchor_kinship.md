# Compute anchor kinship scores

For each test individual and each candidate family, compute the mean IBS
(Identity-By-State) similarity with the family's anchor individuals.

## Usage

``` r
anchor_kinship(genotype, anchors, test_ids, method = "ibs")
```

## Arguments

- genotype:

  Numeric matrix (individuals x SNPs).

- anchors:

  An `aapa_anchors` object (from
  [`read_anchors()`](https://luansheng.github.io/AAPA/reference/read_anchors.md)).

- test_ids:

  Character vector of test individual IDs.

- method:

  Kinship estimation method. Currently only `"ibs"` (proportion of IBS
  matches) is supported.

## Value

A numeric matrix (n_test x n_families) with kinship scores.

## See also

Other scoring:
[`composite_score()`](https://luansheng.github.io/AAPA/reference/composite_score.md),
[`mendelian_conflict()`](https://luansheng.github.io/AAPA/reference/mendelian_conflict.md)
