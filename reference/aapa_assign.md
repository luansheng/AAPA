# AAPA family assignment

Perform anchor-assisted pedigree assignment for test individuals. This
is the main user-facing function that orchestrates the full pipeline:
Mendelian conflict scoring, anchor kinship scoring, composite scoring,
top-k pruning, and rejection filtering.

## Usage

``` r
aapa_assign(
  genotype,
  parents,
  anchors,
  test_ids = NULL,
  alpha = 1,
  beta = 1,
  top_k = 5,
  tau_conf = -0.05,
  tau_rej = 0.1,
  max_conflict = 0.1
)
```

## Arguments

- genotype:

  Numeric matrix (individuals x SNPs) with dosage encoding 0/1/2/NA.

- parents:

  An `aapa_parents` object or a named list of families (see
  [`read_parents()`](https://luansheng.github.io/AAPA/reference/read_parents.md)).

- anchors:

  An `aapa_anchors` object (see
  [`read_anchors()`](https://luansheng.github.io/AAPA/reference/read_anchors.md)).

- test_ids:

  Character vector of individual IDs to assign. If NULL, all individuals
  not in parents or anchors are used.

- alpha:

  Weight for Mendelian conflict penalty (default: 1.0).

- beta:

  Weight for anchor kinship reward (default: 1.0).

- top_k:

  Number of top candidate families to retain (default: 5).

- tau_conf:

  Minimum score threshold for confident assignment (default: -0.05).
  Individuals whose best score falls below this are rejected.

- tau_rej:

  Minimum score gap between top-1 and top-2 families (default: 0.1).
  Individuals with smaller gaps are rejected.

- max_conflict:

  Maximum allowable Mendelian conflict rate (default: 0.1). Individuals
  whose best family exceeds this are rejected.

## Value

An object of class `aapa_result`, a list containing:

- assignment:

  Data frame with columns: individual_id, assigned_family, score,
  confidence, status (ASSIGNED/REJECT).

- topk:

  List of per-individual top-k candidate data frames.

- conflict_matrix:

  Full Mendelian conflict rate matrix.

- kinship_matrix:

  Full anchor kinship score matrix.

- score_matrix:

  Full composite score matrix.

- params:

  List of parameters used.

## See also

Other assignment:
[`print.aapa_result()`](https://luansheng.github.io/AAPA/reference/print.aapa_result.md),
[`summary.aapa_result()`](https://luansheng.github.io/AAPA/reference/summary.aapa_result.md)

## Examples

``` r
# Simulate data and run the full pipeline
sim <- simulate_aapa_data(n_families = 3, n_snps = 100, seed = 1)
#> ℹ Simulating data: 3 families, 100 SNPs
#> ✔ Simulated 47 individuals x 100 markers

# Write temp files for reading functions
tmp_parents <- tempfile(fileext = ".csv")
write.csv(sim$parents, tmp_parents, row.names = FALSE)
parents <- read_parents(tmp_parents, sim$genotype)
#> ℹ Reading parents file: /tmp/Rtmp5VLxv7/file196d2d7dc971.csv
#> ✔ Read 3 candidate families

tmp_anchors <- tempfile(fileext = ".csv")
write.csv(sim$anchors, tmp_anchors, row.names = FALSE)
anchors <- read_anchors(tmp_anchors, sim$genotype)
#> ℹ Reading anchors file: /tmp/Rtmp5VLxv7/file196d70314ae9.csv
#> ✔ Read 6 anchors from 3 families

result <- aapa_assign(sim$genotype, parents, anchors)
#> ℹ Step 1/4: Computing Mendelian conflict rates...
#> ℹ Step 2/4: Computing anchor kinship scores...
#> ℹ Step 3/4: Computing composite scores...
#> ℹ Step 4/4: Top-k pruning and rejection filtering...
#> ✔ Assignment complete: 30 assigned, 5 rejected
print(result)
#> 
#> ── AAPA Assignment Result ──────────────────────────────────────────────────────
#> Total individuals: 35
#> Assigned: 30 (85.7%)
#> Rejected: 5 (14.3%)
#> Families: 3
#> Parameters: alpha=1, beta=1, top_k=5

unlink(c(tmp_parents, tmp_anchors))
```
