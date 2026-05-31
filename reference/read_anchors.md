# Read anchors table

Read a table of anchor individuals with known family assignments.

## Usage

``` r
read_anchors(file, genotype_matrix, sep = ",")
```

## Arguments

- file:

  Path to a CSV/TSV file with columns: individual_id, family_id, and
  optionally weight.

- genotype_matrix:

  A genotype matrix that includes anchor genotypes.

- sep:

  Field separator (default: comma).

## Value

A data.frame of class `aapa_anchors` with columns: individual_id,
family_id, weight, and a `geno` attribute containing the anchor genotype
sub-matrix.

## See also

Other data-io:
[`read_genotype()`](https://luansheng.github.io/AAPA/reference/read_genotype.md),
[`read_parents()`](https://luansheng.github.io/AAPA/reference/read_parents.md),
[`simulate_aapa_data()`](https://luansheng.github.io/AAPA/reference/simulate_aapa_data.md)
