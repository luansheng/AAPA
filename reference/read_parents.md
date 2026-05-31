# Read parents table

Read a table of candidate families with parent genotypes.

## Usage

``` r
read_parents(file, genotype_matrix, sep = ",")
```

## Arguments

- file:

  Path to a CSV/TSV file with columns: family_id, sire_id, dam_id.

- genotype_matrix:

  A genotype matrix (from
  [`read_genotype()`](https://luansheng.github.io/AAPA/reference/read_genotype.md))
  that includes sire and dam genotypes.

- sep:

  Field separator (default: comma).

## Value

A list of class `aapa_parents`, where each element is a list with
`family_id`, `sire_id`, `dam_id`, `sire_geno`, `dam_geno`.

## See also

Other data-io:
[`read_anchors()`](https://luansheng.github.io/AAPA/reference/read_anchors.md),
[`read_genotype()`](https://luansheng.github.io/AAPA/reference/read_genotype.md),
[`simulate_aapa_data()`](https://luansheng.github.io/AAPA/reference/simulate_aapa_data.md)
