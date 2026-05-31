# Simulate AAPA example data

Generate a small simulated dataset for testing and demonstration.

## Usage

``` r
simulate_aapa_data(
  n_families = 10,
  n_snps = 500,
  n_offspring_per_family = 10,
  n_anchors_per_family = 2,
  n_unknown = 5,
  missing_rate = 0.01,
  error_rate = 0.001,
  seed = 42
)
```

## Arguments

- n_families:

  Number of candidate families (default: 10).

- n_snps:

  Number of SNP markers (default: 500).

- n_offspring_per_family:

  Number of offspring per family (default: 10).

- n_anchors_per_family:

  Number of anchors per family (default: 2).

- n_unknown:

  Number of unknown-family test individuals (default: 5).

- missing_rate:

  Proportion of missing genotypes (default: 0.01).

- error_rate:

  Genotyping error rate (default: 0.001).

- seed:

  Random seed for reproducibility.

## Value

A list with components: `genotype` (matrix), `parents` (data.frame),
`anchors` (data.frame), `true_labels` (named character vector).

## See also

Other data-io:
[`read_anchors()`](https://luansheng.github.io/AAPA/reference/read_anchors.md),
[`read_genotype()`](https://luansheng.github.io/AAPA/reference/read_genotype.md),
[`read_parents()`](https://luansheng.github.io/AAPA/reference/read_parents.md)

## Examples

``` r
sim <- simulate_aapa_data(n_families = 3, n_snps = 100)
#> ℹ Simulating data: 3 families, 100 SNPs
#> ✔ Simulated 47 individuals x 100 markers
str(sim, max.level = 1)
#> List of 4
#>  $ genotype   : int [1:47, 1:100] 1 2 1 1 2 2 1 2 2 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ parents    :'data.frame': 3 obs. of  3 variables:
#>  $ anchors    :'data.frame': 6 obs. of  3 variables:
#>  $ true_labels: Named chr [1:41] "FAM001" "FAM001" "FAM001" "FAM001" ...
#>   ..- attr(*, "names")= chr [1:41] "FAM001_ANC1" "FAM001_ANC2" "FAM001_OFF1" "FAM001_OFF2" ...
```
