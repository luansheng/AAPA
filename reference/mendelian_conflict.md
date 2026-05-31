# Compute Mendelian conflict rate

For each test individual and each candidate family, calculate the
proportion of loci where the individual's genotype is incompatible with
the expected offspring genotypes given the parental cross.

## Usage

``` r
mendelian_conflict(genotype, parents, test_ids = NULL)
```

## Arguments

- genotype:

  Numeric matrix (individuals x SNPs), dosage 0/1/2/NA.

- parents:

  An `aapa_parents` object (from
  [`read_parents()`](https://luansheng.github.io/AAPA/reference/read_parents.md)
  or a named list of families with `sire_geno` and `dam_geno`).

- test_ids:

  Character vector of individual IDs to test. If NULL, all non-parent
  individuals in genotype are used.

## Value

A numeric matrix of dimension (n_test x n_families) with Mendelian
conflict rates in `[0, 1]`.

## See also

Other scoring:
[`anchor_kinship()`](https://luansheng.github.io/AAPA/reference/anchor_kinship.md),
[`composite_score()`](https://luansheng.github.io/AAPA/reference/composite_score.md)

## Examples

``` r
sim <- simulate_aapa_data(n_families = 3, n_snps = 100)
#> ℹ Simulating data: 3 families, 100 SNPs
#> ✔ Simulated 47 individuals x 100 markers
tmp <- tempfile(fileext = ".csv")
write.csv(sim$parents, tmp, row.names = FALSE)
parents <- read_parents(tmp, genotype_matrix = sim$genotype)
#> ℹ Reading parents file: /tmp/Rtmp5VLxv7/file196d79549161.csv
#> ✔ Read 3 candidate families
cm <- mendelian_conflict(sim$genotype, parents)
unlink(tmp)
```
