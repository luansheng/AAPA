# Quality control filter for genotype data

Filter SNP markers and individuals based on missing data rates, minor
allele frequency, and other quality metrics.

## Usage

``` r
qc_filter(
  genotype,
  max_snp_missing = 0.1,
  max_ind_missing = 0.2,
  min_maf = 0.01,
  verbose = TRUE
)
```

## Arguments

- genotype:

  Numeric matrix (individuals x SNPs), dosage 0/1/2/NA.

- max_snp_missing:

  Maximum allowable missing rate per SNP (default: 0.1). SNPs exceeding
  this are removed.

- max_ind_missing:

  Maximum allowable missing rate per individual (default: 0.2).
  Individuals exceeding this are removed.

- min_maf:

  Minimum minor allele frequency (default: 0.01). SNPs below this are
  removed.

- verbose:

  Logical; print QC summary? Default TRUE.

## Value

A list with components:

- genotype:

  Filtered genotype matrix.

- removed_snps:

  Names of removed SNPs.

- removed_inds:

  Names of removed individuals.

- summary:

  QC summary statistics.

## Examples

``` r
sim <- simulate_aapa_data(n_families = 3, n_snps = 100)
#> ℹ Simulating data: 3 families, 100 SNPs
#> ✔ Simulated 47 individuals x 100 markers
qc_result <- qc_filter(sim$genotype, verbose = FALSE)
dim(qc_result$genotype)
#> [1] 47 98
```
