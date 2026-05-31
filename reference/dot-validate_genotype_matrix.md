# Validate genotype matrix contract

Check that a genotype matrix satisfies the AAPA object contract: unique
sample IDs, unique marker IDs, and dosage encoding 0/1/2/NA.

## Usage

``` r
.validate_genotype_matrix(genotype, arg = "genotype")
```

## Arguments

- genotype:

  A genotype matrix.

- arg:

  Argument name used in error messages.

## Value

The input matrix, invisibly.
