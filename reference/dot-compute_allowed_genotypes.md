# Compute allowed offspring genotypes per locus

Compute allowed offspring genotypes per locus

## Usage

``` r
.compute_allowed_genotypes(sire_geno, dam_geno)
```

## Arguments

- sire_geno:

  Integer vector of sire dosages (0/1/2/NA).

- dam_geno:

  Integer vector of dam dosages (0/1/2/NA).

## Value

A logical matrix (n_loci x 3) where columns correspond to dosage 0, 1, 2
indicating whether each is an allowed offspring genotype.
