# Align parent genotype vectors to genotype markers

Reorder each family's `sire_geno` and `dam_geno` vectors to the supplied
marker order after checking that all markers are present.

## Usage

``` r
.align_parents_to_markers(parents, marker_ids)
```

## Arguments

- parents:

  An `aapa_parents` object or compatible named list.

- marker_ids:

  Character vector of marker IDs defining the target order.

## Value

A parents object aligned to `marker_ids`.
