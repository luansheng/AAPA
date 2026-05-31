# Validate anchors object contract

Check that an anchors object has required columns, stable IDs, and an
aligned genotype matrix in `attr(anchors, "geno")`.

## Usage

``` r
.validate_anchors_object(anchors, family_ids = NULL)
```

## Arguments

- anchors:

  An `aapa_anchors` object.

- family_ids:

  Optional candidate family IDs used to validate anchor family
  membership.

## Value

The input anchors object, invisibly.
