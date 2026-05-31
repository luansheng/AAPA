# Align anchor genotype matrix to genotype markers

Reorder the anchor genotype submatrix to the supplied marker order after
checking that all required anchor IDs and marker IDs are present.

## Usage

``` r
.align_anchors_to_markers(anchors, marker_ids)
```

## Arguments

- anchors:

  An `aapa_anchors` object.

- marker_ids:

  Character vector of marker IDs defining the target order.

## Value

An anchors object aligned to `marker_ids`.
