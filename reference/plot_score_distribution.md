# Plot score distribution

Visualize the distribution of composite scores across all test
individuals and candidate families, highlighting assigned vs. rejected
individuals.

## Usage

``` r
plot_score_distribution(result, type = "histogram")
```

## Arguments

- result:

  An `aapa_result` object.

- type:

  One of `"histogram"` or `"density"`. Default: `"histogram"`.

## Value

A ggplot2 object (if ggplot2 is available), otherwise a base R plot is
produced and NULL is returned invisibly.

## See also

Other visualization:
[`plot_rejection_diagnostics()`](https://luansheng.github.io/AAPA/reference/plot_rejection_diagnostics.md),
[`plot_topk()`](https://luansheng.github.io/AAPA/reference/plot_topk.md)
