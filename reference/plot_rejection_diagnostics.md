# Plot rejection diagnostics

Visualize the relationship between score, confidence (gap), and conflict
rate, with rejection threshold lines.

## Usage

``` r
plot_rejection_diagnostics(result)
```

## Arguments

- result:

  An `aapa_result` object.

## Value

A ggplot2 object (if available), otherwise base R plot.

## See also

Other visualization:
[`plot_score_distribution()`](https://luansheng.github.io/AAPA/reference/plot_score_distribution.md),
[`plot_topk()`](https://luansheng.github.io/AAPA/reference/plot_topk.md)
