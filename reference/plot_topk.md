# Plot top-k candidates for an individual

Show the top-k candidate families and their scores for a specific test
individual, useful for diagnosing ambiguous or rejected cases.

## Usage

``` r
plot_topk(result, individual_id)
```

## Arguments

- result:

  An `aapa_result` object.

- individual_id:

  Character; the ID of the individual to plot.

## Value

A ggplot2 object (if available), otherwise base R plot.

## See also

Other visualization:
[`plot_rejection_diagnostics()`](https://luansheng.github.io/AAPA/reference/plot_rejection_diagnostics.md),
[`plot_score_distribution()`](https://luansheng.github.io/AAPA/reference/plot_score_distribution.md)
