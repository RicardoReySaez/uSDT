# Plot the prior and posterior of tested quantities

Draws, for each hypothesis of a
[`usdt_bf()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_bf.md)
result, the posterior density (solid) and the prior density (dashed) of
the tested quantity, which names the x axis of its panel: \\\Delta d'\\,
\\\rho\\, or \\\beta_0\\ and \\\beta_1\\ for the intercept and slope of
the latent regression. The panel title says whether the hypothesis is
the null (H0), as a point is, or the alternative (H1), as an inequality
or a region is, and the panel reports the Bayes factor of that
hypothesis against the other: \\\mathrm{BF}\_{01}\\ for a point,
\\\mathrm{BF}\_{10}\\ otherwise. For a point hypothesis, a dotted line
marks the tested value and two points mark the heights whose ratio is
the Savage-Dickey Bayes factor. For a directional or interval
hypothesis, the shaded area is the posterior probability of H1. Each
panel has a density scale of its own, so the panels share one y axis,
labelled Density, drawn on the first panel of each row and without
values.

## Usage

``` r
# S3 method for class 'usdt_bf'
plot(x, caption = TRUE, ...)
```

## Arguments

- x:

  A `usdt_bf` object from
  [`usdt_bf()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_bf.md).

- caption:

  Logical. Show the caption that explains the plot under it (default
  `TRUE`).

- ...:

  Ignored.

## Value

A `ggplot` object. Its underlying data frame is stored in `$data`.

## Details

[`usdt_bf()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_bf.md)
draws this figure by default. Call
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) to draw it
again, or to draw some of its rows, e.g. `plot(b[3, ])`.

## See also

[`usdt_bf()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_bf.md)
