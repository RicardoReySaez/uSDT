# Dichotomize response times into binary choices

Splits response times (or other continuous measures) at each subject's
overall median, following the preprocessing approach of Meyen et al.
(2022). The median is calculated across all trials for each participant
without distinguishing between stimulus conditions or other covariates.
This produces a binary outcome that allows response times to be mapped
onto a Signal Detection Theory sensitivity metric (\\d'\\).

## Usage

``` r
meyen_split(x, by, signal = c("faster", "slower"), ties = c("noise", "random"))
```

## Arguments

- x:

  Numeric vector of continuous values, typically response times.

- by:

  Vector identifying the subject for each observation in `x`. Medians
  are computed independently for each participant.

- signal:

  Character string specifying which side of the median will be treated
  as the "signal" response under an SDT framework. Use `"faster"` when
  the target condition speeds up responses (e.g., facilitatory priming,
  spatial cueing) or `"slower"` when it slows responses down (e.g.,
  interference, Stroop-like effects).

- ties:

  How to handle trials that match the subject's median exactly.
  `"noise"` assigns them to the noise category (0). `"random"` breaks
  ties at random, keeping cell proportions as balanced as possible.

## Value

An integer vector of `0`s (noise response) and `1`s (signal response)
matching the length of `x`. Missing values (`NA`) are preserved.

## Details

With an odd number of trials (\\n\\), a dataset cannot be split into two
equal halves because the median falls exactly on an observed trial.
Setting `ties = "noise"` assigns this middle trial to noise, producing a
signal proportion of \\(n - 1) / (2\cdot n)\\ and slightly shifting the
response criterion. In practice, this difference (\\1 / (2\cdot n)\\) is
negligible, but setting `ties = "random"` resolves ties
probabilistically to avoid any systematic directional bias.

## References

Meyen, S., Zerweck, I. A., Amado, C., von Luxburg, U., & Franz, V. H.
(2022). Advancing research on unconscious priming: When can scientists
claim an indirect task advantage? *Journal of Experimental Psychology:
General*, 151(1), 65–81. \<10.1037/xge0001065\>

## Examples

``` r
rt   <- c(320, 410, 295, 500, 380, 450)
subj <- rep(c("s1", "s2"), each = 3)
meyen_split(rt, by = subj)
#> [1] 0 0 1 0 1 0
```
