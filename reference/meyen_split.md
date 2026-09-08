# Turn a continuous measure into a binary response

Splits a continuous measure, usually response times, at each subject's
own median. The median uses all the trials of that subject, without
separating the conditions. The binary result can then be analysed on the
same sensitivity scale as a direct task that already gives binary
responses. The procedure follows Meyen et al. (2022).

## Usage

``` r
meyen_split(x, by, signal = c("faster", "slower"), ties = c("noise", "random"))
```

## Arguments

- x:

  Numeric vector with the continuous measure, usually response times.

- by:

  Vector identifying the subject of each value in `x`. Every subject
  receives their own median.

- signal:

  Which side of the median counts as a signal response. Use `"faster"`
  when the signal condition speeds responses up, as in priming and
  cueing tasks. Use `"slower"` when the signal condition slows responses
  down, as in interference tasks.

- ties:

  What to do with trials that fall exactly on the median. `"noise"`
  gives them the noise response. `"random"` assigns them at random,
  which keeps the split as close to even as the data allow.

## Value

An integer vector as long as `x`, with `1` for signal responses and `0`
for noise responses. Missing values in `x` stay missing.

## Details

A subject with an odd number of trials cannot be split into two equal
halves, because the median is one of the observed values. With
`ties = "noise"` that subject gets a proportion of `(n - 1) / (2n)`
signal responses instead of 0.5, so the criterion moves slightly away
from zero. The difference is `1 / (2n)` and rarely matters. Setting
`ties = "random"` removes it, because rounding up or down at random is
unbiased across subjects.

## References

Meyen, S., Zerweck, I. A., Amado, C., von Luxburg, U., & Franz, V. H.
(2022). Advancing research on unconscious priming: When can scientists
claim an indirect task advantage? *Journal of Experimental Psychology:
General*, 151(1), 65-81.
[doi:10.1037/xge0001065](https://doi.org/10.1037/xge0001065)

## Examples

``` r
rt   <- c(320, 410, 295, 500, 380, 450)
subj <- rep(c("s1", "s2"), each = 3)
meyen_split(rt, by = subj)
#> [1] 0 0 1 0 1 0
```
