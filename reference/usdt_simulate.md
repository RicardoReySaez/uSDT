# Simulate data from a hierarchical SDT model

Generates trial-level data for a direct and an indirect task. Every
subject has one sensitivity in each task, and the two sensitivities
correlate across subjects. This is the structure that
[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md)
estimates, so the function is useful to check an analysis before
collecting data, or to study the power of a planned design.

## Usage

``` r
usdt_simulate(
  n_subj = 50,
  n_trials = 100,
  gamma_D = 0.8,
  gamma_I = 0.3,
  sd_D = 0.5,
  sd_I = 0.3,
  rho = 0.5,
  crit_D = 0,
  sd_crit = 0.3,
  crit_I = 0,
  sd_crit_I = 0.3,
  rt = FALSE
)
```

## Arguments

- n_subj:

  Number of subjects.

- n_trials:

  Number of trials per subject and task. Half of them belong to the
  signal condition and half to the noise condition.

- gamma_D, gamma_I:

  Average sensitivity (d') of the direct and the indirect task.

- sd_D, sd_I:

  How much each sensitivity varies between subjects.

- rho:

  Correlation between the two sensitivities across subjects.

- crit_D, sd_crit:

  Model intercept of the direct task and how much it varies between
  subjects. The classical criterion has the opposite sign,
  `c = -crit_D`. See Details.

- crit_I, sd_crit_I:

  The same two values for the indirect task. The function sets both to
  zero when `rt = TRUE`, because a median split leaves no criterion to
  estimate. For a binary indirect task, keep `sd_crit_I` above zero. A
  value of zero gives every subject the same criterion, and the model
  then becomes singular when it tries to estimate that variation.

- rt:

  If `TRUE`, the indirect task returns response times instead of binary
  responses, so that the median split of
  [`meyen_split()`](https://ricardoreysaez.github.io/uSDT/reference/meyen_split.md)
  has something to work on. The times fall as the evidence for signal
  rises, so faster responses correspond to signal. See Details.

## Value

A data frame with one row per trial and the columns `subj`, `task`
(`"D"` or `"I"`), `cond` (`1` for signal, `0` for noise) and `response`.
With `rt = TRUE` it also has a numeric `rt` column, and `response` is
`NA` in the indirect task.

## Details

The linear predictor is `crit_D + d' * S`, so `crit_D` is the model
intercept, which
[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md)
reports as `c_D`. The classical criterion has the opposite sign under
deviation coding, `c = -crit_D`. That is the value
[`sdt_moments()`](https://ricardoreysaez.github.io/uSDT/reference/sdt_moments.md)
returns and the value the printed summary shows. Simulating
`crit_D = 0.5` therefore gives a criterion of `-0.5`.

## Simulated response times

With `rt = TRUE` the indirect task returns `exp(6.2 - 0.25 * e)`, where
`e` is the latent evidence of that trial. The logarithm of the time is
normal, so the time itself follows a lognormal distribution, and it
decreases as the evidence grows. A faster response is therefore the
signal response.

The two constants only set the scale in milliseconds. They place the
median at `exp(6.2)`, around 493 ms, with the middle 95% of times
between 302 and 806 ms. A median split gives the same result under any
transformation that preserves the order of the values, so the recovered
sensitivity does not depend on them.

The simulated times carry no variation beyond the evidence itself. This
makes `gamma_I` exactly the sensitivity of the dichotomized measure.
Real response times also vary for reasons unrelated to the
discrimination, and that extra variation would lower the sensitivity
recovered from the split.

## See also

[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md),
[`sdt_moments()`](https://ricardoreysaez.github.io/uSDT/reference/sdt_moments.md)

## Examples

``` r
set.seed(1)
df <- usdt_simulate(n_subj = 20, n_trials = 60)
head(df)
#>   subj task cond response
#> 1    1    D    0        0
#> 2    1    D    0        1
#> 3    1    D    0        0
#> 4    1    D    0        1
#> 5    1    D    0        1
#> 6    1    D    0        0
table(df$task, df$cond)
#>    
#>       0   1
#>   D 600 600
#>   I 600 600
```
