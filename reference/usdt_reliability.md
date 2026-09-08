# Reliability of the direct and indirect measures

Estimates how much of the spread in `d'` comes from real differences
between subjects, and how much comes from the noise of a limited number
of trials. A value close to one means that the task separates subjects
well. A value close to zero means that most of the observed spread is
measurement error.

## Usage

``` r
usdt_reliability(object)

# S3 method for class 'usdt_reliability'
summary(object, ...)

# S3 method for class 'usdt_reliability'
print(x, ...)
```

## Arguments

- object:

  An `hsdt` object from
  [`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md),
  which may also carry the results of
  [`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_boot.md).
  The summary method takes the `usdt_reliability` object that this
  function returns.

- ...:

  Ignored.

- x:

  A `usdt_reliability` object.

## Value

An object of class `usdt_reliability`. Its `tasks` data frame gives one
reliability per task, with the two variances it comes from. Its
`subjects` data frame gives the `d'`, the measurement variance and the
reliability of every subject. Bootstrap intervals appear in both when
they are available.

## Details

For subject `j` in task `t` the reported value is
`tau2_t / (tau2_t + v_tj)`. The first term, `tau2_t`, is the variance of
the sensitivity between subjects, which the model estimates. The second
term, `v_tj`, is the variance of the sensitivity that the trials of that
subject can support on their own.

`v_tj` depends on the number of trials and on the position of the
subject on the response curve. A cell contributes
`n * dnorm(eta)^2 / (p * (1 - p))`, so two subjects with the same number
of trials can differ in precision. The calculation also discounts the
information spent on estimating the criterion. When the model has fixed
a criterion to zero there is nothing to discount, and the same formula
applies.

This measure looks at each task alone and never uses the other task, so
it describes the information the data of one subject actually carry. The
estimates drawn by
[`plot.hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/plot.hsdt.md)
are more precise than this, because the model there does borrow
information across subjects and tasks.

The value reported for a whole task replaces `v_tj` by its average. For
a subject drawn at random, the variance of a single measurement is
`tau2_t + E(v_tj)`, so the ratio says which share of that spread is
real. Averaging the reliabilities of the individual subjects would
answer a different question. The average over subjects treats them all
as equally important, which suits a sample that represents the
population of interest.

`v_tj` is the classical sampling variance of `d'` of Gourevitch and
Galanter (1967), which
[`sdt_moments()`](https://ricardoreysaez.github.io/uSDT/reference/sdt_moments.md)
reports as `var_gg`. It also equals the standard error that a probit
regression would give for the sensitivity of that subject alone. The
three are the same formula, evaluated at different points. `var_gg` uses
the observed rates of the subject, while `v_tj` uses the rates that the
hierarchical model predicts.

When `object` carries a usable bootstrap, the function recomputes the
reliability in every retained replicate, each one with its own variance
and its own subject estimates. Percentile intervals remain available at
a boundary. Normal and basic intervals work on the logit scale, so they
are missing when a reliability reaches zero or one.

## References

Gourevitch, V., & Galanter, E. (1967). A significance test for one
parameter isosensitivity functions. *Psychometrika*.

## See also

[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md),
[`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_boot.md),
[`sdt_moments()`](https://ricardoreysaez.github.io/uSDT/reference/sdt_moments.md)

## Examples

``` r
# \donttest{
set.seed(1)
df <- usdt_simulate(n_subj = 30, n_trials = 80)
data <- usdt_data_long(
  df,
  task_col = "task",
  task_levels = c(direct = "D", indirect = "I"),
  subject_col = "subj",
  condition_col = "cond",
  condition_levels = c(signal = 1, noise = 0),
  response_col = "response",
  response_levels = c(signal = 1, noise = 0)
)
usdt_reliability(hsdt(data))
#> Warning: H2 and H3 are not estimable with reliable standard errors.
#>   At the bound: subj.d_I.
#>   Reason: the sensitivity covariance is on the boundary
#>   Use usdt_boot() for a parametric bootstrap.
#> ── Reliability summary ───────────────────────────────────────────────────────── 
#> 
#>   Task         Subjects Group-level estimate   By-subject estimate median [min, max]
#>   Direct             30                 .765   .767 [.727, .783]
#>   Indirect           30                 .017   .017 [.016, .018]
# }
```
