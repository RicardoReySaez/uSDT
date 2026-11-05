# Priors for the Bayesian hierarchical SDT model

Specifies the prior distributions of the Bayesian model. Each prior is
written as a string in Stan notation, `"family(arguments)"`, with the
families and parameterisations of Stan, e.g. `"normal(0, 1)"`,
`"student_t(4, 0, 0.5)"` or `"gamma(2, 4)"`. Each argument takes the
families that cover the values its parameter can take:

- means of \\d'\\ and of the criterion can take any real value, so their
  prior must be defined on the whole real line: `normal`, `std_normal`,
  `student_t`, `cauchy`, `double_exponential`, `logistic`, `gumbel`,
  `skew_normal`, `exp_mod_normal` or `skew_double_exponential`.

- standard deviations take any of those, truncated at zero (with
  location 0, a `student_t` is a half-Student-\\t\\), or a family on the
  positive values: `lognormal`, `exponential`, `gamma`, `inv_gamma`,
  `weibull`, `frechet`, `rayleigh`, `chi_square`, `inv_chi_square`,
  `scaled_inv_chi_square`, `loglogistic`, or `pareto_type_2` with
  `mu <= 0`.

- correlations take `"lkj(eta)"`, the LKJ prior on a correlation matrix
  and the only correlation prior brms accepts. On a 2 x 2 matrix the
  correlation follows a beta(eta, eta) distribution stretched from \\(0,
  1)\\ to \\(-1, 1)\\: `lkj(1)` is uniform, and larger values of eta
  concentrate it around zero.

## Usage

``` r
usdt_priors(
  dprime = "normal(0, 1)",
  sd_dprime = "student_t(4, 0, 0.5)",
  cor_dprime = "lkj(1)",
  criterion = "normal(0, 0.5)",
  sd_criterion = "student_t(4, 0, 0.5)",
  cor_criterion = "lkj(1)",
  sd_signal = "lognormal(0, 0.5)",
  intercept = "normal(0, 1)"
)

# S3 method for class 'usdt_priors'
print(x, ...)
```

## Arguments

- dprime:

  Prior of the mean \\d'\\ of each task.

- sd_dprime:

  Prior of the between-subject standard deviation of \\d'\\ in each
  task.

- cor_dprime:

  Prior of the correlation between awareness and influence \\d'\\ across
  subjects.

- criterion:

  Prior of the mean criterion \\c\\ of each task, on the classical scale
  where positive values are conservative.

- sd_criterion:

  Prior of the between-subject standard deviation of the criterion in
  each task.

- cor_criterion:

  Prior of the correlation between the two criteria. It is used only
  when both criteria are estimated.

- sd_signal:

  Prior of the standard deviation of the signal distribution, relative
  to the noise distribution. It is used only under unequal variances.

- intercept:

  Prior of the intercept of the latent regression of influence on
  awareness \\d'\\. It is used only by the latent regression
  parameterization, where it replaces the prior of the influence mean
  \\d'\\.

- x:

  A `usdt_priors` object.

- ...:

  Ignored.

## Value

An object of class `usdt_priors`: a data frame with one row per prior,
giving the parameter, the task, the family, its arguments (a list
column) and the prior as a string.

## Details

`dprime`, `sd_dprime`, `criterion`, `sd_criterion` and `sd_signal` take
one prior for both tasks or one per task as
`list(awareness = ..., influence = ...)`. The correlations take a single
prior.

The defaults are centred at zero and weakly informative on the probit
scale: `normal(0, 1)` for the mean \\d'\\, `normal(0, 0.5)` for the mean
criterion, `student_t(4, 0, 0.5)` for every between-subject standard
deviation and a uniform `lkj(1)` for every correlation. With them, the
hit and false-alarm rates of a new subject spread over \\(0, 1)\\
without piling up at either end. See
[`plot.usdt_priors()`](https://ricardoreysaez.github.io/uSDT/reference/plot.usdt_priors.md).
With the same prior for both means, the prior of the difference between
them (H1) is symmetric around zero, and every `lkj` prior is centred at
zero (H2).

The model describes the two sensitivities by their means, standard
deviations and correlation, so the slope and intercept of the latent
regression (H3) have no priors of their own: their priors follow from
those of the means, standard deviations and correlation. With the
defaults, the intercept prior is symmetric around zero but has heavy
tails. It does not change when both standard deviation priors are
rescaled by the same factor, because the slope depends on their ratio.

The latent regression parameterization moves the prior of the influence
mean to the intercept \\\beta_0\\ of the latent regression, so that both
the intercept (H3) and the correlation (H2) have the priors you give
them. It keeps the priors of the awareness mean, of both standard
deviations and of the correlation, and the influence mean follows from
them as \\\mu_I = \beta_0 + \rho \\ (\sigma_I / \sigma_A) \\ \mu_A\\.
`dprime` then applies to the awareness task only. The intercept takes
`normal(0, 1)` by default. Inspect its implications with
[`plot.usdt_priors()`](https://ricardoreysaez.github.io/uSDT/reference/plot.usdt_priors.md).
Two things follow:

- the difference between the means (H1) has an induced prior,
  \\\beta_0 + (\rho \sigma_I / \sigma_A - 1) \mu_A\\, with the heavy
  tails the intercept has under the correlation parameterization.

- the intercept and the correlation are independent a priori but not a
  posteriori. The data fix the influence mean, the sum of \\\beta_0\\
  and \\\rho (\sigma_I / \sigma_A) \mu_A\\, so a narrow intercept prior
  also narrows the posterior of \\\rho\\, mostly when \\\sigma_A\\ is
  small. The Bayes factor of a correlation hypothesis then depends on
  the intercept prior too.

## Examples

``` r
# Default priors
usdt_priors()
#> ── Priors ────────────────────────────────────────────────────────────────────── 
#> 
#>   Parameter    Task       Prior
#>   mean d'      awareness  normal(0, 1)
#>   mean d'      influence  normal(0, 1)
#>   sd(d')       awareness  student_t(4, 0, 0.5)
#>   sd(d')       influence  student_t(4, 0, 0.5)
#>   cor(d')      both       lkj(1)
#>   mean c       awareness  normal(0, 0.5)
#>   mean c       influence  normal(0, 0.5)
#>   sd(c)        awareness  student_t(4, 0, 0.5)
#>   sd(c)        influence  student_t(4, 0, 0.5)
#>   cor(c)       both       lkj(1)
#>   sd(signal)   awareness  lognormal(0, 0.5)
#>   sd(signal)   influence  lognormal(0, 0.5)
#>   intercept    both       normal(0, 1)
#> 
#>   Standard deviation priors are truncated at zero. cor(c) applies when
#>   both criteria are estimated, and sd(signal) under unequal variances.
#>   Under the latent regression, intercept replaces the influence mean d'.

# A narrower prior for the mean influence sensitivity
usdt_priors(dprime = list(awareness = "normal(0, 1)",
                          influence = "normal(0, 0.5)"))
#> ── Priors ────────────────────────────────────────────────────────────────────── 
#> 
#>   Parameter    Task       Prior
#>   mean d'      awareness  normal(0, 1)
#>   mean d'      influence  normal(0, 0.5)
#>   sd(d')       awareness  student_t(4, 0, 0.5)
#>   sd(d')       influence  student_t(4, 0, 0.5)
#>   cor(d')      both       lkj(1)
#>   mean c       awareness  normal(0, 0.5)
#>   mean c       influence  normal(0, 0.5)
#>   sd(c)        awareness  student_t(4, 0, 0.5)
#>   sd(c)        influence  student_t(4, 0, 0.5)
#>   cor(c)       both       lkj(1)
#>   sd(signal)   awareness  lognormal(0, 0.5)
#>   sd(signal)   influence  lognormal(0, 0.5)
#>   intercept    both       normal(0, 1)
#> 
#>   Standard deviation priors are truncated at zero. cor(c) applies when
#>   both criteria are estimated, and sd(signal) under unequal variances.
#>   Under the latent regression, intercept replaces the influence mean d'.

# A correlation prior that concentrates near zero
usdt_priors(cor_dprime = "lkj(2)")
#> ── Priors ────────────────────────────────────────────────────────────────────── 
#> 
#>   Parameter    Task       Prior
#>   mean d'      awareness  normal(0, 1)
#>   mean d'      influence  normal(0, 1)
#>   sd(d')       awareness  student_t(4, 0, 0.5)
#>   sd(d')       influence  student_t(4, 0, 0.5)
#>   cor(d')      both       lkj(2)
#>   mean c       awareness  normal(0, 0.5)
#>   mean c       influence  normal(0, 0.5)
#>   sd(c)        awareness  student_t(4, 0, 0.5)
#>   sd(c)        influence  student_t(4, 0, 0.5)
#>   cor(c)       both       lkj(1)
#>   sd(signal)   awareness  lognormal(0, 0.5)
#>   sd(signal)   influence  lognormal(0, 0.5)
#>   intercept    both       normal(0, 1)
#> 
#>   Standard deviation priors are truncated at zero. cor(c) applies when
#>   both criteria are estimated, and sd(signal) under unequal variances.
#>   Under the latent regression, intercept replaces the influence mean d'.

# Standard deviations of d' expected around 0.3
usdt_priors(sd_dprime = "student_t(4, 0.3, 0.2)")
#> ── Priors ────────────────────────────────────────────────────────────────────── 
#> 
#>   Parameter    Task       Prior
#>   mean d'      awareness  normal(0, 1)
#>   mean d'      influence  normal(0, 1)
#>   sd(d')       awareness  student_t(4, 0.3, 0.2)
#>   sd(d')       influence  student_t(4, 0.3, 0.2)
#>   cor(d')      both       lkj(1)
#>   mean c       awareness  normal(0, 0.5)
#>   mean c       influence  normal(0, 0.5)
#>   sd(c)        awareness  student_t(4, 0, 0.5)
#>   sd(c)        influence  student_t(4, 0, 0.5)
#>   cor(c)       both       lkj(1)
#>   sd(signal)   awareness  lognormal(0, 0.5)
#>   sd(signal)   influence  lognormal(0, 0.5)
#>   intercept    both       normal(0, 1)
#> 
#>   Standard deviation priors are truncated at zero. cor(c) applies when
#>   both criteria are estimated, and sd(signal) under unequal variances.
#>   Under the latent regression, intercept replaces the influence mean d'.

# Any family that covers the positive values, e.g. a gamma distribution
usdt_priors(sd_dprime = "gamma(2, 4)")
#> ── Priors ────────────────────────────────────────────────────────────────────── 
#> 
#>   Parameter    Task       Prior
#>   mean d'      awareness  normal(0, 1)
#>   mean d'      influence  normal(0, 1)
#>   sd(d')       awareness  gamma(2, 4)
#>   sd(d')       influence  gamma(2, 4)
#>   cor(d')      both       lkj(1)
#>   mean c       awareness  normal(0, 0.5)
#>   mean c       influence  normal(0, 0.5)
#>   sd(c)        awareness  student_t(4, 0, 0.5)
#>   sd(c)        influence  student_t(4, 0, 0.5)
#>   cor(c)       both       lkj(1)
#>   sd(signal)   awareness  lognormal(0, 0.5)
#>   sd(signal)   influence  lognormal(0, 0.5)
#>   intercept    both       normal(0, 1)
#> 
#>   Standard deviation priors are truncated at zero. cor(c) applies when
#>   both criteria are estimated, and sd(signal) under unequal variances.
#>   Under the latent regression, intercept replaces the influence mean d'.
```
