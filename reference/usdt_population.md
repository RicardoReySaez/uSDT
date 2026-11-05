# Population values for planning a uSDT study

Describes the population a planned study samples from: how awareness and
influence sensitivity (\\d'\\) vary and covary across subjects, and how
the response criteria vary. The arguments take the names of
[`usdt_priors()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_priors.md),
with numbers instead of priors, and the population can be written in
either parameterization of
[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md).

## Usage

``` r
usdt_population(
  dprime = NULL,
  sd_dprime = NULL,
  cor_dprime = NULL,
  criterion = 0,
  sd_criterion = 0.3,
  cor_criterion = 0,
  sd_signal = 1,
  intercept = NULL,
  from = NULL
)

# S3 method for class 'usdt_population'
print(x, ...)

# S3 method for class 'usdt_population'
update(object, ...)
```

## Arguments

- dprime:

  Mean \\d'\\ of each task. Under the latent regression parameterization
  (when `intercept` is given), the mean of the awareness task only.

- sd_dprime:

  Between-subject standard deviation of \\d'\\ in each task.

- cor_dprime:

  Correlation between awareness and influence \\d'\\ across subjects.

- criterion:

  Mean criterion \\c\\ of each task, on the classical scale where
  positive values are conservative.

- sd_criterion:

  Between-subject standard deviation of the criterion in each task. Zero
  gives every subject the same criterion.

- cor_criterion:

  Correlation between the two criteria across subjects.

- sd_signal:

  Standard deviation of the signal distribution of each task, relative
  to the noise distribution. The default, 1, is the equal-variance model
  that
  [`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md)
  fits.

- intercept:

  Intercept of the latent regression of influence on awareness \\d'\\.
  Giving it selects the latent regression parameterization, in which it
  replaces the influence mean.

- from:

  A frequentist
  [`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md)
  fit, such as a pilot study, whose estimates give every value not set
  in the call.

- x:

  A `usdt_population` object.

- ...:

  For [`update()`](https://rdrr.io/r/stats/update.html), the arguments
  of `usdt_population()` to change. A task argument may name one task
  only, which keeps the value of the other.

- object:

  A `usdt_population` object.

## Value

An object of class `usdt_population`: a data frame with one row per
scenario and the columns `scenario`, `mu_A`, `mu_I`, `sigma_A`,
`sigma_I`, `rho`, the tested quantities `diff`, `intercept` and `slope`,
and `criterion_A`, `criterion_I`, `sd_criterion_A`, `sd_criterion_I`,
`cor_criterion`, `sd_signal_A` and `sd_signal_I`. Its `latent` attribute
is `"correlation"` or `"regression"`.

## Values for each task

`dprime`, `sd_dprime`, `criterion`, `sd_criterion` and `sd_signal` take
one value for both tasks or one for each, as
`c(awareness = , influence = )` or `list(awareness = , influence = )`.
The correlations and the intercept take a single value.

## The two parameterizations

The correlation parameterization gives the bivariate normal distribution
of \\d'\\: its means \\\mu_A\\ and \\\mu_I\\, standard deviations
\\\sigma_A\\ and \\\sigma_I\\ and correlation \\\rho\\. The latent
regression parameterization gives the intercept \\\beta_0\\ instead of
\\\mu_I\\, which then follows as \\\mu_I = \beta_0 + \rho \\ (\sigma_I /
\sigma_A) \\ \mu_A\\, as in `hsdt(latent = "regression")`. Both describe
the same distribution, and the printed population reports the tested
quantities of either: the difference \\\mu_I - \mu_A\\ (H1), the
correlation (H2) and the intercept and slope of the latent regression
(H3).

The parameterization decides what stays fixed when a value changes, in a
scenario or in [`update()`](https://rdrr.io/r/stats/update.html). Under
the correlation parameterization a new correlation keeps both means and
moves the intercept (unless \\\mu_A = 0\\). Under the latent regression
it keeps the intercept and moves the influence mean, and with it the
difference.

## Criteria

The criteria matter only through the precision with which each subject's
trials measure \\d'\\. A task split at the median has its criterion
fixed at zero by the split, whatever its value here. The criterion
correlation applies when both criteria are estimated.

## Scenarios

Several values in an argument give one scenario for each, and several
arguments with several values give one scenario for each combination. A
value shared by both tasks keeps them equal in every scenario:
`sd_dprime = c(0.2, 0.3)` gives two scenarios, while
`sd_dprime = list(awareness = c(0.2, 0.3), influence = c(0.2, 0.3))`
gives four.

## Values from a fitted model

With `from`, the estimates of a frequentist
[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md) fit
give the means, standard deviations and correlations, and the arguments
set in the call replace them, so
`usdt_population(from = m, cor_dprime = 0.3)` keeps every estimate of
`m` except the correlation. Setting `intercept` switches to the latent
regression parameterization and drops the estimated influence mean. A
criterion fixed at zero by the median split is zero with no variation.

## See also

[`usdt_priors()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_priors.md),
[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md)

## Examples

``` r
# The correlation parameterization: means, SDs and correlation of d'
pop <- usdt_population(dprime    = c(awareness = 0.15, influence = 0.30),
                       sd_dprime = c(awareness = 0.35, influence = 0.30),
                       cor_dprime = 0.3)
pop
#> ── Population ─────────────────────────────────────────────────── correlation ── 
#> 
#>   Parameter     Awareness  Influence
#>   mean d'           0.150      0.300
#>   sd(d')            0.350      0.300
#>   cor(d')           0.300           
#>   mean c            0.000      0.000
#>   sd(c)             0.300      0.300
#> 
#>        Tested quantity                 Value
#>   H1   Δd' (influence - awareness)     0.150
#>   H2   correlation                     0.300
#>   H3   intercept                       0.261
#>        slope                           0.257
#> 
#>   A task split at the median has its criterion fixed at 0.

# The latent regression parameterization: the intercept replaces the
# influence mean
usdt_population(dprime = 0.15, intercept = 0.25,
                sd_dprime = c(awareness = 0.35, influence = 0.30),
                cor_dprime = 0.3)
#> ── Population ───────────────────────────────────────────── latent regression ── 
#> 
#>   Parameter     Awareness  Influence
#>   mean d'           0.150     0.289*
#>   sd(d')            0.350      0.300
#>   cor(d')           0.300           
#>   intercept         0.250           
#>   mean c            0.000      0.000
#>   sd(c)             0.300      0.300
#> 
#>        Tested quantity                 Value
#>   H1   Δd' (influence - awareness)     0.139
#>   H2   correlation                     0.300
#>   H3   intercept                       0.250
#>        slope                           0.257
#> 
#>   * The influence mean follows from the intercept: mu_I = b0 + rho (sd_I / sd_A) mu_A.
#>   A task split at the median has its criterion fixed at 0.

# Three scenarios for the correlation
update(pop, cor_dprime = c(0.1, 0.3, 0.5))
#> ── Population ───────────────────────────────────── correlation · 3 scenarios ── 
#> 
#>   Parameter     Awareness  Influence
#>   mean d'           0.150      0.300
#>   sd(d')            0.350      0.300
#>   cor(d')          varies           
#>   mean c            0.000      0.000
#>   sd(c)             0.300      0.300
#> 
#>   Scenario    cor(d')       Δd' intercept     slope
#>   1             0.100     0.150     0.287     0.086
#>   2             0.300     0.150     0.261     0.257
#>   3             0.500     0.150     0.236     0.429
#> 
#>   A task split at the median has its criterion fixed at 0.

# The estimates of a pilot study, with a weaker correlation
d <- usdt_data_tasks(
  awareness = vadillo_awareness,
  influence = vadillo_cuing,
  subject_col      = "subj",
  condition_col    = "condition",
  condition_levels = c(signal = "old", noise = "new"),
  response_col     = list(awareness = "judged.old", influence = "rt"),
  response_levels  = list(awareness = c(signal = 1, noise = 0),
                          influence = c(signal = "faster",
                                        noise = "slower")),
  dichotomize      = list(awareness = FALSE, influence = TRUE)
)
usdt_population(from = hsdt(d), cor_dprime = 0.3)
#> ── Population ─────────────────────────────────────────────────── correlation ── 
#> 
#>   Parameter     Awareness  Influence
#>   mean d'           0.235      0.128
#>   sd(d')            0.122      0.088
#>   cor(d')           0.300           
#>   mean c           -0.036      0.000
#>   sd(c)             0.247      0.000
#> 
#>        Tested quantity                 Value
#>   H1   Δd' (influence - awareness)    -0.107
#>   H2   correlation                     0.300
#>   H3   intercept                       0.077
#>        slope                           0.217
#> 
#>   A task split at the median has its criterion fixed at 0.
```
