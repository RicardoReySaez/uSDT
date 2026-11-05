# Reliability of awareness and influence task measures

Estimates the reliability of sensitivity (\\d'\\) by separating true
variance across participants from sampling noise caused by finite trial
counts. Values close to 1 indicate that the measure reliably separates
participants, whereas values near 0 indicate that observed differences
are mostly measurement error.

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
  frequentist or Bayesian. A frequentist fit may also contain bootstrap
  results from
  [`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_boot.md).

- ...:

  Ignored.

- x:

  A `usdt_reliability` object.

## Value

An object of class `usdt_reliability` containing:

- `$tasks`: Overall reliability and variance components for each task.

- `$subjects`: Participant-level \\d'\\, error variances, and individual
  reliabilities.

- Bootstrap intervals for both components when available in `object`.

For a Bayesian fit, both tables hold posterior means and add the columns
`sd`, `conf.low` and `conf.high`, and `$posterior$tasks` keeps the draws
of each task's reliability.

## Details

For participant \\i\\ in task \\j\\, reliability is defined as:
\$\$\frac{\tau_j^2}{\tau_j^2 + v\_{ij}}\$\$ where \\\tau_j^2\\ is the
true variance in sensitivity across participants from the model's random
effects, and \\v\_{ij}\\ is the squared standard error of \\d'\\ for
that participant. This error variance reflects how precisely their
trials determine sensitivity, accounting for trial count, performance
level on the probit curve, and uncertainty in the criterion.

This variance is calculated using the large-sample Fisher information
formula from Gourevitch and Galanter (1967). Unlike
[`sdt_moments()`](https://ricardoreysaez.github.io/uSDT/reference/sdt_moments.md),
which evaluates that formula at raw empirical proportions (`var_gg`),
`usdt_reliability()` evaluates it at the response probabilities
predicted by the fitted hierarchical model.

Overall task reliability averages \\v\_{ij}\\ across participants,
representing the expected proportion of true variance for a participant
drawn at random from the sample.

When `object` includes bootstrap replicates from
[`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_boot.md),
confidence intervals for reliability are computed automatically across
all retained samples.

For a Bayesian fit, reliability is computed draw by draw: each draw
gives \\\tau_j^2\\ and, at that draw's subject sensitivities and
criteria, each \\v\_{ij}\\. The tables report posterior means, posterior
standard deviations and central credible intervals. Under unequal
variances the information of the signal trials is scaled by the signal
standard deviation \\\sigma_s\\, so that \\v\_{ij} = \sigma_s^2 /
w\_{\mathrm{signal}} + 1 / w\_{\mathrm{noise}}\\ with an estimated
criterion, which reduces to the Gourevitch and Galanter formula when
\\\sigma_s = 1\\.

## References

Gourevitch, V., & Galanter, E. (1967). A significance test for one
parameter isosensitivity functions. *Psychometrika*, 32(1), 25–33.
[doi:10.1007/BF02289402](https://doi.org/10.1007/BF02289402)

## See also

[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md),
[`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_boot.md),
[`sdt_moments()`](https://ricardoreysaez.github.io/uSDT/reference/sdt_moments.md)

## Examples

``` r
# Contextual cuing data from Vadillo et al. (2025)
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

r <- usdt_reliability(hsdt(d))
r
#> ── Reliability summary ───────────────────────────────────────────────────────── 
#> 
#>   Task         Subjects Group-level estimate   By-subject estimate median [min, max]
#>   Awareness         104                 .129   .130 [.119, .131]
#>   Influence         104                 .323   .323 [.322, .323]

# Participant-level estimates (one row per subject and task)
head(r$subjects)
#>        task subj    dprime   variance reliability
#> 1 Awareness 2001 0.2707158 0.09899329   0.1305303
#> 2 Awareness 2002 0.2485344 0.10427183   0.1247468
#> 3 Awareness 2003 0.2172227 0.10017561   0.1291888
#> 4 Awareness 2004 0.3175901 0.09908306   0.1304275
#> 5 Awareness 2005 0.1444059 0.10037191   0.1289687
#> 6 Awareness 2006 0.1733233 0.09926175   0.1302232
```
