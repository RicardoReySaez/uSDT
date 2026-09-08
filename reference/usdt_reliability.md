# Reliability of direct and indirect task measures

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
  which may also contain bootstrap results from
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
# \donttest{
# Contextual cuing data from Vadillo et al. (2025)
d <- usdt_data_tasks(
  direct   = vadillo_awareness,
  indirect = vadillo_cuing,
  subject_col      = "subj",
  condition_col    = "condition",
  condition_levels = c(signal = "old", noise = "new"),
  response_col     = list(direct = "judged.old", indirect = "rt"),
  response_levels  = list(direct   = c(signal = 1, noise = 0),
                          indirect = c(signal = "faster", noise = "slower")),
  dichotomize      = list(direct = FALSE, indirect = TRUE)
)

r <- usdt_reliability(hsdt(d))
r
#> ── Reliability summary ───────────────────────────────────────────────────────── 
#> 
#>   Task         Subjects Group-level estimate   By-subject estimate median [min, max]
#>   Direct            104                 .129   .130 [.119, .131]
#>   Indirect          104                 .323   .323 [.322, .323]

# Participant-level estimates (one row per subject and task)
head(r$subjects)
#>     task subj    dprime   variance reliability
#> 1 Direct 2001 0.2706988 0.09899319   0.1304247
#> 2 Direct 2002 0.2485231 0.10427127   0.1246457
#> 3 Direct 2003 0.2172299 0.10017551   0.1290841
#> 4 Direct 2004 0.3175052 0.09908257   0.1303224
#> 5 Direct 2005 0.1444723 0.10037192   0.1288640
#> 6 Direct 2006 0.1733643 0.09926182   0.1301177
# }
```
