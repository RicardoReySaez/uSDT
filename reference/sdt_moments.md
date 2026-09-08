# Signal detection measures for each subject

Computes empirical hit rates, false-alarm rates, sensitivity (\\d'\\),
and response criteria for each participant without fitting a model. It
can also calculate sampling variances, standard errors, and expected
values for \\d'\\.

## Usage

``` r
sdt_moments(
  data,
  subject_col = NULL,
  condition_col = NULL,
  condition_levels = NULL,
  response_col = NULL,
  response_levels = NULL,
  coding = c("deviation", "treatment"),
  correction = c("hautus", "none"),
  variances = FALSE
)
```

## Arguments

- data:

  A `usdt_data` object or a standard trial-level data frame. When given
  a `usdt_data` object, the function processes both tasks and includes a
  `task` column in the output.

- subject_col, condition_col, response_col:

  Column names for subject, condition, and response variables. Only
  required when `data` is a plain data frame.

- condition_levels, response_levels:

  Named vectors mapping condition and response labels, like
  `c(signal = "old", noise = "new")`. Required for a plain data frame. A
  `usdt_data` object supplies its own roles and needs neither.

- coding:

  Criterion definition to report: `"deviation"` measures the criterion
  from the midpoint between the signal and noise distributions, whereas
  `"treatment"` measures it from the noise distribution. A `usdt_data`
  object supplies its own coding.

- correction:

  Handling of extreme rates (0 or 1) that make \\d'\\ infinite.
  `"hautus"` adds 0.5 to all four cell counts for affected participants.
  `"none"` leaves infinite values in place.

- variances:

  Logical. If `TRUE`, computes the sampling variance of \\d'\\ from
  Gourevitch and Galanter (1967) and Miller (1996), each with its own
  standard error, as well as the expected value of \\d'\\ under Miller's
  distribution.

## Value

A data frame with one row per subject (or per subject and task for
`usdt_data` inputs) containing:

- `hit`, `miss`, `fa`, `cr`: Raw response counts.

- `hr`, `far`: Observed hit and false-alarm rates.

- `zhr`, `zfar`: Probit-transformed rates.

- `dprime`, `criterion`: Descriptive SDT estimates.

- `corrected`: Logical flag indicating whether the participant received
  an edge correction.

- `var_gg`, `se_gg`: Asymptotic variance and standard error from
  Gourevitch and Galanter (1967), present when `variances = TRUE`.

- `var_miller`, `se_miller`, `expected_dprime`: Moments from Miller
  (1996), present when `variances = TRUE`.

## Details

Edge corrections apply only to participants with extreme rates (0 or 1)
rather than the whole sample, leaving well-defined rates unchanged.

When requested, the sampling variance of \\d'\\ is estimated using the
asymptotic approximation of Gourevitch and Galanter (1967) and the
binomial-distribution approach of Miller (1996). See Suero et al. (2017)
for a comparison between the two approaches.

## References

Gourevitch, V., & Galanter, E. (1967). A significance test for one
parameter isosensitivity functions. *Psychometrika*, 32(1), 25–33.
[doi:10.1007/BF02289402](https://doi.org/10.1007/BF02289402)

Hautus, M. J. (1995). Corrections for extreme proportions and their
biasing effects on estimated values of \\d'\\. *Behavior Research
Methods, Instruments, & Computers*, 27(1), 46–51.
[doi:10.3758/BF03203619](https://doi.org/10.3758/BF03203619)

Miller, J. (1996). The sampling distribution of \\d'\\. *Perception &
Psychophysics*, 58(1), 65–72.
[doi:10.3758/BF03205476](https://doi.org/10.3758/BF03205476)

Suero, M., Privado, J., & Botella, J. (2017). Methods to estimate the
variance of some indices of the signal detection theory: A simulation
study. *Psicologica*, 38(1), 77–109.

## See also

[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md)

## Examples

``` r
# 1. From a prepared usdt_data object (both tasks at once)
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

head(sdt_moments(d))
#>     task subj hit miss fa cr      hr     far        zhr        zfar      dprime
#> 1 Direct 2001  19   13 15 17 0.59375 0.46875 0.23720211 -0.07841241  0.31561452
#> 2 Direct 2002  24    8 21 11 0.75000 0.65625 0.67448975  0.40225007  0.27223968
#> 3 Direct 2003  20   12 19 13 0.62500 0.59375 0.31863936  0.23720211  0.08143725
#> 4 Direct 2004  22   10 10 22 0.68750 0.31250 0.48877641 -0.48877641  0.97755282
#> 5 Direct 2005  17   15 23  9 0.53125 0.71875 0.07841241  0.57913216 -0.50071975
#> 6 Direct 2006  17   15 20 12 0.53125 0.62500 0.07841241  0.31863936 -0.24022695
#>     criterion corrected
#> 1 -0.07939485     FALSE
#> 2 -0.53836991     FALSE
#> 3 -0.27792074     FALSE
#> 4  0.00000000     FALSE
#> 5 -0.32877229     FALSE
#> 6 -0.19852589     FALSE

# 2. From raw trials with sampling variances and standard errors
head(sdt_moments(vadillo_awareness,
                 subject_col      = "subj",
                 condition_col    = "condition",
                 condition_levels = c(signal = "old", noise = "new"),
                 response_col     = "judged.old",
                 response_levels  = c(signal = 1, noise = 0),
                 variances        = TRUE))
#>   subj hit miss fa cr      hr     far        zhr        zfar      dprime
#> 1 2001  19   13 15 17 0.59375 0.46875 0.23720211 -0.07841241  0.31561452
#> 2 2002  24    8 21 11 0.75000 0.65625 0.67448975  0.40225007  0.27223968
#> 3 2003  20   12 19 13 0.62500 0.59375 0.31863936  0.23720211  0.08143725
#> 4 2004  22   10 10 22 0.68750 0.31250 0.48877641 -0.48877641  0.97755282
#> 5 2005  17   15 23  9 0.53125 0.71875 0.07841241  0.57913216 -0.50071975
#> 6 2006  17   15 20 12 0.53125 0.62500 0.07841241  0.31863936 -0.24022695
#>     criterion corrected     var_gg     se_gg var_miller se_miller
#> 1 -0.07939485     FALSE 0.09930004 0.3151191  0.1049949 0.3240291
#> 2 -0.53836991     FALSE 0.11009703 0.3318087  0.1207283 0.3474598
#> 3 -0.27792074     FALSE 0.10104010 0.3178681  0.1073747 0.3276808
#> 4  0.00000000     FALSE 0.10713629 0.3273168  0.1160479 0.3406581
#> 5 -0.32877229     FALSE 0.10470577 0.3235827  0.1128152 0.3358798
#> 6 -0.19852589     FALSE 0.10013446 0.3164403  0.1061461 0.3258007
#>   expected_dprime
#> 1      0.32409622
#> 2      0.28275954
#> 3      0.08381293
#> 4      1.00624073
#> 5     -0.51642519
#> 6     -0.24693939
```
