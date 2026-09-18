# Parametric bootstrap intervals for hierarchical SDT models

Simulates new datasets from a fitted model using
[`lme4::bootMer()`](https://rdrr.io/pkg/lme4/man/bootMer.html), refits
the model to each replicate, and computes bootstrap confidence intervals
for the three core hypotheses (H1, H2, H3). This is especially useful
when asymptotic Wald intervals are unreliable due to singular or
boundary fits.

## Usage

``` r
usdt_boot(
  object,
  nsim = 1000,
  ncores = 1L,
  max_attempts = 2 * nsim,
  seed = NULL,
  progress = interactive(),
  level = 0.95,
  type = c("perc", "norm", "basic")
)
```

## Arguments

- object:

  An `hsdt` object fitted by
  [`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md).

- nsim:

  Target number of successful replicates (at least 500).

- ncores:

  Number of CPU cores for parallel processing. Values above 1 create a
  temporary cluster that works across Windows, macOS, and Linux, and
  automatically stops when finished.

- max_attempts:

  Maximum number of refits to attempt. Defaults to `2 * nsim`.

- seed:

  Random seed for reproducibility. With the default `NULL`, no seed is
  set and the bootstrap continues the current random number stream.

- progress:

  Logical. Display a progress bar during fitting (defaults to `TRUE` in
  interactive sessions).

- level:

  Confidence level for intervals (default is 0.95).

- type:

  Type of bootstrap interval: `"perc"` (percentile), `"norm"` (normal
  approximation with bias correction), or `"basic"` (empirical basic).
  For correlations, `"norm"` and `"basic"` operate on the Fisher-\\z\\
  scale; if the sample correlation lies on the boundary (-1 or 1), these
  types return `NA`, whereas `"perc"` remains available.

## Value

An updated `hsdt` object where interval columns in `$tests` are replaced
by bootstrap estimates. A new `$boot` element contains:

- `$t`: Matrix of replicates for the three hypotheses.

- `$variance`: Replicates of sensitivity variances.

- `$population`: Replicates of average criteria and sensitivities.

- `$subjects`: Replicates of individual-level parameters.

- Run diagnostics and convergence counts (`usable`, `attempted`,
  `retained`, `failures`).

## Details

Replicates are dropped only if the model fails to fit or does not
converge. Singular fits and boundary estimates are intentionally
retained because discarding them artificially narrows intervals in
constrained settings.

When an attempted run finishes with fewer than 500 usable replicates,
the original Wald intervals are preserved, a warning is issued, and raw
attempt diagnostics are stored in `$boot`.

Point estimates remain identical to the original model fit. Two-sided
bootstrap \\p\\-values compare the observed test statistic against the
centered bootstrap distribution using standard finite-sample adjustment
(\\(k + 1) / (B + 1)\\), ensuring \\p\\-values never equal zero.

## See also

[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md),
[`usdt_tests()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_hypotheses.md)

## Examples

``` r
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

m <- hsdt(d)

# \donttest{
# Run parametric bootstrap with 500 replicates. Refitting this model 500
# times takes several minutes
b <- usdt_boot(m, nsim = 500, seed = 1)

# Inspect updated summary with bootstrap intervals and p-values
summary(b)
#> ── Model summary ─────────────────────────────────────────────────────────────── 
#> 
#>   Subjects:       104
#>   Observations:   416 aggregated rows (46,592 trials)
#>   Family:         binomial (probit)
#>   Coding:         deviation
#>   Criteria:       Direct estimated, Indirect fixed to 0 (Meyen split, mean |c| = 0.0000)
#>   Estimation:     lme4::glmer (bobyqa)
#>   Convergence:    TRUE
#>   Bootstrap:      500 usable replicates (508 attempts; 190 at the boundary)
#> 
#> ── Fixed effects ───────────────────────────────────────────────────────────────
#> 
#>   Parameter      Task       Estimate       SE  95% CI                   z   p-value
#>   criterion      Direct      -0.0357   0.0289  [ -0.092,  0.021]    -1.24      .215
#>   d'             Direct       0.2354   0.0335  [  0.170,  0.301]     7.03     <.001
#>   d'             Indirect     0.1283   0.0153  [  0.098,  0.158]     8.41     <.001
#> 
#> ── Random effects ──────────────────────────────────────────────────────────────
#> 
#>   Parameter      Task       Estimate       SE  95% CI            
#>   sd(criterion)  Direct       0.2473   0.0246  [  0.203,  0.301]
#>   sd(d')         Direct       0.1219   0.0666  [  0.042,  0.356]
#>   sd(d')         Indirect     0.0883   0.0190  [  0.058,  0.135]
#>   cor(d')        both         0.4912   0.5202  [ -1.000,  1.000]
#> 
#> ── Hypotheses ──────────────────────────────────────────────────────────────────
#> 
#> H1: Group-level sensitivity difference (Δd' = Indirect d' - Direct d')
#>   Parameter         Estimate  Boot SE  Boot 95% CI (percentile)  Boot p-value
#>   Δd' (I - D)       -0.1071   0.0346  [ -0.174, -0.038]                 .004
#> 
#> H2: Correlation between sensitivities across tasks
#>   Parameter         Estimate  Boot SE  Boot 95% CI (percentile)  Boot p-value
#>   rho                 0.4912   0.5202  [ -1.000,  1.000]                 .483
#> 
#> H3: Latent regression of Indirect d' on Direct d'
#>   Parameter         Estimate  Boot SE  Boot 95% CI (percentile)  Boot p-value
#>   Intercept           0.0445   0.8949  [ -1.125,  0.670]                 .685
#>   Slope               0.3561   3.9094  [ -2.224,  5.438]                 .441

# Check fit diagnostics across bootstrap replicates
b$boot[c("usable", "attempted", "retained", "failures")]
#> $usable
#> [1] 500
#> 
#> $attempted
#> [1] 508
#> 
#> $retained
#> singular boundary 
#>      187      190 
#> 
#> $failures
#>    non_finite non_converged 
#>             7             1 
#> 
# }
```
