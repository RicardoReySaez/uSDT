# Fit a hierarchical signal detection theory model

Fits a bivariate hierarchical SDT model using
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) and evaluates
the core unconscious processing hypotheses. The model estimates
task-specific sensitivities (\\d'\\) and response criteria (\\c\\) as
fixed effects, while estimating their variation and correlation across
participants via random effects.

## Usage

``` r
hsdt(
  data,
  estimation = c("frequentist"),
  fix_criteria = c("auto", "none"),
  level = 0.95,
  optimizer = "bobyqa",
  ...
)

# S3 method for class 'hsdt'
summary(object, ...)

# S3 method for class 'hsdt'
print(x, ...)
```

## Arguments

- data:

  A `usdt_data` object from
  [`usdt_data_tasks()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_data.md)
  or
  [`usdt_data_long()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_data.md).

- estimation:

  Estimation framework. Currently only `"frequentist"` (maximum
  likelihood via Laplace approximation) is supported.

- fix_criteria:

  How to handle response criteria. `"auto"` fixes to zero any criterion
  that is zero by design (such as a task split at the median under
  deviation coding). `"none"` estimates all criteria.

- level:

  Confidence level for Wald intervals (default is 0.95).

- optimizer:

  Primary optimizer passed to
  [`lme4::glmerControl()`](https://rdrr.io/pkg/lme4/man/lmerControl.html).
  Alternative optimizers are automatically evaluated if the default
  fails to converge or produces a singular fit.

- ...:

  Additional arguments passed to
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html). Model
  formula, family, and data inputs remain managed by the package.

- object:

  An `hsdt` object.

- x:

  An `hsdt` object.

## Value

An object of class `hsdt` containing:

- `$fit`: The underlying `glmerMod` object from `lme4`.

- `$tests`: Summary table for hypotheses H1, H2, and H3.

- `$pars`: Model parameter estimates on the SDT scale.

- `$design`: Summary of the model specification and formula.

- `$diagnostics`: Convergence flags and singular fit indicators.

## Details

The model fits trial counts with a binomial probit link, directly
mapping coefficients to standard Signal Detection Theory parameters.
Fixed effects capture population sensitivities and criteria, while
random effects estimate participant variation and the latent correlation
between direct and indirect sensitivity.

Hypotheses evaluated by default:

- **H1:** Mean sensitivity difference between tasks.

- **H2:** Latent correlation of sensitivities across participants.

- **H3:** Latent regression of indirect on direct sensitivity. Its
  intercept reflects expected indirect performance when direct awareness
  is zero (\\d'\_{\mathrm{Direct}} = 0\\).

When sample sizes or trial counts are low, variance components can reach
singular boundaries. In these cases, the function issues a warning, and
parametric bootstrap intervals can be calculated using
[`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_boot.md).

## See also

[`usdt_data_tasks()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_data.md),
[`usdt_tests()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_hypotheses.md),
[`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_boot.md),
[`plot.hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/plot.hsdt.md)

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

m <- hsdt(d)

# Full summary table with SDT parameters and hypothesis tests
summary(m)
#> ── Model summary ─────────────────────────────────────────────────────────────── 
#> 
#>   Subjects:       104
#>   Observations:   416 aggregated rows (46,592 trials)
#>   Family:         binomial (probit)
#>   Coding:         deviation
#>   Criteria:       Direct estimated, Indirect fixed to 0 (Meyen split, mean |c| = 0.0000)
#>   Estimation:     lme4::glmer (bobyqa)
#>   Convergence:    TRUE
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
#>   cor(d')        both         0.4912   0.5190  [ -0.666,  0.954]
#> 
#> ── Hypotheses ──────────────────────────────────────────────────────────────────
#> 
#> H1: Group-level sensitivity difference (Δd' = Indirect d' - Direct d')
#>   Parameter       Estimate       SE  95% CI                   z   p-value
#>   Δd' (I - D)     -0.1071   0.0354  [ -0.176, -0.038]    -3.03      .002
#> 
#> H2: Correlation between sensitivities across tasks
#>   Parameter       Estimate       SE  95% CI                   z   p-value
#>   rho               0.4912   0.5190  [ -0.666,  0.954]     1.01      .313
#> 
#> H3: Latent regression of Indirect d' on Direct d'
#>   Parameter       Estimate       SE  95% CI                   z   p-value
#>   Intercept         0.0445   0.1160  [ -0.183,  0.272]     0.38      .701
#>   Slope             0.3561   0.4870  [ -0.599,  1.311]     1.01      .313
#> 
#> ── Notes ───────────────────────────────────────────────────────────────────────
#> 
#>   Use usdt_boot(fit, nsim = 1000, ncores = 4) for bootstrap CIs.

# Inspect the model formula (indirect criterion omitted by default)
m$design$formula
#> cbind(y, n - y) ~ 0 + c_D + d_D + d_I + (0 + c_D | subj) + (0 + 
#>     d_D + d_I | subj)
#> <environment: 0x557b9b17f2d8>
# }
```
