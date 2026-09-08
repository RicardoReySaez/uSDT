# Fit a hierarchical signal detection theory model

Fits the model to prepared data and tests the three hypotheses of the
unconscious processing design. Every subject has one sensitivity in each
task, and the model allows the two sensitivities to correlate across
subjects. Estimation uses
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html).

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
  [`usdt_data_long()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_data.md)
  or
  [`usdt_data_tasks()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_data.md).

- estimation:

  Estimation method. Only `"frequentist"` is available, and it fits the
  model by maximum likelihood.

- fix_criteria:

  What to do with the criteria. `"auto"` fixes to zero every criterion
  that the data make zero by construction, which happens after a median
  split under deviation coding. `"none"` estimates all of them.

- level:

  Confidence level.

- optimizer:

  Optimizer given to
  [`lme4::glmerControl()`](https://rdrr.io/pkg/lme4/man/lmerControl.html).
  The function tries other optimizers when this one does not converge or
  reaches a singular fit.

- ...:

  Named arguments passed on to every
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) call. The
  formula, the data, the family, the optimizer and `nAGQ = 1` stay under
  the control of the package. The print and summary methods ignore this
  argument.

- object:

  An `hsdt` object.

- x:

  An `hsdt` object.

## Value

An object of class `hsdt`. It holds the fitted model in `fit`, the three
hypothesis tests in `tests`, the estimates they are built from in
`pars`, a description of the model formula in `design`, and the fitting
diagnostics in `diagnostics`. Use
[`summary()`](https://rdrr.io/r/base/summary.html) to read it.

## Details

The model works on counts of signal responses. It uses a probit link, so
its parameters keep the usual signal detection meaning. The fixed
effects give the average sensitivity of each task, and any criterion
that is estimated. The random effects give the departure of each subject
from those averages, and the two sensitivities share one covariance,
which is what makes the latent correlation available.

The three hypotheses come out of that covariance and the two averages.
H1 is the difference between the average sensitivities. H2 is their
correlation across subjects. H3 is the regression of the indirect
sensitivity on the direct one, and its intercept is the sensitivity
expected in the indirect task from a subject whose direct sensitivity is
zero. See
[`usdt_tests()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_hypotheses.md).

Some datasets do not carry enough information for H2 and H3. The
function still returns the fit and warns, and
[`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_boot.md)
can then provide intervals by simulation.

## See also

[`usdt_data_long()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_data.md),
[`usdt_tests()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_hypotheses.md),
[`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_boot.md),
[`plot.hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/plot.hsdt.md)

## Examples

``` r
# \donttest{
set.seed(1)
df <- usdt_simulate(n_subj = 40, n_trials = 100)
d  <- usdt_data_long(df, task_col = "task",
                     task_levels      = c(direct = "D", indirect = "I"),
                     subject_col      = "subj",
                     condition_col    = "cond",
                     condition_levels = c(signal = 1, noise = 0),
                     response_col     = "response",
                     response_levels  = c(signal = 1, noise = 0))
m <- hsdt(d)
summary(m)
#> ── Model summary ─────────────────────────────────────────────────────────────── 
#> 
#>   Subjects:       40
#>   Observations:   160 aggregated rows (8,000 trials)
#>   Family:         binomial (probit)
#>   Coding:         deviation
#>   Criteria:       Direct estimated, Indirect estimated
#>   Estimation:     lme4::glmer (bobyqa)
#>   Convergence:    TRUE
#> 
#> ── Fixed effects ───────────────────────────────────────────────────────────────
#> 
#>   Parameter      Task       Estimate       SE  95% CI                   z   p-value
#>   criterion      Direct      -0.0285   0.0399  [ -0.107,  0.050]    -0.71      .475
#>   criterion      Indirect     0.0840   0.0479  [ -0.010,  0.178]     1.75      .080
#>   d'             Direct       0.8201   0.0869  [  0.650,  0.991]     9.43     <.001
#>   d'             Indirect     0.3528   0.0646  [  0.226,  0.479]     5.46     <.001
#> 
#> ── Random effects ──────────────────────────────────────────────────────────────
#> 
#>   Parameter      Task       Estimate       SE  95% CI            
#>   sd(criterion)  Direct       0.2150   0.0336  [  0.158,  0.292]
#>   sd(criterion)  Indirect     0.2743   0.0378  [  0.209,  0.359]
#>   cor(c)         both         0.0710   0.2041  [ -0.319,  0.441]
#>   sd(d')         Direct       0.4816   0.0702  [  0.362,  0.641]
#>   sd(d')         Indirect     0.3165   0.0600  [  0.218,  0.459]
#>   cor(d')        both         0.4680   0.2007  [  0.004,  0.766]
#> 
#> ── Hypotheses ──────────────────────────────────────────────────────────────────
#> 
#> H1: Group-level sensitivity difference (Δd' = Direct d' - Indirect d')
#>   Parameter       Estimate       SE  95% CI                   z   p-value
#>   Δd' (D - I)      0.4673   0.0903  [  0.290,  0.644]     5.17     <.001
#> 
#> H2: Correlation between sensitivities across tasks
#>   Parameter       Estimate       SE  95% CI                   z   p-value
#>   rho               0.4680   0.2007  [  0.004,  0.766]     1.93      .054
#> 
#> H3: Latent regression of Indirect d' on Direct d'
#>   Parameter       Estimate       SE  95% CI                   z   p-value
#>   Intercept         0.1006   0.1344  [ -0.163,  0.364]     0.75      .454
#>   Slope             0.3075   0.1460  [  0.021,  0.594]     1.93      .054
#> 
#> ── Notes ───────────────────────────────────────────────────────────────────────
#> 
#>   Use usdt_boot(fit, nsim = 1000, ncores = 4) for bootstrap CIs.
# }
```
