# Fit a hierarchical signal detection theory model

Models paired awareness and influence sensitivities and criteria, their
population means and between-subject variation. Frequentist estimation
uses [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html).
Bayesian estimation uses brms and Stan.

## Usage

``` r
hsdt(
  data,
  estimation = c("frequentist", "bayesian"),
  fix_criteria = c("auto", "none"),
  level = 0.95,
  optimizer = "bobyqa",
  priors = usdt_priors(),
  unequal_variances = FALSE,
  latent = c("correlation", "regression"),
  backend = getOption("brms.backend", "rstan"),
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

  Estimation framework: `"frequentist"` (maximum likelihood via Laplace
  approximation with lme4) or `"bayesian"` (Hamiltonian Monte Carlo with
  Stan, through brms).

- fix_criteria:

  How to handle response criteria. `"auto"` fixes to zero a dichotomized
  task's criterion when its mean absolute descriptive subject criterion
  is below 0.02. `"none"` estimates all criteria.

- level:

  Confidence level for Wald intervals, or credible level for the
  posterior intervals of a Bayesian fit (default is 0.95).

- optimizer:

  Primary optimizer passed to
  [`lme4::glmerControl()`](https://rdrr.io/pkg/lme4/man/lmerControl.html).
  Alternative optimizers are automatically evaluated if the default
  fails to converge or produces a singular fit. Frequentist estimation
  only.

- priors:

  A `usdt_priors` object from
  [`usdt_priors()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_priors.md).
  Bayesian estimation only.

- unequal_variances:

  Logical. Estimate the standard deviation of the signal distribution of
  each task instead of fixing it to that of the noise distribution.
  Bayesian estimation only.

- latent:

  How the Bayesian model writes the relation between the two latent
  sensitivities, which decides where it places the prior on the location
  of the influence sensitivity: on its mean, by their latent correlation
  (`"correlation"`, the default), or on the intercept of the latent
  regression of the influence on the awareness sensitivity
  (`"regression"`). Both keep the priors on the awareness mean, the
  standard deviations and the correlation. Bayesian estimation only.

- backend:

  Stan interface brms uses for Bayesian estimation: `"rstan"` (available
  from CRAN) or `"cmdstanr"` (installed from the Stan R-universe
  together with CmdStan). As in brms, the default is the `brms.backend`
  option, or `"rstan"` when it is not set. Bayesian estimation only.

- ...:

  For frequentist estimation, additional arguments passed to
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html). The model
  formula, family and data remain managed by the package. For Bayesian
  estimation, sampling options: `chains` (default 4), `iter` (iterations
  per chain including warmup, default 6000), `warmup` (default 1000, so
  each chain keeps 5000 draws), `cores` (chains run at once, default
  `min(parallel::detectCores() - 2, 4)`, at least 1), `seed`,
  `control = list(adapt_delta, max_treedepth)` (default 0.95 and 10) and
  `refresh` (iterations between progress updates of each chain, default
  1% of `iter`, and `0` hides the progress bar).

- object:

  An `hsdt` object.

- x:

  An `hsdt` object.

## Value

An object of class `hsdt` containing:

- `$fit`: The underlying fit: a `glmerMod` object from lme4, or a
  `brmsfit` object from brms (see
  [`brms::prior_summary()`](https://mc-stan.org/rstantools/reference/prior_summary.html)
  and
  [`brms::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html)).

- `$tests`: Summary table for hypotheses H1, H2, and H3.

- `$design`: Summary of the model specification.

- `$diagnostics`: Convergence diagnostics.

A frequentist fit also contains `$pars`, the parameter estimates and
their joint covariance. A Bayesian fit also contains `$priors`. Its
posterior draws are those of `$fit`, where `b_dA_Intercept` and
`b_dI_Intercept` are the mean \\d'\\ of the awareness and influence
tasks, `sd_subj__dA_Intercept` and `sd_subj__dI_Intercept` their SDs,
`cor_subj__dA_Intercept__dI_Intercept` their correlation,
`b_cA_Intercept` and `b_cI_Intercept` the mean of an estimated criterion
(on the classical scale), and `b_sA_Intercept` and `b_sI_Intercept` the
signal SDs under unequal variances.

## Details

Hypotheses evaluated by default:

- **H1:** Mean sensitivity difference between tasks.

- **H2:** Latent correlation of sensitivities across participants.

- **H3:** Latent regression of influence on awareness sensitivity. Its
  intercept reflects expected influence performance when awareness
  sensitivity is zero (\\d'\_{\mathrm{Awareness}} = 0\\).

When sample sizes or trial counts are low, variance components can reach
singular boundaries. In these cases, the function issues a warning, and
parametric bootstrap intervals can be calculated using
[`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_boot.md),
or the model can be fitted with `estimation = "bayesian"`.

## Bayesian estimation

The Bayesian model places its priors on the mean, standard deviation and
correlation of the two sensitivities and on the criteria. See
[`usdt_priors()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_priors.md).
The slope and intercept of the latent regression are derived from them
draw by draw. The hypothesis table reports the posterior mean and
standard deviation, a central credible interval and a two-sided
posterior p-value, \\2 \min\\P(\theta \> 0 \mid y), P(\theta \< 0 \mid
y)\\\\. The p-value falls below `1 - level` exactly when zero lies
outside the credible interval, so with the default `level = 0.95` a
p-value below .05 means that the 95% credible interval excludes zero.

With `latent = "regression"`, the prior of the influence mean moves to
the intercept \\\beta_0\\ of the latent regression, and the influence
mean follows as \\\mu_I = \beta_0 + \rho \\ (\sigma_I / \sigma_A) \\
\mu_A\\. The intercept tested by H3 then has its specified prior, and
the correlation (H2) keeps its own. The difference between the means
(H1) has an induced, heavy-tailed prior instead, and because the data
fix \\\mu_I\\, the intercept prior also informs the posterior of
\\\rho\\. See
[`usdt_priors()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_priors.md)
and
[`plot.usdt_priors()`](https://ricardoreysaez.github.io/uSDT/reference/plot.usdt_priors.md).
The fit samples the same parameters and reports the same hypotheses
either way. Near-zero awareness SDs can produce difficult posterior
geometry and divergences. A prior with less weight near zero may help,
if substantively justified. Inspect its implications with
[`plot.usdt_priors()`](https://ricardoreysaez.github.io/uSDT/reference/plot.usdt_priors.md).

With `unequal_variances = TRUE`, the noise distribution of each task is
\\N(-d/2, 1)\\ and the signal distribution \\N(d/2, \sigma_s^2)\\, so
\\d\\ is the distance between their means in units of the noise standard
deviation, and H1-H3 compare these \\d\\. When \\\sigma_s\\ differs
between the tasks, \\d\\ no longer measures the same discriminability in
both: the measure comparable with an equal-variance \\d'\\ is \\d_a = d
\sqrt{2 / (1 + \sigma_s^2)}\\, and the mean \\d_a\\ of a task is its
mean \\d\\ times that factor, which can be computed from the draws of
`b_dA_Intercept` and `b_sA_Intercept` (and their influence counterparts)
in `$fit`. The correlation (H2) is the same on either scale.

The model is a brms model of the aggregated counts, with one prior per
parameter and task. It is compiled the first time it is used with a
given set of priors, which needs a C++ toolchain, and stays in a
per-user cache directory for later sessions. Inspect the reported
divergences, tree-depth warnings, R-hat and effective sample sizes
before interpreting results.

## See also

[`usdt_data_tasks()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_data.md),
[`usdt_priors()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_priors.md),
[`usdt_tests()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_hypotheses.md),
[`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_boot.md),
[`plot.hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/plot.hsdt.md)

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

m <- hsdt(d)

# Full summary table with SDT parameters and hypothesis tests
summary(m)
#> ── Model summary ─────────────────────────────────────────────────────────────── 
#> 
#>   Subjects:       104
#>   Observations:   416 aggregated rows (46,592 trials)
#>   Family:         binomial (probit)
#>   Coding:         deviation
#>   Criteria:       Awareness estimated, Influence fixed to 0 (Meyen split, mean |c| = 0.0000)
#>   Estimation:     lme4::glmer (Nelder_Mead), tried bobyqa -> Nelder_Mead
#>   Convergence:    TRUE
#> 
#> ── Fixed effects ───────────────────────────────────────────────────────────────
#> 
#>   Parameter      Task       Estimate       SE  95% CI                   z   p-value
#>   criterion      Awareness   -0.0357   0.0289  [ -0.092,  0.021]    -1.24      .215
#>   d'             Awareness    0.2354   0.0335  [  0.170,  0.301]     7.03     <.001
#>   d'             Influence    0.1283   0.0153  [  0.098,  0.158]     8.41     <.001
#> 
#> ── Random effects ──────────────────────────────────────────────────────────────
#> 
#>   Parameter      Task       Estimate       SE  95% CI            
#>   sd(criterion)  Awareness    0.2473   0.0246  [  0.204,  0.301]
#>   sd(d')         Awareness    0.1219   0.0666  [  0.042,  0.356]
#>   sd(d')         Influence    0.0883   0.0190  [  0.058,  0.135]
#>   cor(d')        both         0.4910   0.5188  [ -0.665,  0.954]
#> 
#> ── Hypotheses ──────────────────────────────────────────────────────────────────
#> 
#> H1: Group-level sensitivity difference (Δd' = Influence d' - Awareness d')
#>   Parameter       Estimate       SE  95% CI                   z   p-value
#>   Δd' (I - A)     -0.1071   0.0354  [ -0.176, -0.038]    -3.03      .002
#> 
#> H2: Correlation between sensitivities across tasks
#>   Parameter       Estimate       SE  95% CI                   z   p-value
#>   rho               0.4910   0.5188  [ -0.665,  0.954]     1.01      .313
#> 
#> H3: Latent regression of Influence d' on Awareness d'
#>   Parameter       Estimate       SE  95% CI                   z   p-value
#>   Intercept         0.0445   0.1159  [ -0.183,  0.272]     0.38      .701
#>   Slope             0.3558   0.4865  [ -0.598,  1.309]     1.01      .313
#> 
#> ── Notes ───────────────────────────────────────────────────────────────────────
#> 
#>   Use usdt_boot(fit, nsim = 1000, ncores = 2) for bootstrap CIs.

# Inspect the model formula (influence criterion omitted by default)
m$design$formula
#> cbind(y, n - y) ~ 0 + c_A + d_A + d_I + (0 + c_A | subj) + (0 + 
#>     d_A + d_I | subj)
#> <environment: 0x55799ea17fd8>

if (FALSE) { # \dontrun{
# Requires a C++ toolchain and the BH and RcppEigen headers for rstan.
# The first call compiles the model.
mb <- hsdt(d, estimation = "bayesian", cores = 2, seed = 1)
summary(mb)
} # }
```
