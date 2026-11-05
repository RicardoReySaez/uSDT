# Test the three core hypotheses of a hierarchical SDT model

Evaluates the difference between average task sensitivities (H1), their
correlation across subjects (H2), and the regression of influence
sensitivity on awareness sensitivity (H3). `usdt_tests()` computes all
three together, while individual functions compute them separately.

## Usage

``` r
usdt_tests(fit, awareness = "d_A", influence = "d_I", level = 0.95)

sensitivity_diff(fit, awareness = "d_A", influence = "d_I", level = 0.95)

latent_cor(fit, awareness = "d_A", influence = "d_I", level = 0.95)

latent_regression(fit, awareness = "d_A", influence = "d_I", level = 0.95)
```

## Arguments

- fit:

  A fitted model: an `hsdt` object, frequentist or Bayesian, a
  `glmerMod` from
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) or a
  `brmsfit` from
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).

- awareness, influence:

  Character strings naming the sensitivity terms in the model. Defaults
  match the internal naming of
  [`hsdt()`](https://ricardoreysaez.github.io/uSDT/dev/reference/hsdt.md).
  For custom models, both terms must be fixed effects and share a common
  random-effects grouping by subject.

- level:

  Confidence level for intervals (default is 0.95).

## Value

A data frame with columns `term`, `estimate`, `se`, `statistic`,
`p.value`, `conf.low`, `conf.high`, and `ci_method`. Columns `status`
and `reason` flag unsupported estimates (e.g., singular fits).
`usdt_tests()` includes an extra `hypothesis` column (`H1`, `H2`, `H3`).

For a Bayesian `hsdt` fit the columns are `term`, `estimate` (posterior
mean), `est.error` (posterior standard deviation), `conf.low` and
`conf.high` (central credible interval), `ci_method`, `p.value`
(two-sided posterior p-value, below `1 - level` exactly when the
interval excludes zero), `rhat` and `ess_bulk`.

## Details

These tests run automatically inside
[`hsdt()`](https://ricardoreysaez.github.io/uSDT/dev/reference/hsdt.md)
and appear in its summary. Calling them directly is especially useful
when fitting custom models with
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) or
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html),
allowing you to test these hypotheses while controlling for additional
covariates (e.g., set size, experimental groups).

## The three hypotheses

- **H1 (Mean difference):** Tests whether average sensitivity differs
  between the awareness and influence tasks.

- **H2 (Correlation):** Tests the correlation between task sensitivities
  across participants using a Fisher-\\z\\ transformed interval.

- **H3 (Latent regression):** Regresses influence sensitivity onto
  awareness sensitivity. The intercept represents expected influence
  performance when awareness sensitivity is zero
  (\\d'\_{\mathrm{Awareness}} = 0\\), testing for unconscious
  processing.

With positive sensitivity variances, H2 and the H3 slope share the same
zero-covariance null. The analytical frequentist output uses the same
covariance test for both. The correlation interval uses Fisher-\\z\\.

## Custom models with covariates

To adjust tests for additional factors, specify the model directly using
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) or
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html). As
long as the two sensitivity terms are included as fixed effects and
correlated across subjects via random slopes, `usdt_tests()` will
compute the latent tests conditional on those covariates. A brms model
is summarised from its posterior, like a Bayesian
[`hsdt()`](https://ricardoreysaez.github.io/uSDT/dev/reference/hsdt.md)
fit, and keeps the term names that lme4 would give it. Aggregating the
trials with
[`usdt_aggregate()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_aggregate.md)
gives either model the same estimates and makes it much faster to fit.

## See also

[`hsdt()`](https://ricardoreysaez.github.io/uSDT/dev/reference/hsdt.md),
[`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_boot.md),
[`usdt_aggregate()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_aggregate.md)

## Examples

``` r
# 1. Standard model via hsdt()
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

# All three tests at once
usdt_tests(m)
#>   hypothesis                          term   estimate         se  statistic
#> 1         H1 d'(influence) - d'(awareness) -0.1071347 0.03538288 -3.0278700
#> 2         H2                   correlation  0.4909823 0.51880441  1.0089996
#> 3         H3                     intercept  0.0445278 0.11585861  0.3843288
#> 4         H3                         slope  0.3557591 0.48650940  1.0089996
#>      p.value   conf.low   conf.high ci_method status reason
#> 1 0.00246284 -0.1764839 -0.03778559      Wald     ok   <NA>
#> 2 0.31297483 -0.6654131  0.95424003  Fisher-z     ok   <NA>
#> 3 0.70073479 -0.1825509  0.27160650     delta     ok   <NA>
#> 4 0.31297483 -0.5977818  1.30930002      Wald     ok   <NA>

# Or one test at a time
sensitivity_diff(m)   # H1
#>                            term   estimate         se statistic    p.value
#> 1 d'(influence) - d'(awareness) -0.1071347 0.03538288  -3.02787 0.00246284
#>     conf.low   conf.high ci_method status reason
#> 1 -0.1764839 -0.03778559      Wald     ok   <NA>
latent_cor(m)         # H2
#>          term  estimate        se statistic   p.value   conf.low conf.high
#> 1 correlation 0.4909823 0.5188044     1.009 0.3129748 -0.6654131   0.95424
#>   ci_method status reason
#> 1  Fisher-z     ok   <NA>
latent_regression(m)  # H3
#>        term  estimate        se statistic   p.value   conf.low conf.high
#> 1 intercept 0.0445278 0.1158586 0.3843288 0.7007348 -0.1825509 0.2716065
#> 2     slope 0.3557591 0.4865094 1.0089996 0.3129748 -0.5977818 1.3093000
#>   ci_method status reason
#> 1     delta     ok   <NA>
#> 2      Wald     ok   <NA>

# \donttest{
# 2. Custom model with covariates via glmer()
# Controlling for display set size across both tasks
trials <- rbind(
  data.frame(vadillo_awareness[c("subj", "condition", "set.size")],
             task = "A", resp = vadillo_awareness$judged.old),
  data.frame(vadillo_cuing[c("subj", "condition", "set.size")],
             task = "I", resp = meyen_split(vadillo_cuing$rt,
                                            by = vadillo_cuing$subj))
)

# Deviation contrasts (-0.5 vs 0.5). `awareness` flags the awareness task
trials <- within(trials, {
  cond   <- ifelse(condition == "old", 0.5, -0.5)
  size   <- ifelse(set.size == "set size 16", 0.5, -0.5)
  awareness <- as.integer(task == "A")
})

# Counting responses per subject, task and condition keeps the likelihood
# and makes the fit much faster than one row per trial
counts <- usdt_aggregate(trials, response = "resp",
                         by = c("subj", "task", "cond", "size", "awareness"))

# Standard glmer formula: influence criterion is omitted (fixed at 0
# by the median split). Random effects estimate the awareness criterion
# and correlated task sensitivities across subjects.
fit <- lme4::glmer(
  cbind(y, n - y) ~ 0 + awareness + task:size + task:cond +
    (0 + awareness | subj) + (0 + task:cond | subj),
  data = counts, family = binomial("probit"),
  control = lme4::glmerControl(optimizer = "bobyqa")
)

# Check the names lme4 assigned to the sensitivity terms
names(lme4::fixef(fit))
#> [1] "awareness"  "taskA:size" "taskI:size" "taskA:cond" "taskI:cond"

# Evaluate hypotheses conditional on set size
usdt_tests(fit, awareness = "taskA:cond", influence = "taskI:cond")
#>   hypothesis                          term    estimate         se  statistic
#> 1         H1 d'(influence) - d'(awareness) -0.11544701 0.03729773 -3.0952826
#> 2         H2                   correlation  0.38184383 0.39474844  0.9804973
#> 3         H3                     intercept  0.07364629 0.07325229  1.0053788
#> 4         H3                         slope  0.24001692 0.28643292  0.9804973
#>       p.value    conf.low   conf.high ci_method status reason
#> 1 0.001966254 -0.18854922 -0.04234481      Wald     ok   <NA>
#> 2 0.326840696 -0.46489648  0.86376154  Fisher-z     ok   <NA>
#> 3 0.314714491 -0.06992555  0.21721814     delta     ok   <NA>
#> 4 0.326840696 -0.32138129  0.80141513      Wald     ok   <NA>
# }

# \donttest{
# 3. The same model fitted with brms. Compiling it needs a C++ toolchain.
fit_b <- brms::brm(
  y | trials(n) ~ 0 + awareness + task:size + task:cond +
    (0 + awareness | subj) + (0 + task:cond | subj),
  data = counts, family = binomial("probit"),
  prior = brms::prior(normal(0, 1), class = b), seed = 1, refresh = 0
)
#> Compiling Stan program...
#> Start sampling
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess

# The same term names. The tests now summarise the posterior
usdt_tests(fit_b, awareness = "taskA:cond", influence = "taskI:cond")
#>   hypothesis                          term    estimate  est.error   conf.low
#> 1         H1 d'(influence) - d'(awareness) -0.11403474 0.03856358 -0.1887904
#> 2         H2                   correlation  0.25928434 0.40702225 -0.7085880
#> 3         H3                     intercept  0.13843497 1.79305979 -0.4160614
#> 4         H3                         slope -0.01880521 7.04631649 -2.5555180
#>     conf.high ci_method p.value     rhat  ess_bulk
#> 1 -0.03896435  quantile  0.0025 1.001837 4293.1610
#> 2  0.92546160  quantile  0.4725 1.009576  337.8252
#> 3  0.75631276  quantile  0.3695 1.008695  428.8847
#> 4  2.22775595  quantile  0.4725 1.008806  426.4601

# With a proper prior on the sensitivities, also Bayes factors
set.seed(1)
usdt_bf(fit_b, awareness = "taskA:cond", influence = "taskI:cond")

#>  uSDT hypothesis ─ Savage-Dickey density ratio test
#> ════════════════════════════════════════════════════════════════════════════════
#>   Hypothesis (H1)             : Δd' != 0
#>   Hypothesis (H0)             : Δd' = 0
#>   Posterior mean              : -0.114
#>   Posterior SD                : 0.039
#>   95% CrI                     : [-0.189, -0.039]
#>   log BF10                    : 0.706
#>   log BF01                    : -0.706
#>   Evidence                    : anecdotal for H1
#> ────────────────────────────────────────────────────────────────────────────────
#> 
#>  uSDT hypothesis ─ Savage-Dickey density ratio test
#> ════════════════════════════════════════════════════════════════════════════════
#>   Hypothesis (H1)             : rho != 0
#>   Hypothesis (H0)             : rho = 0
#>   Posterior mean              : 0.259
#>   Posterior SD                : 0.407
#>   95% CrI                     : [-0.709, 0.925]
#>   log BF10                    : -0.299
#>   log BF01                    : 0.299
#>   Evidence                    : anecdotal for H0
#> ────────────────────────────────────────────────────────────────────────────────
#> 
#>  uSDT hypothesis ─ Savage-Dickey density ratio test
#> ════════════════════════════════════════════════════════════════════════════════
#>   Hypothesis (H1)             : intercept != 0
#>   Hypothesis (H0)             : intercept = 0
#>   Posterior mean              : 0.138
#>   Posterior SD                : 1.793
#>   95% CrI                     : [-0.416, 0.756]
#>   log BF10                    : -1.762
#>   log BF01                    : 1.762
#>   Evidence                    : moderate for H0
#> ────────────────────────────────────────────────────────────────────────────────
# }
```
