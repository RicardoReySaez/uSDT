# Test the three core hypotheses of a hierarchical SDT model

Evaluates the difference between average task sensitivities (H1), their
correlation across subjects (H2), and the regression of indirect
sensitivity on direct sensitivity (H3). `usdt_tests()` computes all
three together, while individual functions compute them separately.

## Usage

``` r
usdt_tests(fit, direct = "d_D", indirect = "d_I", level = 0.95)

sensitivity_diff(fit, direct = "d_D", indirect = "d_I", level = 0.95)

latent_cor(fit, direct = "d_D", indirect = "d_I", level = 0.95)

latent_regression(fit, direct = "d_D", indirect = "d_I", level = 0.95)
```

## Arguments

- fit:

  A fitted model: an `hsdt` object or a `glmerMod` from
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html).

- direct, indirect:

  Character strings naming the sensitivity terms in the model. Defaults
  match the internal naming of
  [`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md).
  For custom models, both terms must be fixed effects and share a common
  random-effects grouping by subject.

- level:

  Confidence level for intervals (default is 0.95).

## Value

A data frame with columns `term`, `estimate`, `se`, `statistic`,
`p.value`, `conf.low`, `conf.high`, and `ci_method`. Columns `status`
and `reason` flag unsupported estimates (e.g., singular fits).
`usdt_tests()` includes an extra `hypothesis` column (`H1`, `H2`, `H3`).

## Details

These tests run automatically inside
[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md) and
appear in its summary. Calling them directly is especially useful when
fitting custom models with
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html), allowing you
to test these hypotheses while controlling for additional covariates
(e.g., set size, experimental groups).

## The three hypotheses

- **H1 (Mean difference):** Tests whether average sensitivity differs
  between the direct and indirect tasks.

- **H2 (Correlation):** Tests the correlation between task sensitivities
  across participants using a Fisher-\\z\\ transformed interval.

- **H3 (Latent regression):** Regresses indirect sensitivity onto direct
  sensitivity. The intercept represents expected indirect performance
  when direct sensitivity is zero (\\d'\_{\mathrm{Direct}} = 0\\),
  testing for unconscious processing.

Because both the correlation (H2) and regression slope (H3) are zero if
and only if the covariance between sensitivities is zero, they evaluate
the same association and share identical test statistics.

## Custom models with covariates

To adjust tests for additional factors, specify the model directly using
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html). As long as
the two sensitivity terms are included as fixed effects and correlated
across subjects via random slopes, `usdt_tests()` will compute the
latent tests conditional on those covariates.

## See also

[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md),
[`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_boot.md)

## Examples

``` r
# \donttest{
# 1. Standard model via hsdt()
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

# All three tests at once
usdt_tests(m)
#>   hypothesis                      term    estimate         se  statistic
#> 1         H1 d'(indirect) - d'(direct) -0.10712457 0.03537881 -3.0279299
#> 2         H2               correlation  0.49116444 0.51904281  1.0090130
#> 3         H3                 intercept  0.04446212 0.11597753  0.3833684
#> 4         H3                     slope  0.35605319 0.48703535  1.0090130
#>       p.value   conf.low   conf.high ci_method status reason
#> 1 0.002462352 -0.1764658 -0.03778337      Wald     ok   <NA>
#> 2 0.312968398 -0.6657983  0.95434470  Fisher-z     ok   <NA>
#> 3 0.701446664 -0.1828497  0.27177390     delta     ok   <NA>
#> 4 0.312968398 -0.5985186  1.31062494      Wald     ok   <NA>

# Or one test at a time
sensitivity_diff(m)   # H1
#>                        term   estimate         se statistic     p.value
#> 1 d'(indirect) - d'(direct) -0.1071246 0.03537881  -3.02793 0.002462352
#>     conf.low   conf.high ci_method status reason
#> 1 -0.1764658 -0.03778337      Wald     ok   <NA>
latent_cor(m)         # H2
#>          term  estimate        se statistic   p.value   conf.low conf.high
#> 1 correlation 0.4911644 0.5190428  1.009013 0.3129684 -0.6657983 0.9543447
#>   ci_method status reason
#> 1  Fisher-z     ok   <NA>
latent_regression(m)  # H3
#>        term   estimate        se statistic   p.value   conf.low conf.high
#> 1 intercept 0.04446212 0.1159775 0.3833684 0.7014467 -0.1828497 0.2717739
#> 2     slope 0.35605319 0.4870354 1.0090130 0.3129684 -0.5985186 1.3106249
#>   ci_method status reason
#> 1     delta     ok   <NA>
#> 2      Wald     ok   <NA>


# 2. Custom model with covariates via glmer()
# Controlling for display set size across both tasks
trials <- rbind(
  data.frame(vadillo_awareness[c("subj", "condition", "set.size")],
             task = "D", resp = vadillo_awareness$judged.old),
  data.frame(vadillo_cuing[c("subj", "condition", "set.size")],
             task = "I", resp = meyen_split(vadillo_cuing$rt,
                                            by = vadillo_cuing$subj))
)

# Deviation contrasts (-0.5 vs 0.5); `direct` flags the direct task
trials <- within(trials, {
  cond   <- ifelse(condition == "old", 0.5, -0.5)
  size   <- ifelse(set.size == "set size 16", 0.5, -0.5)
  direct <- as.integer(task == "D")
})

# Standard glmer formula: indirect criterion is omitted (fixed at 0
# by the median split). Random effects estimate the direct criterion
# and correlated task sensitivities across subjects.
fit <- lme4::glmer(
  resp ~ 0 + direct + task:size + task:cond +
    (0 + direct | subj) + (0 + task:cond | subj),
  data = trials, family = binomial("probit"),
  control = lme4::glmerControl(optimizer = "bobyqa")
)

# Check the names lme4 assigned to the sensitivity terms
names(lme4::fixef(fit))
#> [1] "direct"     "taskD:size" "taskI:size" "taskD:cond" "taskI:cond"

# Evaluate hypotheses conditional on set size
usdt_tests(fit, direct = "taskD:cond", indirect = "taskI:cond")
#>   hypothesis                      term    estimate         se  statistic
#> 1         H1 d'(indirect) - d'(direct) -0.11543872 0.03708874 -3.1125005
#> 2         H2               correlation  0.38197036 0.39480466  0.9806745
#> 3         H3                 intercept  0.07360717 0.07328775  1.0043583
#> 4         H3                     slope  0.24018109 0.28657808  0.9806745
#>       p.value    conf.low   conf.high ci_method status reason
#> 1 0.001855097 -0.18813131 -0.04274613      Wald     ok   <NA>
#> 2 0.326753304 -0.46496183  0.86385791  Fisher-z     ok   <NA>
#> 3 0.315205936 -0.07003419  0.21724852     delta     ok   <NA>
#> 4 0.326753304 -0.32150163  0.80186381      Wald     ok   <NA>
# }
```
