# Equivalence, non-inferiority and minimal-effect tests

Tests interval, one-sided and point hypotheses about the quantities of
the three uSDT hypotheses in a frequentist
[`hsdt()`](https://ricardoreysaez.github.io/uSDT/dev/reference/hsdt.md)
fit or a custom model fitted with
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html): the
difference between mean sensitivities (`diff`, H1, reported as \\\Delta
d'\\, influence minus awareness), their correlation (`rho`, H2) and the
slope and intercept of the latent regression (`slope`, `intercept`, H3).
The hypotheses are written as in
[`usdt_bf()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_bf.md),
so a frequentist and a Bayesian analysis can state the same bounds.

## Usage

``` r
usdt_equivalence(
  fit,
  hypothesis,
  awareness = "d_A",
  influence = "d_I",
  alpha = 0.05
)

# S3 method for class 'usdt_equivalence'
print(x, digits = 3L, width = 80L, ...)
```

## Arguments

- fit:

  A frequentist `hsdt` object, with or without the replicates of
  [`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_boot.md),
  or a `glmerMod` from
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html).

- hypothesis:

  Character vector of hypotheses. It has no default: the bounds are the
  smallest effects of interest, which only the analysis can set. Each
  one is `"<quantity> <op> <value>"`, where the quantity is `diff`,
  `rho`, `slope` or `intercept` and the operator is `=`, `<` or `>`
  (`<=` and `>=` are read as `<` and `>`), or `"<quantity> in [a, b]"`
  and `"<quantity> out [a, b]"` for a region and its complement.

- awareness, influence:

  For a `glmerMod`, the fixed effects of the awareness and influence
  sensitivities, as in
  [`usdt_tests()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_hypotheses.md).
  Defaults match the internal naming of
  [`hsdt()`](https://ricardoreysaez.github.io/uSDT/dev/reference/hsdt.md).

- alpha:

  Significance level of every test, between 0 and 0.5 (default 0.05).

- x:

  A `usdt_equivalence` object.

- digits:

  Number of decimals printed.

- width:

  Width of the printed report, in characters.

- ...:

  Ignored.

## Value

An object of class `usdt_equivalence`: a data frame with one row per
hypothesis and the columns `hypothesis`, `test` (`"Two-sided test"`,
`"One-sided test"`, `"Equivalence test (TOST)"` or
`"Minimal-effect test"`), `H1`, `H0`, `estimate`, `se`, `conf.level`,
`conf.low` and `conf.high` (the interval that decides the test),
`ci_method`, `p.lower` and `p.upper` (the one-sided p-values at the two
bounds of a region), `p.value`, `decision`, `status` and `reason`. Rows
that cannot be estimated have `status = "not estimable"` and give the
`reason`. It prints one report per hypothesis.

## The tests

A point statement names the null hypothesis H0. An inequality or a
region names the alternative H1, as in
[`usdt_bf()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_bf.md).
Each test rejects H0 at level `alpha`:

- `q = c`: the two-sided test of H0: \\q = c\\.

- `q > c` and `q < c`: one-sided tests. With a margin on the
  unfavourable side of the reference value, `diff > -0.2` is a
  non-inferiority test and `diff < 0.2` a non-superiority test.

- `q in [a, b]`: an equivalence test by two one-sided tests (TOST,
  Schuirmann, 1987. Lakens et al., 2018). H0, \\q \le a\\ or \\q \ge
  b\\, is rejected when both one-sided tests reject their part of it, so
  `p.value` is the larger of `p.lower` and `p.upper`.

- `q out [a, b]`: a minimal-effect test. H0, \\a \le q \le b\\, is
  rejected when either one-sided test rejects it at `alpha / 2`, so
  `p.value` is twice the smaller of `p.lower` and `p.upper`. As the
  region shrinks to a point it becomes the two-sided test of that point.

Each test is decided by the position of a confidence interval, and the
row reports that interval: of level `1 - 2 * alpha` for the one-sided
and equivalence tests, and `1 - alpha` for the two-sided and
minimal-effect tests, which split `alpha` between two tails. With the
default `alpha = 0.05`, an equivalence test reports a 90% interval, and
H0 is rejected exactly when that interval lies inside the region.

## The correlation

The bounds of `rho` are given, and its estimate and interval reported,
on the correlation scale. The test runs on the Fisher-\\z\\ scale, where
the sampling distribution of the estimate is closer to normal: \$\$z =
\frac{\mathrm{atanh}(\hat\rho) - \mathrm{atanh}(c)}{SE / (1 -
\hat\rho^2)}.\$\$ Because
[`atanh()`](https://rdrr.io/r/base/Hyperbolic.html) is increasing,
\\\rho \> c\\ and \\\mathrm{atanh}(\rho) \> \mathrm{atanh}(c)\\ are the
same statement, so the p-value and the decision are those of the
hypothesis about `rho`, and the interval, transformed back, stays inside
(-1, 1).

Every p-value comes from the interval reported with it.
[`usdt_tests()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_hypotheses.md)
instead tests a null association on the covariance that `rho` and the
slope share, so its p-values for `rho = 0` and `slope = 0` can differ
slightly from these.

## The intercept and the slope

Both divide by the variance of the awareness sensitivity. When the
awareness sensitivity barely varies across subjects, their standard
errors are large and their Wald intervals inaccurate, so their
equivalence tests need many subjects to succeed, and the bootstrap of
[`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_boot.md)
gives more accurate intervals.

## Boundary fits and the bootstrap

The correlation, slope and intercept need the full covariance matrix of
the fit, which is unavailable when the sensitivity covariance is on the
boundary (a variance estimated as zero, or a correlation of -1 or 1).
Their rows are then marked `"not estimable"`, while the difference,
which only needs the fixed effects, is still tested. After
[`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_boot.md),
every test uses the bootstrap replicates and the interval type chosen
there. The p-value at a bound is the share of the bootstrap distribution
beyond it, counted as \\(k + 1) / (B + 1)\\ as in
[`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_boot.md),
so the decision agrees with the reported interval up to that correction.

## References

Lakens, D., Scheel, A. M., & Isager, P. M. (2018). Equivalence testing
for psychological research: A tutorial. *Advances in Methods and
Practices in Psychological Science*, 1(2), 259-269.
[doi:10.1177/2515245918770963](https://doi.org/10.1177/2515245918770963)

Schuirmann, D. J. (1987). A comparison of the two one-sided tests
procedure and the power approach for assessing the equivalence of
average bioavailability. *Journal of Pharmacokinetics and
Biopharmaceutics*, 15(6), 657-680.
[doi:10.1007/BF01068419](https://doi.org/10.1007/BF01068419)

## See also

[`usdt_tests()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_hypotheses.md),
[`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_boot.md),
[`usdt_bf()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_bf.md)

## Examples

``` r
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

# Mean sensitivities equivalent within 0.2 d' units, and an unconscious
# effect at zero awareness smaller than 0.2
usdt_equivalence(m, c("diff in [-0.2, 0.2]", "intercept < 0.2"))
#>  uSDT hypothesis ─ Equivalence test (TOST)
#> ════════════════════════════════════════════════════════════════════════════════
#>   Hypothesis (H1)             : Δd' in [-0.2, 0.2]
#>   Hypothesis (H0)             : Δd' out [-0.2, 0.2]
#>   Estimate                    : -0.107
#>   Standard error              : 0.035
#>   90% CI, Wald                : [-0.165, -0.049]
#>   p-value, lower bound        : .004
#>   p-value, upper bound        : <.001
#>   p-value                     : .004
#>   Decision (alpha = 0.05)     : H0 rejected, in favour of H1
#> ────────────────────────────────────────────────────────────────────────────────
#> 
#>  uSDT hypothesis ─ One-sided test
#> ════════════════════════════════════════════════════════════════════════════════
#>   Hypothesis (H1)             : intercept < 0.2
#>   Hypothesis (H0)             : intercept >= 0.2
#>   Estimate                    : 0.045
#>   Standard error              : 0.116
#>   90% CI, Wald                : [-0.146, 0.235]
#>   p-value                     : .090
#>   Decision (alpha = 0.05)     : H0 not rejected
#> ────────────────────────────────────────────────────────────────────────────────

# The bounds of the correlation are correlations
usdt_equivalence(m, "rho in [-0.3, 0.3]")
#>  uSDT hypothesis ─ Equivalence test (TOST)
#> ════════════════════════════════════════════════════════════════════════════════
#>   Hypothesis (H1)             : rho in [-0.3, 0.3]
#>   Hypothesis (H0)             : rho out [-0.3, 0.3]
#>   Estimate                    : 0.491
#>   Standard error              : 0.519
#>   90% CI, Fisher-z            : [-0.528, 0.930]
#>   p-value, lower bound        : .108
#>   p-value, upper bound        : .631
#>   p-value                     : .631
#>   Decision (alpha = 0.05)     : H0 not rejected
#> ────────────────────────────────────────────────────────────────────────────────
```
