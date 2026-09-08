# Test the three hypotheses of a hierarchical SDT model

Computes the difference between the average sensitivities of the two
tasks, their correlation across subjects, and the regression of the
indirect sensitivity on the direct one. `usdt_tests()` returns the three
together, and the other three functions return one each.

## Usage

``` r
usdt_tests(fit, direct = "d_D", indirect = "d_I", level = 0.95)

sensitivity_diff(fit, direct = "d_D", indirect = "d_I", level = 0.95)

latent_cor(fit, direct = "d_D", indirect = "d_I", level = 0.95)

latent_regression(fit, direct = "d_D", indirect = "d_I", level = 0.95)
```

## Arguments

- fit:

  A fitted model, either an `hsdt` object or a `glmerMod` from
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html).

- direct, indirect:

  Names of the two sensitivity terms in the model.

- level:

  Confidence level.

## Value

A data frame with one row per quantity. The columns are `term`,
`estimate`, `se`, `statistic`, `p.value`, `conf.low`, `conf.high` and
`ci_method`. Two further columns, `status` and `reason`, mark the
results that the data cannot support and explain why. `usdt_tests()`
adds a `hypothesis` column with the values `H1`, `H2` and `H3`.

## Details

[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md)
already runs these tests, so most users read them in its summary.
Calling them directly is useful for a model fitted by hand, because they
accept any `glmerMod` in which the two sensitivities are fixed effects
and share a random-effects term.

H1 compares the two average sensitivities. A clear difference means that
the direct task measures more than the indirect one, or the reverse.

H2 gives the correlation between the two sensitivities across subjects.
It asks whether the subjects who are sensitive in one task are also the
sensitive ones in the other.

H3 regresses the indirect sensitivity on the direct one. Its intercept
is the sensitivity expected in the indirect task from a subject whose
direct sensitivity is zero, which is the test for unconscious
processing.

H1 and the two regression terms use Wald intervals. The correlation uses
a Fisher-z interval, so its limits stay between -1 and 1.

The slope of H3 is zero exactly when the covariance between the two
sensitivities is zero, and so is the correlation of H2. The two
therefore state the same null hypothesis, and both report the same test
on that covariance.

## See also

[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md),
[`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_boot.md)

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
usdt_tests(m$fit)
#>   hypothesis                      term  estimate         se statistic
#> 1         H1 d'(direct) - d'(indirect) 0.4673093 0.09034982 5.1722222
#> 2         H2               correlation 0.4679562 0.20065656 1.9258246
#> 3         H3                 intercept 0.1006160 0.13438542 0.7487121
#> 4         H3                     slope 0.3075323 0.14596629 1.9258246
#>        p.value     conf.low conf.high ci_method status reason
#> 1 2.313262e-07  0.290226956 0.6443917      Wald     ok   <NA>
#> 2 5.412628e-02  0.003902169 0.7661747  Fisher-z     ok   <NA>
#> 3 4.540307e-01 -0.162774588 0.3640066     delta     ok   <NA>
#> 4 5.412628e-02  0.021443616 0.5936210      Wald     ok   <NA>
# }
```
