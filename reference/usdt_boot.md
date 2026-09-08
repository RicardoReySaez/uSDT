# Bootstrap intervals for a hierarchical SDT model

Simulates many datasets from the fitted model, refits the model to each
one, and builds the intervals of the three hypotheses from the results.
This is useful when the ordinary intervals are unavailable or hard to
trust, which happens when the model reaches a boundary. The work is done
by [`lme4::bootMer()`](https://rdrr.io/pkg/lme4/man/bootMer.html).

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

  An `hsdt` object from
  [`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md).

- nsim:

  Number of usable replicates to reach. It must be at least 500.

- ncores:

  Number of cores to use. Values above one run the replicates in
  parallel through the `parallel` package, which comes with R. The
  temporary cluster behaves the same way on Windows, macOS and Linux,
  and it closes when the bootstrap ends.

- max_attempts:

  Largest number of replicates to fit. The default allows two attempts
  for every usable replicate requested.

- seed:

  Seed for the simulated datasets, so the result can be reproduced.

- progress:

  Show a progress bar. It appears by default in interactive sessions.

- level:

  Confidence level.

- type:

  Type of interval. `"perc"` takes the percentiles of the replicates,
  `"norm"` builds a normal interval around the bias-corrected estimate,
  and `"basic"` reflects the percentiles around the estimate. The last
  two work on the Fisher-z scale for the correlation, which keeps their
  limits inside its range. They need a finite centre on that scale, so a
  correlation that sits on the boundary reports them as missing.
  `"perc"` stays available in that case.

## Value

The `hsdt` object, with the interval columns of its `tests` table
replaced by the bootstrap results. The new `boot` element holds the
replicates of the three hypotheses in `t`, the sensitivity variances in
`variance`, the average task parameters in `population`, and the
estimates of every subject in `subjects`. It also holds the counts and
diagnostics of the run.

## Details

The function drops a replicate only when the model fails to fit or fails
to converge. It keeps singular and boundary replicates, because they are
the answer the model gives for a difficult dataset, and removing them
would make the intervals narrower than they should be. `boot$retained`
reports how many there were. The run continues until it reaches `nsim`
usable replicates or `max_attempts` fitted samples.

An incomplete run still returns the object, with every attempt and its
diagnostics in `boot`, and it gives a warning. Bootstrap summaries
replace the original intervals only from 500 usable replicates onwards.

The point estimates do not change. A bootstrap describes how much an
estimate would vary from sample to sample, and the estimate itself
remains the one the model produced. The bootstrap p-value compares the
fitted estimate in absolute value with the centred distribution of the
replicates. The count adds one to the numerator and the denominator, so
a finite simulation never returns a p-value of zero.

The `population` and `subjects` components keep four parameters from
every attempted replicate, the two criterion intercepts `c_D` and `c_I`
and the two sensitivities `d_D` and `d_I`. A criterion that the model
fixed is stored as zero, and the values of a subject combine the
refitted average with that subject's own departure from it. Their first
dimension follows `boot$ok`, so the same usable replicates can be
selected again.

## See also

[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md),
[`usdt_tests()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_hypotheses.md)

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
b <- usdt_boot(m, nsim = 500)
summary(b)
#> ── Model summary ─────────────────────────────────────────────────────────────── 
#> 
#>   Subjects:       40
#>   Observations:   160 aggregated rows (8,000 trials)
#>   Family:         binomial (probit)
#>   Coding:         deviation
#>   Criteria:       Direct estimated, Indirect estimated
#>   Estimation:     lme4::glmer (bobyqa)
#>   Convergence:    TRUE
#>   Bootstrap:      500 usable replicates (503 attempts; 8 at the boundary)
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
#>   cor(d')        both         0.4680   0.2211  [  0.017,  0.926]
#> 
#> ── Hypotheses ──────────────────────────────────────────────────────────────────
#> 
#> H1: Group-level sensitivity difference (Δd' = Direct d' - Indirect d')
#>   Parameter         Estimate  Boot SE  Boot 95% CI (percentile)  Boot p-value
#>   Δd' (D - I)        0.4673   0.0875  [  0.307,  0.637]                 .002
#> 
#> H2: Correlation between sensitivities across tasks
#>   Parameter         Estimate  Boot SE  Boot 95% CI (percentile)  Boot p-value
#>   rho                 0.4680   0.2211  [  0.017,  0.926]                 .050
#> 
#> H3: Latent regression of Indirect d' on Direct d'
#>   Parameter         Estimate  Boot SE  Boot 95% CI (percentile)  Boot p-value
#>   Intercept           0.1006   0.1454  [ -0.188,  0.370]                 .477
#>   Slope               0.3075   0.1573  [  0.010,  0.649]                 .068
# }
```
