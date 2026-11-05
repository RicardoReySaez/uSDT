# Power of a planned uSDT study

Computes the power of the tests of the three uSDT hypotheses, the
difference between mean sensitivities (H1), their correlation (H2) and
the latent regression (H3), for a study that samples a given number of
subjects, each with a given number of trials per task, from a population
described by
[`usdt_population()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_population.md).

## Usage

``` r
usdt_power(
  population,
  subjects,
  trials,
  hypothesis = c("diff = 0", "rho = 0", "intercept = 0"),
  alpha = 0.05,
  signal = 0.5,
  dichotomize = c("none", "influence", "awareness", "both"),
  combine = c("each", "all"),
  method = c("approximate", "simulation"),
  nsim = 1000,
  seed = NULL,
  ncores = 1L,
  progress = interactive()
)

# S3 method for class 'usdt_power'
print(x, digits = 3L, ...)

# S3 method for class 'usdt_power'
summary(object, ...)

# S3 method for class 'summary.usdt_power'
print(x, digits = 3L, ...)
```

## Arguments

- population:

  A `usdt_population` object from
  [`usdt_population()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_population.md).
  With several scenarios, the power of each. It can also be a
  `usdt_plan` from
  [`usdt_plan()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_plan.md),
  whose designs, hypotheses and settings are then used, so that
  `usdt_power(plan, method = "simulation")` checks a plan. `subjects`
  and `trials` are then left out.

- subjects:

  Number of subjects. Several values give one design each.

- trials:

  Number of trials per subject in each task, signal and noise trials
  together: one value for both tasks, or one for each as
  `c(awareness = , influence = )` or `list(awareness = , influence = )`.
  Several values give one design each. A value shared by both tasks
  keeps them equal, and values given for each task are crossed.

- hypothesis:

  Character vector of hypotheses, written as in
  [`usdt_equivalence()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_equivalence.md):
  `"<quantity> <op> <value>"` with the quantity `diff`, `rho`, `slope`
  or `intercept` and the operator `=`, `<` or `>`, or
  `"<quantity> in [a, b]"` and `"<quantity> out [a, b]"`. The default
  tests the three point nulls that
  [`hsdt()`](https://ricardoreysaez.github.io/uSDT/dev/reference/hsdt.md)
  reports.

- alpha:

  Significance level of every test, between 0 and 0.5 (default 0.05).

- signal:

  Proportion of signal trials in each task (default 0.5), one value for
  both tasks or one for each. Each task must have a whole number of
  signal and of noise trials.

- dichotomize:

  Which tasks are split at the median, as in
  [`usdt_data_tasks()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_data.md):
  `"none"`, `"influence"`, `"awareness"`, `"both"`, or one logical value
  for each task. The split fixes the criterion of a task at zero, so its
  criterion values in `population` are not used.

- combine:

  `"each"` reports the power of each hypothesis. `"all"` also reports
  the probability that every one of them rejects H0 in the same study.

- method:

  `"approximate"` (the default) computes power from the sampling
  distribution of the estimates, in seconds. `"simulation"` simulates
  whole studies and analyses each one with
  [`hsdt()`](https://ricardoreysaez.github.io/uSDT/dev/reference/hsdt.md),
  which takes minutes. See Details.

- nsim:

  Number of simulated studies for each scenario and design (default
  1000). Simulation only.

- seed:

  Random seed of the simulated studies. With the default `NULL`, the
  seed is drawn from the current random number stream. Simulation only.

- ncores:

  Number of CPU cores. Values above 1 fit the simulated studies in
  parallel on a temporary cluster, which works on Windows, macOS and
  Linux and stops when the run ends. Simulation only.

- progress:

  Logical. Show the expected duration and a progress bar (default `TRUE`
  in interactive sessions). Simulation only.

- x:

  A `usdt_power` object.

- digits:

  Number of decimals printed.

- ...:

  Ignored.

- object:

  A `usdt_power` object from a simulation.

## Value

An object of class `usdt_power`: a data frame with one row per scenario,
design and hypothesis, and the columns `scenario`, `subjects`,
`trials_awareness`, `trials_influence`, `hypothesis`, `test`, `H1`,
`H0`, `true` (the population value of the tested quantity), `se` (its
expected standard error), `power`, `estimable` (the share of studies in
which the analysis can estimate it), `in_null` (whether the population
lies in H0, so that `power` is the type I error rate), and
`reliability_awareness` and `reliability_influence`. A simulation adds
`mc_se`, the Monte Carlo standard error of `power`, and `power.low` and
`power.high`, its 95% Wilson interval. With `combine = "all"`, a row per
scenario and design whose `hypothesis` is `"all"` gives the joint power.
It prints one report for a single design and a table for several.

For a simulation, [`summary()`](https://rdrr.io/r/base/summary.html)
reports the precision of each quantity across the simulated studies that
estimate it: the mean estimate, its bias, the empirical standard error
(the SD of the estimates), the mean standard error the fits report, the
mean width of the intervals and their coverage of the population value.

## The approximation

Each subject's estimated \\d'\\ is treated as its true \\d'\\ plus a
measurement error. The variance of that error, \\v\\, is the one
[`usdt_reliability()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_reliability.md)
reports, averaged over the subjects of the population: it falls with the
number of trials and grows with extreme rates. The observed \\d'\\ of
the subjects then have covariance \\S = \Sigma + \mathrm{diag}(v_A,
v_I)\\, where \\\Sigma\\ is the covariance of the population, and the
reliability of each task is \\\sigma^2 / (\sigma^2 + v)\\. Trials act on
the tests only through \\v\\, and subjects through their number \\N\\.

With \\N\\ subjects, the sample mean of their \\d'\\ is normal with
covariance \\S / N\\ and their sample covariance has a Wishart
distribution. `usdt_power()` draws both 20,000 times, subtracts the
error variance to estimate \\\Sigma\\, and in each draw runs every test
as the fitted model would, with the standard errors it would compute at
those estimates. Power is the share of draws that reject H0, with a
Monte Carlo error of at most 0.004. The draws are the same for every
design, so power changes smoothly with the design, and they do not
change the random number stream of the session.

A draw whose estimated covariance of the sensitivities is not positive
definite (a variance of zero or a correlation of -1 or 1) is a fit on
the boundary. There, the fitted model reports the difference but marks
the correlation, slope and intercept as not estimable, so the draw
counts as a failure to reject them, and `estimable` reports how often
that happens. With an awareness task of low reliability \\R_A\\, this
probability is close to \\\Phi(-R_A \sqrt{N / 2})\\, and more awareness
trials lower it faster than more subjects.

## The tests

Each hypothesis is decided as the analysis would decide it: `"rho = 0"`
and `"slope = 0"` by the test on the covariance of
[`usdt_tests()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_hypotheses.md),
and every other statement by the tests of
[`usdt_equivalence()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_equivalence.md),
with the correlation on the Fisher-\\z\\ scale. `se` is the standard
error at the population values. When the population value lies in H0,
`power` is the type I error rate, and the report says so.

## Limits of the approximation

The approximation assumes the equal-variance model that
[`hsdt()`](https://ricardoreysaez.github.io/uSDT/dev/reference/hsdt.md)
fits, so every `sd_signal` of `population` must be 1, and it assumes the
estimated standard deviations of the criteria stay away from zero. It
follows the fitted model closely when the reliability of both tasks is
moderate or high. When the awareness task is unreliable, boundary fits
are more frequent in the fitted model than here, and the power of H2 and
H3 is somewhat optimistic.

## Simulation

With `method = "simulation"`, each study draws its subjects from the
population and their trials from their d' and criteria. A task split at
the median draws one latent score per trial and is split with
[`meyen_split()`](https://ricardoreysaez.github.io/uSDT/dev/reference/meyen_split.md),
as response times would be. Any other task draws binomial counts. The
study is then analysed by
[`hsdt()`](https://ricardoreysaez.github.io/uSDT/dev/reference/hsdt.md)
with the confidence level `1 - alpha`, and each hypothesis is decided as
in the approximation. A fit that fails, or that cannot estimate the
quantity, counts as a failure to reject. Each study has its own random
number stream, derived from `seed`, so a run gives the same studies
whatever `ncores`, and every design uses the same streams, which makes
the differences between designs less noisy than their Monte Carlo errors
suggest. The simulation also accepts a population with unequal variances
(`sd_signal` other than 1), whose studies are then analysed with the
equal-variance model. `se` and the reliabilities are then `NA`.

A simulation is slow, about a second per study with 100 subjects, so it
is best kept for checking a design chosen with the approximation, with
enough studies to bound the Monte Carlo error: 1000 studies give an
error of at most 0.016.

## See also

[`usdt_population()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_population.md),
[`usdt_tests()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_hypotheses.md),
[`usdt_equivalence()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_equivalence.md),
[`usdt_reliability()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_reliability.md)

## Examples

``` r
pop <- usdt_population(dprime    = c(awareness = 0.15, influence = 0.30),
                       sd_dprime = c(awareness = 0.35, influence = 0.30),
                       cor_dprime = 0.3, criterion = 0.1)

# 80 subjects, 64 awareness trials and 200 influence trials split at the
# median
usdt_power(pop, subjects = 80, trials = c(awareness = 64, influence = 200),
           dichotomize = "influence")
#> ── Power ──────────────────────────────────────────────────────── approximate ── 
#> 
#>   Design       80 subjects; trials per subject: 64 Awareness, 200 Influence (split at the median)
#>   Reliability  Awareness 0.543, Influence 0.738
#> 
#>   Hypothesis        Test              True      SE  Estimable   Power
#>   H1 Δd' = 0        Two-sided test   0.150   0.060      1.000   0.712
#>   H2 rho = 0        Two-sided test   0.300   0.168      1.000   0.379
#>   H3 intercept = 0  Two-sided test   0.261   0.045      1.000   0.992
#> 
#>   Approximate power: the share of 20,000 draws of the subjects' sample
#>   mean and covariance in which the test rejects H0, at alpha = 0.05.

# Several numbers of subjects, and one-sided and equivalence tests
usdt_power(pop, subjects = c(40, 80, 160),
           trials = c(awareness = 64, influence = 200),
           hypothesis = c("diff > 0", "intercept > 0", "rho in [-0.1, 0.1]"),
           dichotomize = "influence")
#> ── Power ──────────────────────────────────────────── approximate · 3 designs ── 
#> 
#>   H1: Δd' > 0, one-sided test
#>   H3: intercept > 0, one-sided test
#>   H2: rho in [-0.1, 0.1], equivalence test (TOST)
#> 
#>   Subjects  Awareness  Influence  Δd' > 0   intercept > 0   rho in [-0.1, 0.1] 
#>         40         64        200    0.562           0.935                0.000^
#>         80         64        200    0.809           0.995                0.000^
#>        160         64        200    0.972           1.000                0.000^
#> 
#>   ^ The population lies in H0: the rejection rate is the type I error rate.
#>   Approximate power: the share of 20,000 draws of the subjects' sample
#>   mean and covariance in which the test rejects H0, at alpha = 0.05.

# The same design under weaker correlations
usdt_power(update(pop, cor_dprime = c(0.1, 0.2)), subjects = 80,
           trials = c(awareness = 64,
                      influence = 200), dichotomize = "influence")
#> ── Power ────────────────────────────────────────── approximate · 2 scenarios ── 
#> 
#>   H1: Δd' = 0, two-sided test
#>   H2: rho = 0, two-sided test
#>   H3: intercept = 0, two-sided test
#> 
#>   Scenario  Subjects  Awareness  Influence  Δd' = 0  rho = 0  intercept = 0
#>          1        80         64        200    0.656    0.080          0.998
#>          2        80         64        200    0.684    0.188          0.996
#> 
#>   Approximate power: the share of 20,000 draws of the subjects' sample
#>   mean and covariance in which the test rejects H0, at alpha = 0.05.

# \donttest{
# Check the design by simulating and analysing 50 studies (under a
# minute). A real check needs more, e.g. nsim = 1000 with ncores = 4.
sims <- usdt_power(pop, subjects = 80,
                   trials = c(awareness = 64, influence = 200),
                   dichotomize = "influence", method = "simulation",
                   nsim = 50, seed = 1)
sims
#> ── Power ───────────────────────────────────────────────────────── simulation ── 
#> 
#>   Design       80 subjects; trials per subject: 64 Awareness, 200 Influence (split at the median)
#>   Reliability  Awareness 0.543, Influence 0.738
#> 
#>   Hypothesis        Test              True      SE  Estimable   Power   MC SE
#>   H1 Δd' = 0        Two-sided test   0.150   0.060      1.000   0.700   0.065
#>   H2 rho = 0        Two-sided test   0.300   0.168      1.000   0.400   0.069
#>   H3 intercept = 0  Two-sided test   0.261   0.045      1.000   0.980   0.020
#> 
#>   Simulated power: the share of 50 simulated studies, each analysed
#>   with hsdt(), that reject H0 at alpha = 0.05 (Monte Carlo SE at most 0.069).
#>   summary() reports the precision of the estimates.

# Bias, empirical and model standard errors, interval width and coverage
summary(sims)
#> ── Precision ──────────────────────────────── 50 simulated studies per design ── 
#> 
#>   80 subjects; 64 Awareness and 200 Influence trials
#> 
#>   Quantity    True   Mean    Bias  Emp. SE  Mean SE  Width  Coverage  Estimated
#>   Δd'        0.150  0.143  -0.007    0.064    0.059  0.230     0.900         50
#>   rho        0.300  0.322   0.022    0.186    0.172  0.654     0.940         50
#>   intercept  0.261  0.255  -0.006    0.052    0.049  0.193     0.980         50
#>   slope      0.257  0.288   0.031    0.194    0.177  0.695     0.960         50
#> 
#>   Width and coverage are those of the 95% intervals of usdt_tests(),
#>   over the studies that estimate each quantity.
# }
```
