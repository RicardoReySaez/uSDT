# Plan the subjects and trials of a uSDT study

Finds the design, the number of subjects and of trials per task, that
reaches a given power for the tests of the uSDT hypotheses at the lowest
cost, the design with the highest power within a budget, or the design
that estimates a quantity with a given precision. What is free is given
as a range and what is fixed as a single number, so the same function
finds the subjects for given trials, the trials for given subjects, or
both.

## Usage

``` r
usdt_plan(
  population,
  hypothesis = c("diff = 0", "rho = 0", "intercept = 0"),
  power = 0.8,
  subjects = c(10, 1000),
  trials,
  cost = NULL,
  budget = NULL,
  precision = NULL,
  combine = c("each", "all"),
  alpha = 0.05,
  signal = 0.5,
  dichotomize = c("none", "influence", "awareness", "both")
)

# S3 method for class 'usdt_plan'
print(x, digits = 3L, ...)
```

## Arguments

- population:

  A `usdt_population` object from
  [`usdt_population()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_population.md).
  With several scenarios, one plan for each.

- hypothesis:

  Character vector of hypotheses, written as in
  [`usdt_equivalence()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_equivalence.md):
  `"<quantity> <op> <value>"` with the quantity `diff`, `rho`, `slope`
  or `intercept` and the operator `=`, `<` or `>`, or
  `"<quantity> in [a, b]"` and `"<quantity> out [a, b]"`. The default
  tests the three point nulls that
  [`hsdt()`](https://ricardoreysaez.github.io/uSDT/dev/reference/hsdt.md)
  reports.

- power:

  Power to reach (default 0.80).

- subjects:

  Number of subjects: one number fixes it, and a range `c(min, max)`
  lets the plan choose it (default `c(10, 1000)`).

- trials:

  Number of trials per subject in each task, signal and noise trials
  together: one number fixes it, and a range `c(min, max)` lets the plan
  choose it. One value or range for both tasks, which keeps them equal,
  or one for each as `c(awareness = , influence = )` or
  `list(awareness = , influence = )`.

- cost:

  Cost of a design, which the plan minimises: a named vector with the
  cost of each `subject` (a session, recruitment) and of each
  `awareness` and `influence` trial, in any unit (missing names cost 0),
  so that a design costs
  `subjects * (subject + awareness * trials_awareness + influence * trials_influence)`,
  or a function of `subjects`, `awareness` and `influence` that returns
  the cost of each design. Needed when more than one of the subjects and
  the trials is free. With only one free, the plan minimises it.

- budget:

  Largest cost allowed. With a budget, the plan finds the design with
  the highest power that costs no more, and `power` is not used.

- precision:

  Named vector of the widths that the confidence intervals of
  usdt_tests() (level `1 - alpha`) should not exceed, by quantity:
  `diff`, `rho`, `intercept` or `slope`, e.g. `c(intercept = 0.2)`. With
  it, the plan finds the cheapest design whose expected interval is that
  narrow, and `hypothesis` and `power` are not used.

- combine:

  `"each"` plans a design for each hypothesis. `"all"` plans one design
  in which every hypothesis rejects H0 in the same study with the given
  power.

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

- x:

  A `usdt_plan` object.

- digits:

  Number of decimals printed.

- ...:

  Ignored.

## Value

An object of class `usdt_plan`: a list with `$design`, a data frame with
one row per scenario and hypothesis (or quantity, for a precision
target) and the columns `scenario`, `hypothesis`, `subjects`,
`trials_awareness`, `trials_influence`, `power` (approximate, as in
[`usdt_power()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_power.md)),
`estimable`, `width` (for a precision target), `cost` and `limit` (the
free argument whose limit the design reaches, if any). `$frontier`
records the smallest number of subjects that reaches the target with
each number of trials, from the closed form. `$checked` holds the
designs the search checked with the draws of the sample mean and
covariance, with the number of subjects and cost of each. `$settings`
stores the planning options. It prints one report per design.

## How the design is found

The search follows the two approximations of
[`usdt_power()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_power.md).
With each number of trials in the ranges, on a grid of at most 400
values per task, the closed form gives the smallest number of subjects
that reaches the power, and with it the cost of each design. The
cheapest designs are then checked with the draws of the sample mean and
covariance, which see the fits on the boundary, and the number of
subjects of each is raised or lowered until those draws reach the power.
The search stops when the next design costs more in the closed form than
the best design checked so far. A design whose subjects or trials reach
the limit of their range, or come within 2% of it, is marked: the range
may be too narrow.

The design chosen this way rests on the approximation. Simulate it with
`usdt_power(plan, method = "simulation")` before relying on it, as the
approximation is somewhat optimistic for H2 and H3 when the awareness
task is unreliable.

## Costs

Many designs reach the same power, trading subjects for trials, so the
choice among them depends on what each costs. With linear costs and the
difference between means (H1), the cheapest number of trials of each
task is close to \\\sqrt{(c_s / c_t) (a_t / \tau^2)}\\, where \\c_s\\ is
the cost of a subject, \\c_t\\ that of a trial of the task, \\a_t\\ the
error variance of a single trial (about \\2\pi\\ near \\d' = 0\\) and
\\\tau^2\\ the variance of the difference between the two \\d'\\ across
subjects, whatever the power.

## See also

[`usdt_power()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_power.md),
[`usdt_mdes()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_mdes.md),
[`usdt_population()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_population.md)

## Examples

``` r
pop <- usdt_population(dprime    = c(awareness = 0.15, influence = 0.30),
                       sd_dprime = c(awareness = 0.35, influence = 0.30),
                       cor_dprime = 0.3, criterion = 0.1)

# Subjects for given trials
usdt_plan(pop, hypothesis = "diff = 0",
          trials = c(awareness = 64, influence = 200),
          dichotomize = "influence")
#> ── Plan ───────────────────────────────────────────── approximate · power 0.8 ── 
#> 
#>   Subjects     10 to 1000
#>   Trials       Awareness 64 (fixed); Influence 200 (fixed)
#> 
#>   Target   Subjects  Awareness  Influence  Power  Estimable     Cost
#>   Δd' = 0       100         64        200  0.803      1.000      100
#> 
#>   Approximate power. Check the design with usdt_power(plan, method = "simulation").

# Subjects and trials together, with the session time of each in seconds
plan <- usdt_plan(pop, hypothesis = "diff = 0", power = 0.9,
                  subjects = c(20, 300),
                  trials = list(awareness = c(32, 160),
                                influence = c(100, 300)),
                  cost = c(subject = 600, awareness = 4, influence = 1.5),
                  dichotomize = "influence")
plan
#> ── Plan ───────────────────────────────────────────── approximate · power 0.9 ── 
#> 
#>   Subjects     20 to 300
#>   Trials       Awareness 32 to 160; Influence 100 to 300
#>   Cost         600 per subject + 4 per Awareness trial + 1.5 per Influence trial
#> 
#>   Target   Subjects  Awareness  Influence  Power  Estimable     Cost
#>   Δd' = 0       131         80        132  0.902      1.000  146,458
#> 
#>   Approximate power. Check the design with usdt_power(plan, method = "simulation").

# The most powerful design within 30 hours of sessions
usdt_plan(pop, hypothesis = "rho = 0", subjects = c(20, 400),
          trials = list(awareness = c(32, 256), influence = c(100, 800)),
          cost = c(subject = 600, awareness = 4, influence = 1.5),
          budget = 30 * 3600, dichotomize = "influence")
#> ── Plan ───────────────────────────────────────── approximate · budget 108000 ── 
#> 
#>   Subjects     20 to 400
#>   Trials       Awareness 32 to 256; Influence 100 to 800
#>   Cost         600 per subject + 4 per Awareness trial + 1.5 per Influence trial
#> 
#>   Target   Subjects  Awareness  Influence  Power  Estimable     Cost
#>   rho = 0        81        106        206  0.456      1.000  107,973
#> 
#>   Approximate power. Check the design with usdt_power(plan, method = "simulation").

# The intercept within an interval 0.2 wide
usdt_plan(pop, precision = c(intercept = 0.2),
          trials = c(awareness = 64, influence = 200),
          dichotomize = "influence")
#> ── Plan ───────────────────────────────────────────── approximate · precision ── 
#> 
#>   Subjects     10 to 1000
#>   Trials       Awareness 64 (fixed); Influence 200 (fixed)
#> 
#>   Target                  Subjects  Awareness  Influence  Width  Estimable     Cost
#>   intercept width <= 0.2        63         64        200  0.199      0.999       63
#> 
#>   Expected width of the 95% interval of usdt_tests() at the population values.

# \donttest{
# The intercept, whose search checks more designs (several seconds)
usdt_plan(pop, hypothesis = "intercept > 0", power = 0.9,
          subjects = c(20, 300),
          trials = list(awareness = c(32, 160), influence = c(100, 300)),
          cost = c(subject = 600, awareness = 4, influence = 1.5),
          dichotomize = "influence")
#> ── Plan ───────────────────────────────────────────── approximate · power 0.9 ── 
#> 
#>   Subjects     20 to 300
#>   Trials       Awareness 32 to 160; Influence 100 to 300
#>   Cost         600 per subject + 4 per Awareness trial + 1.5 per Influence trial
#> 
#>   Target         Subjects  Awareness  Influence  Power  Estimable     Cost
#>   intercept > 0        25        114        200  0.901      0.994   33,900
#> 
#>   ! intercept > 0: the design is at the limit of `subjects`; a wider range may give a cheaper one.
#>   Approximate power. Check the design with usdt_power(plan, method = "simulation").

# Check the chosen design by simulation
usdt_power(plan, method = "simulation", nsim = 50, seed = 1)
#> ── Power ───────────────────────────────────────────────────────── simulation ── 
#> 
#>   Design       131 subjects; trials per subject: 80 Awareness, 132 Influence (split at the median)
#>   Reliability  Awareness 0.597, Influence 0.650
#> 
#>   Hypothesis  Test              True      SE  Estimable   Power   MC SE
#>   H1 Δd' = 0  Two-sided test   0.150   0.046      1.000   0.920   0.038
#> 
#>   Simulated power: the share of 50 simulated studies, each analysed
#>   with hsdt(), that reject H0 at alpha = 0.05 (Monte Carlo SE at most 0.038).
#>   summary() reports the precision of the estimates.
# }
```
