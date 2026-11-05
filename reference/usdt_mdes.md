# Smallest effect a uSDT study detects

Finds, for each hypothesis, the smallest effect that a design detects
with the given power: the population value of the tested quantity, moved
away from H0 with the other values of the population fixed, at which the
approximate power of
[`usdt_power()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_power.md)
reaches it. For an equivalence test (`in`), it finds the narrowest
equivalence region, around the same centre, that the design can
establish.

## Usage

``` r
usdt_mdes(
  population,
  subjects,
  trials,
  hypothesis = c("diff = 0", "rho = 0", "intercept = 0"),
  power = 0.8,
  alpha = 0.05,
  signal = 0.5,
  dichotomize = c("none", "influence", "awareness", "both"),
  vary = NULL,
  range = NULL
)

# S3 method for class 'usdt_mdes'
print(x, digits = 3L, ...)
```

## Arguments

- population:

  A `usdt_population` object, or a `usdt_plan` from
  [`usdt_plan()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_plan.md),
  whose designs, hypotheses and settings are then used.

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
  [`usdt_equivalence()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_equivalence.md):
  `"<quantity> <op> <value>"` with the quantity `diff`, `rho`, `slope`
  or `intercept` and the operator `=`, `<` or `>`, or
  `"<quantity> in [a, b]"` and `"<quantity> out [a, b]"`. The default
  tests the three point nulls that
  [`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md)
  reports.

- power:

  Power to reach (default 0.80).

- alpha:

  Significance level of every test, between 0 and 0.5 (default 0.05).

- signal:

  Proportion of signal trials in each task (default 0.5), one value for
  both tasks or one for each. Each task must have a whole number of
  signal and of noise trials.

- dichotomize:

  Which tasks are split at the median, as in
  [`usdt_data_tasks()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_data.md):
  `"none"`, `"influence"`, `"awareness"`, `"both"`, or one logical value
  for each task. The split fixes the criterion of a task at zero, so its
  criterion values in `population` are not used.

- vary:

  Optional name of the population value to move instead of the tested
  quantity: `"mu_A"`, `"mu_I"`, `"sigma_A"`, `"sigma_I"`, `"rho"`,
  `"intercept"`, `"criterion_A"`, `"criterion_I"`, `"sd_criterion_A"` or
  `"sd_criterion_I"`. Needs `range`.

- range:

  With `vary`, the interval of values to search, `c(from, to)`. The
  result is the value nearest to `from` that reaches the power.

- x:

  A `usdt_mdes` object.

- digits:

  Number of decimals printed.

- ...:

  Ignored.

## Value

An object of class `usdt_mdes`: a data frame with one row per scenario,
design and hypothesis, and the columns `scenario`, `subjects`,
`trials_awareness`, `trials_influence`, `hypothesis`, `varies` (the
value moved), `population` (its value in the population),
`power_population` (the power there) and `mdes` (the value at which the
power is reached, or `NA` if it is not reached in the range searched).
For an equivalence test, `varies` is `"margin"` and `mdes` the
half-width of the narrowest region. The attribute `curve` holds the
power over the values searched.

## Details

Moving a quantity keeps the other values of the population fixed as its
parameterization says (see
[`usdt_population()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_population.md)):
the difference moves the influence mean with the awareness mean fixed.
The intercept moves the influence mean with the means and SDs and the
correlation fixed. The correlation (and the slope, through it) keeps
both means under the correlation parameterization and the intercept
under the latent regression. The search runs from the value of H0 (or
the nearest bound) outwards, in the direction of H1, and the power along
the way comes from the same draws of the sample mean and covariance as
in
[`usdt_power()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_power.md),
so it can level off below the target when fits on the boundary are
frequent. The result is then `NA`.

## See also

[`usdt_power()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_power.md),
[`usdt_plan()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_plan.md)

## Examples

``` r
pop <- usdt_population(dprime    = c(awareness = 0.15, influence = 0.30),
                       sd_dprime = c(awareness = 0.35, influence = 0.30),
                       cor_dprime = 0.3, criterion = 0.1)

usdt_mdes(pop, subjects = 80, trials = c(awareness = 64, influence = 200),
          dichotomize = "influence")
#> ── Smallest detectable effect ─────────────────────── approximate · power 0.8 ── 
#> 
#>   Design       80 subjects; trials per subject: 64 Awareness, 200 Influence
#> 
#>   Hypothesis        Varies     Population  Power there  Smallest
#>   H1 Δd' = 0        diff            0.150        0.712     0.168
#>   H2 rho = 0        rho             0.300        0.379     0.499
#>   H3 intercept = 0  intercept       0.261        0.992     0.135
#> 
#>   Approximate power, as in usdt_power(). NA: the power is not reached 
#>   in the range searched.

# The narrowest equivalence margin for the difference
usdt_mdes(update(pop, dprime = c(awareness = 0.15, influence = 0.15)),
          subjects = 80, trials = c(awareness = 64, influence = 200),
          hypothesis = "diff in [-0.2, 0.2]", dichotomize = "influence")
#> ── Smallest detectable effect ─────────────────────── approximate · power 0.8 ── 
#> 
#>   Design       80 subjects; trials per subject: 64 Awareness, 200 Influence
#> 
#>   Hypothesis             Varies  Population  Power there  Smallest
#>   H1 Δd' in [-0.2, 0.2]  margin       0.000        0.913     0.175
#> 
#>   Approximate power, as in usdt_power(). NA: the power is not reached 
#>   in the range searched.
#>   For an equivalence test, the smallest value is the half-width of the
#>   narrowest region the design can establish.

# The SD of the awareness d' the intercept test needs
usdt_mdes(pop, subjects = 80, trials = c(awareness = 64, influence = 200),
          hypothesis = "intercept > 0", vary = "sigma_A",
          range = c(0.05, 0.6), dichotomize = "influence")
#> ── Smallest detectable effect ─────────────────────── approximate · power 0.8 ── 
#> 
#>   Design       80 subjects; trials per subject: 64 Awareness, 200 Influence
#> 
#>   Hypothesis        Varies   Population  Power there  Smallest
#>   H3 intercept > 0  sigma_A       0.350        0.995     0.217
#> 
#>   Approximate power, as in usdt_power(). NA: the power is not reached 
#>   in the range searched.
```
