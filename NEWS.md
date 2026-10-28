# uSDT (development version)

## Bayesian estimation

* `hsdt(estimation = "bayesian")` fits the hierarchical model with brms on
  the aggregated hit and false-alarm counts, and returns the `brmsfit` in
  `$fit`. The priors go on the means, standard deviations and correlation of
  the two sensitivities and on the criteria, one brms prior per parameter and
  task. The model compiles once per set of priors into a per-user cache,
  through rstan or cmdstanr. Equal variances are the default; unequal
  variances are optional. By default four chains keep 5000 draws each, and up
  to four run at once, leaving two cores free. Each chain runs in a background
  session, so sampling shows a single progress bar and none of Stan's
  messages, followed by a summary of divergent transitions, tree depth, R-hat
  and effective sample size. As in brms, the backend defaults to the
  `brms.backend` option. brms, posterior, logspline, callr and cli are now
  imported, so the Bayesian estimation needs no optional package apart from
  `BH` and `RcppEigen` for the rstan backend, or cmdstanr.
* `hsdt(parameterization = "regression")` moves the prior of the indirect
  mean sensitivity to the intercept of the latent regression of the indirect
  on the direct one, so that both the intercept tested by H3 and the
  correlation tested by H2 have the priors you give them. The difference
  tested by H1 then has an induced, heavy-tailed prior, and the intercept
  prior also informs the posterior of the correlation. It samples the same
  parameters and reports the same hypotheses as the default
  parameterization.
* `usdt_priors()` specifies the priors as Stan-style strings, from any Stan
  distribution that covers the values of its parameter (the whole real line
  for the means; the positive values, by truncation at zero if needed, for
  the standard deviations; `lkj(eta)` for the correlations).
  `plot.usdt_priors()` shows what they imply before any data are fitted.
* `usdt_tests()`, `sensitivity_diff()`, `latent_cor()` and
  `latent_regression()` summarise a Bayesian fit with posterior means and
  standard deviations, central credible intervals and a two-sided posterior
  p-value that falls below .05 exactly when the 95% interval excludes zero.
  They also accept custom models fitted with brms (`brmsfit`), named by the
  same terms as an lme4 model.
* `usdt_bf()` computes Bayes factors for point (Savage-Dickey), directional
  and interval hypotheses about the difference, correlation, slope and
  intercept. A quantity with a prior of its own (the correlation, and under
  the latent regression the intercept) gets it exactly; the others, induced
  by the other priors, get Rao-Blackwell averages over 100,000 prior draws,
  whatever the prior families. A comparison of the slope with zero is tested
  on the correlation, which has its sign and a prior of its own. Point
  hypotheses read the posterior density from a logspline fit to the draws,
  used only if it stays close to the draws; directional and
  interval ones read the posterior probability as the share of draws in H1.
  Each hypothesis
  prints as a report, and `usdt_bf()` draws each prior against its posterior,
  in one row for up to three hypotheses, unless `plot = FALSE`. It also takes
  custom brms models, whose priors it reads from `prior_summary()`: proper
  priors from the families of `usdt_priors()` on the sensitivities and an
  LKJ prior on their correlation give the same priors as a Bayesian `hsdt()`
  fit.
* `plot()` draws the four plots of a Bayesian fit with posterior means and
  credible intervals, and `usdt_reliability()` gives the posterior of each
  task's reliability.
* The new article "Bayesian estimation with uSDT" walks through the Bayesian
  analysis of the tutorial data.

## Equivalence and one-sided tests

* `usdt_equivalence()` tests equivalence (TOST), non-inferiority and other
  one-sided, minimal-effect and point hypotheses about the difference,
  correlation, slope and intercept of a frequentist `hsdt()` fit or a custom
  glmer model, written as in `usdt_bf()`. Every hypothesis must be stated:
  the bounds are the smallest effects of interest, so there is no default.
  Each test reports the interval that decides it, of level `1 - 2 * alpha`
  for the one-sided and equivalence tests. The bounds, estimate and interval
  of the correlation are on the correlation scale, and its tests run on the
  Fisher-z scale. On a boundary fit the correlation, slope and intercept are
  marked not estimable, and after `usdt_boot()` every test uses the bootstrap
  replicates. Each hypothesis prints as a report.

## Study planning

* `usdt_population()` describes the population a planned study samples from:
  the means, standard deviations and correlation of the two sensitivities and
  the distribution of the criteria, with the argument names of
  `usdt_priors()`. The intercept of the latent regression can replace the
  indirect mean, as in `hsdt(parameterization = "regression")`, and the
  values can come from a frequentist `hsdt()` fit of a pilot study. Several
  values in an argument give one scenario each, and `update()` changes them.
  It is the first piece of the power analysis.
* `usdt_power()` gives the power of the tests of H1-H3 for a number of
  subjects and of trials per task, or for a grid of them, from such a
  population. Hypotheses are written as in `usdt_equivalence()`, so it covers
  point nulls at any value, one-sided, equivalence (TOST) and minimal-effect
  tests, and `combine = "all"` adds the power of rejecting all of them in the
  same study. The power is approximated from the sampling distribution of the
  subjects' sample mean and covariance of d', whose measurement error is the
  one `usdt_reliability()` reports. It counts the studies that end on the
  boundary, where the correlation and the latent regression are not
  estimable, as failures to reject them, reports how often that happens, and
  flags a population in H0, where the rejection rate is the type I error
  rate.
* `usdt_power(method = "simulation")` checks a design by simulating whole
  studies and analysing each one with `hsdt()`: a task split at the median
  draws latent scores per trial and is split with `meyen_split()`, as
  response times would be. It reports the simulated power with its Monte
  Carlo error and Wilson interval, and `summary()` the bias, empirical and
  model standard errors, interval width and coverage of each quantity. Each
  study has its own random number stream, so `ncores` runs them in parallel
  with the same results, and the simulation also accepts populations with
  unequal variances, analysed with the equal-variance model.
* `usdt_plan()` chooses the design. A range marks what is free and a number
  what is fixed, so it finds the subjects for given trials, the trials for
  given subjects, or both together at the lowest `cost`; with a `budget` it
  finds the most powerful design that costs no more, and with `precision`
  the cheapest design whose expected interval is narrow enough. It searches
  every number of trials with the closed form and checks the cheapest
  designs with the draws of `usdt_power()`, and it marks designs at the limit
  of a range. `usdt_power(plan, method = "simulation")` checks the chosen
  designs.
* `usdt_mdes()` finds the smallest effect a design detects with a given
  power for each hypothesis, moving the tested quantity with the rest of the
  population fixed as its parameterization says, or the narrowest
  equivalence region it can establish; `vary` moves another population
  value instead, such as the SD of the direct d'.

## Data preparation

* `usdt_aggregate()` counts binary trials into binomial `y`/`n` rows by any
  grouping columns, so custom lme4 and brms models with covariates fit on
  counts instead of trials, with the same estimates and much faster.

# uSDT 0.1.0

First release.

## Data preparation and descriptive estimates

* `usdt_data_tasks()` and `usdt_data_long()` prepare paired direct and indirect
  measures and report the column mappings and preparation diagnostics. Every
  column, level and format argument takes either one value for both tasks or
  one per task as `list(direct = ..., indirect = ...)`, so the two tasks may
  differ in the columns they use, in the values those columns take, and in the
  format they arrive in.
* `meyen_split()` dichotomizes continuous responses using a within-subject
  median pooled across conditions, following Meyen et al. (2022).
* `sdt_moments()` computes subject-level SDT estimates and optional sampling
  variances without fitting a hierarchical model.

## Model fitting and inference

* `hsdt()` fits a binomial probit mixed model with correlated random
  sensitivities and returns the three hypothesis tests.
* `sensitivity_diff()`, `latent_cor()` and `latent_regression()` calculate the
  mean sensitivity difference, latent correlation and latent regression.
  `usdt_tests()` collects the analytical tests in one table.
* `usdt_boot()` adds parametric bootstrap inference to the fitted object's
  hypothesis table when enough usable replicates are available.

## Visualization and reliability

* `plot()` for a fitted `hsdt` object draws observed and latent regressions,
  shrinkage, subject intervals and model-implied ROC curves.
* `usdt_reliability()` estimates task reliability from the fitted between-subject
  sensitivity variance and trial-level measurement variance.

## Example data

* `vadillo_awareness` and `vadillo_cuing` provide trial-level data from
  Experiment 2 of Vadillo, Malejka and Shanks (2025).
