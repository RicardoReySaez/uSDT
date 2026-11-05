# uSDT (development version)

## Breaking changes

* Tasks are now named **awareness** and **influence** throughout the API.
  Use these names in `usdt_data_tasks()`, task-specific lists, `task_levels`,
  `labels`, `dichotomize` and sensitivity-term arguments. The old names are
  no longer accepted.
* The awareness code is `"A"`, with model terms `d_A` and `c_A`. Planning
  columns use `trials_awareness`, `trials_influence`,
  `reliability_awareness` and `reliability_influence`.

## Bayesian estimation

* `hsdt(estimation = "bayesian")` fits with brms and Stan and returns the
  `brmsfit` in `$fit`. Both rstan and cmdstanr are supported. Compiled models
  are cached. Defaults are four chains with 1,000 warmup iterations and 5,000
  kept draws per chain. Sampling reports divergences, tree depth, R-hat and
  effective sample sizes.
* `usdt_priors()` specifies priors for task means, SDs and correlations.
  `plot()` checks their predictive and hypothesis implications.
* `hsdt(latent = "regression")` places the influence-location prior on the
  latent intercept instead of the mean. This changes the induced prior of
  the mean difference. The intercept prior also affects inference on the
  correlation.
* `unequal_variances = TRUE` estimates signal SDs for Bayesian fits.
* `usdt_tests()` and the individual hypothesis functions accept Bayesian
  `hsdt` and custom brms fits. They report posterior means, SDs, central
  credible intervals and two-sided posterior sign p-values.
* `usdt_bf()` tests point, directional and interval hypotheses about `diff`,
  `rho`, `intercept` and `slope`. Point nulls use Savage-Dickey density ratios.
  Directions and regions compare posterior with prior odds. Reports give
  both log Bayes factors. By default, prior/posterior plots accompany them.
  Use `plot = FALSE` to suppress drawing and `plot(bf)` to draw later.
* Bayes factor plots label the tested hypothesis and quantity. A point-null
  panel reports BF01. A directional or interval panel reports BF10.
* `plot()` and `usdt_reliability()` support Bayesian fits.
* The [Bayesian tutorial](https://ricardoreysaez.github.io/uSDT/articles/uSDT-bayesian.html) covers priors,
  posterior summaries, Bayes factors and the latent regression prior.

## Equivalence and one-sided tests

* `usdt_equivalence()` accepts frequentist `hsdt` and custom glmer fits.
  Hypotheses use the syntax of `usdt_bf()`, with explicit bounds for
  equivalence (TOST), non-inferiority, one-sided and minimal-effect tests.
  Reports include the interval used for the decision. Bootstrap replicates
  are used when available. Unsupported inference is marked not estimable.

## Study planning

* `usdt_population()` defines population parameters, optionally from a
  frequentist pilot fit. Multiple values create scenarios. `update()`
  changes them. A latent intercept can replace the influence mean.
* `usdt_power()` evaluates designs or grids using a sampling-distribution
  approximation. It supports the hypotheses of `usdt_equivalence()`.
  `combine = "all"` adds joint rejection power. Non-estimable outcomes count
  as failures to reject, and populations in H0 are flagged.
* `usdt_power(method = "simulation")` simulates and fits whole studies,
  repeating the median split for dichotomized tasks. It reports Monte Carlo
  error, Wilson intervals and fitting failures. `summary()` adds bias,
  standard errors, interval widths and coverage. Results are reproducible
  across core counts. Unequal-variance populations can be simulated and
  analysed with the equal-variance model.
* `usdt_plan()` selects subjects, trials or both to meet a power or precision
  target at minimum cost, or maximize power within a budget. Ranges identify
  free values. Single numbers fix them. Check a chosen plan with
  `usdt_power(plan, method = "simulation")`.
* `usdt_mdes()` finds detectable effects or equivalence regions. `vary`
  selects another population parameter to vary.
* `plot()` supports populations, power curves, plans and detectable effects.
  The [planning tutorial](https://ricardoreysaez.github.io/uSDT/articles/uSDT-power.html) illustrates the workflow.

## Data preparation and plots

* `usdt_aggregate()` counts binary trials by grouping columns for custom
  binomial models with covariates.
* Plots share panel styling and captions that wrap to the figure width.
  Plot methods accept `caption = FALSE`.
* Regression plots accept `nrow = 2` to stack the observed and latent panels.
* ROC legends show sensitivity and AUC inside the panel.

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
