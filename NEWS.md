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

## Simulation and example data

* `usdt_simulate()` generates data from the model for simulation studies.
* `vadillo_awareness` and `vadillo_cuing` provide trial-level data from
  Experiment 2 of Vadillo, Malejka and Shanks (2025).
