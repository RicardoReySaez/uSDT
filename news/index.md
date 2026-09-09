# Changelog

## uSDT 0.1.0

First release.

### Data preparation and descriptive estimates

- [`usdt_data_tasks()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_data.md)
  and
  [`usdt_data_long()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_data.md)
  prepare paired direct and indirect measures and report the column
  mappings and preparation diagnostics. Every column, level and format
  argument takes either one value for both tasks or one per task as
  `list(direct = ..., indirect = ...)`, so the two tasks may differ in
  the columns they use, in the values those columns take, and in the
  format they arrive in.
- [`meyen_split()`](https://ricardoreysaez.github.io/uSDT/reference/meyen_split.md)
  dichotomizes continuous responses using a within-subject median pooled
  across conditions, following Meyen et al. (2022).
- [`sdt_moments()`](https://ricardoreysaez.github.io/uSDT/reference/sdt_moments.md)
  computes subject-level SDT estimates and optional sampling variances
  without fitting a hierarchical model.

### Model fitting and inference

- [`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md)
  fits a binomial probit mixed model with correlated random
  sensitivities and returns the three hypothesis tests.
- [`sensitivity_diff()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_hypotheses.md),
  [`latent_cor()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_hypotheses.md)
  and
  [`latent_regression()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_hypotheses.md)
  calculate the mean sensitivity difference, latent correlation and
  latent regression.
  [`usdt_tests()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_hypotheses.md)
  collects the analytical tests in one table.
- [`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_boot.md)
  adds parametric bootstrap inference to the fitted object’s hypothesis
  table when enough usable replicates are available.

### Visualization and reliability

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for a fitted
  `hsdt` object draws observed and latent regressions, shrinkage,
  subject intervals and model-implied ROC curves.
- [`usdt_reliability()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_reliability.md)
  estimates task reliability from the fitted between-subject sensitivity
  variance and trial-level measurement variance.

### Example data

- `vadillo_awareness` and `vadillo_cuing` provide trial-level data from
  Experiment 2 of Vadillo, Malejka and Shanks (2025).
