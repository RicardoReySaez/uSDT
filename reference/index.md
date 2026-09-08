# Package index

## Prepare the data

Build paired task data and dichotomize continuous responses.

- [`usdt_data_tasks()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_data.md)
  [`usdt_data_long()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_data.md)
  [`print(`*`<usdt_data>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/usdt_data.md)
  : Prepare data for hierarchical SDT models
- [`meyen_split()`](https://ricardoreysaez.github.io/uSDT/reference/meyen_split.md)
  : Dichotomize response times into binary choices

## Estimate descriptive SDT measures

Calculate subject-level SDT quantities without fitting a model.

- [`sdt_moments()`](https://ricardoreysaez.github.io/uSDT/reference/sdt_moments.md)
  : Signal detection measures for each subject

## Fit the hierarchical model

Estimate the paired sensitivities with a binomial probit mixed model.

- [`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md)
  [`summary(`*`<hsdt>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md)
  [`print(`*`<hsdt>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md)
  : Fit a hierarchical signal detection theory model

## Test hypotheses and bootstrap uncertainty

sensitivity_diff() compares mean sensitivities; latent_cor() estimates
their latent correlation; latent_regression() estimates the intercept
and slope. usdt_tests() collects all three hypotheses. usdt_boot() adds
parametric bootstrap inference to a fitted uSDT model.

- [`usdt_tests()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_hypotheses.md)
  [`sensitivity_diff()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_hypotheses.md)
  [`latent_cor()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_hypotheses.md)
  [`latent_regression()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_hypotheses.md)
  : Test the three hypotheses of a hierarchical SDT model
- [`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_boot.md)
  : Bootstrap intervals for a hierarchical SDT model

## Visualize the fit and estimate reliability

Compare observed and fitted estimates and quantify measurement
reliability.

- [`plot(`*`<hsdt>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/plot.hsdt.md)
  : Plot a fitted hierarchical SDT model
- [`usdt_reliability()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_reliability.md)
  [`summary(`*`<usdt_reliability>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/usdt_reliability.md)
  [`print(`*`<usdt_reliability>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/usdt_reliability.md)
  : Reliability of the direct and indirect measures

## Simulate data

Generate model-based data for simulation studies and power assessment.

- [`usdt_simulate()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_simulate.md)
  : Simulate data from a hierarchical SDT model

## Example datasets

Trial-level data from Experiment 2 of Vadillo, Malejka and Shanks
(2025).

- [`vadillo_awareness`](https://ricardoreysaez.github.io/uSDT/reference/vadillo_awareness.md)
  : Awareness data from a probabilistic cuing experiment
- [`vadillo_cuing`](https://ricardoreysaez.github.io/uSDT/reference/vadillo_cuing.md)
  : Cuing data from a probabilistic cuing experiment
