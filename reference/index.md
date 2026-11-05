# Package index

## Prepare the data

Build paired task data, dichotomize continuous responses and count
trials for custom binomial models.

- [`usdt_data_tasks()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_data.md)
  [`usdt_data_long()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_data.md)
  [`print(`*`<usdt_data>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/usdt_data.md)
  : Prepare data for hierarchical SDT models
- [`meyen_split()`](https://ricardoreysaez.github.io/uSDT/reference/meyen_split.md)
  : Dichotomize response times into binary choices
- [`usdt_aggregate()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_aggregate.md)
  : Count binary responses by group

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

## Specify the Bayesian priors

Choose the prior distributions of the Bayesian model and check what they
imply.

- [`usdt_priors()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_priors.md)
  [`print(`*`<usdt_priors>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/usdt_priors.md)
  : Priors for the Bayesian hierarchical SDT model
- [`plot(`*`<usdt_priors>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/plot.usdt_priors.md)
  : Prior predictive check for the Bayesian hierarchical SDT model

## Test hypotheses and bootstrap uncertainty

sensitivity_diff() compares mean sensitivities. Latent_cor() estimates
their latent correlation. Latent_regression() estimates the intercept
and slope. usdt_tests() collects all three hypotheses. usdt_boot() adds
parametric bootstrap inference to a fitted uSDT model.
usdt_equivalence() tests equivalence, non-inferiority and minimal-effect
hypotheses about the same quantities.

- [`usdt_tests()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_hypotheses.md)
  [`sensitivity_diff()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_hypotheses.md)
  [`latent_cor()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_hypotheses.md)
  [`latent_regression()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_hypotheses.md)
  : Test the three core hypotheses of a hierarchical SDT model
- [`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_boot.md)
  : Parametric bootstrap intervals for hierarchical SDT models
- [`usdt_equivalence()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_equivalence.md)
  [`print(`*`<usdt_equivalence>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/usdt_equivalence.md)
  : Equivalence, non-inferiority and minimal-effect tests

## Bayes factors

Test point, directional and interval hypotheses about the difference,
correlation and latent regression of a Bayesian fit, and plot their
prior and posterior.

- [`usdt_bf()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_bf.md)
  [`print(`*`<usdt_bf>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/usdt_bf.md)
  : Bayes factors for the hypotheses of a Bayesian hierarchical SDT
  model
- [`plot(`*`<usdt_bf>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/plot.usdt_bf.md)
  : Plot the prior and posterior of tested quantities

## Plan a study

Describe the population a planned study samples from, in either
parameterization of the model or from the estimates of a pilot study.
Compute the power of the three hypotheses for a number of subjects and
of trials per task. Choose the subjects and trials at the lowest cost
and find the smallest effect a design detects.

- [`usdt_population()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_population.md)
  [`print(`*`<usdt_population>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/usdt_population.md)
  [`update(`*`<usdt_population>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/usdt_population.md)
  : Population values for planning a uSDT study
- [`usdt_power()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_power.md)
  [`print(`*`<usdt_power>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/usdt_power.md)
  [`summary(`*`<usdt_power>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/usdt_power.md)
  [`print(`*`<summary.usdt_power>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/usdt_power.md)
  : Power of a planned uSDT study
- [`usdt_plan()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_plan.md)
  [`print(`*`<usdt_plan>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/usdt_plan.md)
  : Plan the subjects and trials of a uSDT study
- [`usdt_mdes()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_mdes.md)
  [`print(`*`<usdt_mdes>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/usdt_mdes.md)
  : Smallest effect a uSDT study detects
- [`plot(`*`<usdt_population>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/plot_power.md)
  [`plot(`*`<usdt_power>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/plot_power.md)
  [`plot(`*`<usdt_mdes>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/plot_power.md)
  [`plot(`*`<usdt_plan>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/plot_power.md)
  : Plot a population, its power, its smallest effects or a plan

## Visualize the fit and estimate reliability

Compare observed and fitted estimates and quantify measurement
reliability.

- [`plot(`*`<hsdt>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/plot.hsdt.md)
  : Diagnostic and analytical plots for hierarchical SDT models
- [`usdt_reliability()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_reliability.md)
  [`summary(`*`<usdt_reliability>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/usdt_reliability.md)
  [`print(`*`<usdt_reliability>`*`)`](https://ricardoreysaez.github.io/uSDT/reference/usdt_reliability.md)
  : Reliability of awareness and influence task measures

## Example datasets

Trial-level data from Experiment 2 of Vadillo, Malejka and Shanks
(2025).

- [`vadillo_awareness`](https://ricardoreysaez.github.io/uSDT/reference/vadillo_awareness.md)
  : Awareness data from a probabilistic cuing experiment
- [`vadillo_cuing`](https://ricardoreysaez.github.io/uSDT/reference/vadillo_cuing.md)
  : Cuing data from a probabilistic cuing experiment
