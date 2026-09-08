# Plot a fitted hierarchical SDT model

Draws one of four plots from a fitted model. They show the relation
between the two tasks, how much the model corrects the observed values,
the sensitivity of each subject, and the implied ROC curves.

## Usage

``` r
# S3 method for class 'hsdt'
plot(
  x,
  type = c("regression", "shrinkage", "caterpillar", "roc"),
  subject_id = NULL,
  band = TRUE,
  population_reference = TRUE,
  ...
)
```

## Arguments

- x:

  A fitted `hsdt` object.

- type:

  Which plot to draw. `"regression"` compares the regression of the
  observed sensitivities with the one the model implies. `"shrinkage"`
  joins the observed `d'` of each subject to the estimate the model
  gives them. `"caterpillar"` shows the sensitivity of every subject
  with its interval, next to zero. `"roc"` draws the ROC curves that
  follow from the two sensitivities.

- subject_id:

  Subject to draw in an ROC plot. Without it the plot shows the curves
  of an average subject.

- band:

  Show the uncertainty band in a regression or ROC plot.

- population_reference:

  Add the average curves as thin dashed lines behind the curves of one
  subject in an ROC plot.

- ...:

  Reserved for future plot types.

## Value

A `ggplot` object. Its `data` component holds the values that the
selected plot shows.

## The regression plot

The plot has two panels that share their axes. The left panel shows the
observed sensitivities and the ordinary regression line through them.
Trial noise in the direct task pulls the slope of that line towards
zero, so it understates the relation between the tasks.

The right panel joins each observed value to the estimate the model
gives that subject, and draws the line
`gamma_I + beta1 * (x - gamma_D)`. Its slope divides the covariance of
the two sensitivities by the variance of the direct one, which removes
the effect of trial noise. The value of this line at `x = 0` is the
intercept that H3 reports, the sensitivity expected in the indirect task
from a subject with no direct sensitivity. The marker at `x = 0` reads
the same band, so it always agrees with the hypothesis table.

Each panel gives its intercept and slope with the corresponding
p-values. The band of the observed line is the usual confidence interval
of a least-squares fit. The band of the model line comes from the delta
method, or from the replicates once
[`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_boot.md)
has run.

## The shrinkage plot

The observed values come from
[`sdt_moments()`](https://ricardoreysaez.github.io/uSDT/reference/sdt_moments.md)
with the Hautus correction, and the model values add the departure of
each subject to the average of their task. A line joins the two values
of every subject. The grey contours describe the observed values and the
coloured contour describes the model values, which makes visible how
much the model pulls the extreme subjects towards the centre.

## The caterpillar plot

The intervals of the observed values use Miller standard errors and a
normal approximation. The model intervals cover the whole sensitivity of
the subject, so they carry the uncertainty of the task average, of the
variance components, and of the departure of that subject, along with
the relations among them. They are therefore only slightly wider than
the departures that
[`lme4::ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html)
returns alone. They become unavailable when the model cannot estimate
its full covariance safely. The percentage in the lower right corner of
each panel describes how many of the intervals shown include zero.

## The ROC plot

The curves follow the equal-variance identity
`HR = pnorm(qnorm(FAR) + dprime)`, and the point on each curve marks the
fitted criterion of that task. The average curves describe a subject
whose departures are zero. The curves of one subject use the estimates
the model gives that subject.

The band of an average curve transforms the interval of the sensitivity,
or the bootstrap replicates when they exist, using the interval type
chosen in
[`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_boot.md).
The band of a subject holds the average and the variance components
fixed. Every curve comes from the model. A binary response gives one
point per task, so these plots do not show an ROC curve measured across
several criteria.

## See also

[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md),
[`sdt_moments()`](https://ricardoreysaez.github.io/uSDT/reference/sdt_moments.md)

## Examples

``` r
# \donttest{
set.seed(1)
df <- usdt_simulate(n_subj = 40, n_trials = 100)
data <- usdt_data_long(
  df,
  task_col = "task",
  task_levels = c(direct = "D", indirect = "I"),
  subject_col = "subj",
  condition_col = "cond",
  condition_levels = c(signal = 1, noise = 0),
  response_col = "response",
  response_levels = c(signal = 1, noise = 0)
)
fit <- hsdt(data)
plot(fit, type = "regression")

plot(fit, type = "shrinkage")

plot(fit, type = "roc")

plot(fit, type = "roc", subject_id = "1")

# }
```
