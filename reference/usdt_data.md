# Prepare data for hierarchical SDT models

Formats direct and indirect task data into a standardized structure for
[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md). Use
`usdt_data_tasks()` when tasks are stored in separate data frames, or
`usdt_data_long()` when both tasks are kept in a single data frame with
a column that identifies each task.

## Usage

``` r
usdt_data_tasks(
  direct,
  indirect,
  subject_col,
  condition_col = NULL,
  condition_levels = NULL,
  response_col = NULL,
  response_levels = NULL,
  successes_col = NULL,
  trials_col = NULL,
  successes_type = c("auto", "counts", "proportions"),
  sdt_cols = NULL,
  dichotomize = c("none", "indirect", "direct", "both"),
  ties = c("noise", "random"),
  coding = c("deviation", "treatment"),
  labels = c(direct = "Direct", indirect = "Indirect")
)

usdt_data_long(
  data,
  task_col,
  task_levels,
  subject_col,
  condition_col = NULL,
  condition_levels = NULL,
  response_col = NULL,
  response_levels = NULL,
  successes_col = NULL,
  trials_col = NULL,
  successes_type = c("auto", "counts", "proportions"),
  sdt_cols = NULL,
  dichotomize = c("none", "indirect", "direct", "both"),
  ties = c("noise", "random"),
  coding = c("deviation", "treatment"),
  labels = NULL
)

# S3 method for class 'usdt_data'
print(x, ...)
```

## Arguments

- direct, indirect:

  Data frames for each task (used in `usdt_data_tasks()`).

- subject_col:

  Name of the column that identifies participants.

- condition_col:

  Name of the column for signal and noise conditions. Not needed when
  using `sdt_cols`.

- condition_levels:

  Named vector mapping condition labels, like
  `c(signal = "old", noise = "new")`. Guessed automatically if left
  empty.

- response_col:

  Name of the column with responses. Can be binary choices or continuous
  values (like response times) to split at the median.

- response_levels:

  Named vector mapping responses, like `c(signal = 1, noise = 0)` or
  `c(signal = "faster", noise = "slower")`. Guessed automatically if
  left empty.

- successes_col, trials_col:

  Names of columns with pre-calculated counts or proportions of signal
  responses and total trials. Use these instead of `response_col`.

- successes_type:

  Format of `successes_col`: `"counts"`, `"proportions"`, or `"auto"`.

- sdt_cols:

  Named vector for SDT table columns, like
  `c(hit = "H", miss = "M", fa = "FA", cr = "CR")`.

- dichotomize:

  Which tasks to split at the median using
  [`meyen_split()`](https://ricardoreysaez.github.io/uSDT/reference/meyen_split.md).
  Use `"none"`, `"direct"`, `"indirect"`, `"both"`, or a list like
  `list(direct = FALSE, indirect = TRUE)`.

- ties:

  How to handle trials that fall exactly on the median. See
  [`meyen_split()`](https://ricardoreysaez.github.io/uSDT/reference/meyen_split.md).

- coding:

  How condition is coded in the model: `"deviation"` (-0.5, 0.5) or
  `"treatment"` (0, 1). See Details.

- labels:

  Optional names for the tasks in printed output.

- data:

  A single data frame with both tasks (used in `usdt_data_long()`).

- task_col:

  Name of the column that identifies the task in `data`.

- task_levels:

  Named vector mapping task labels, like
  `c(direct = "D", indirect = "I")`.

- x:

  A `usdt_data` object.

- ...:

  Ignored.

## Value

An object of class `usdt_data`. The `$agg` table contains the counts
used by
[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md), and
`$meta` contains setup details and summaries.

## Details

The functions count responses for each subject and condition, check that
the same subjects appear in both tasks, and print a summary table so you
can verify the column settings before fitting the model.

## Settings per task

Arguments for columns and levels take either a single value (used for
both tasks) or a list with separate settings for each task:

    condition_col    = list(direct = "cond", indirect = "cue")
    condition_levels = list(direct   = c(signal = "old", noise = "new"),
                            indirect = c(signal = "congruent", noise = "incongruent"))
    dichotomize      = list(direct = FALSE, indirect = TRUE)

This works for all column and level arguments, so you can combine
trial-level data in one task with summary tables in the other.

## Condition coding

Under `"deviation"` coding (-0.5 vs. +0.5), the intercept is \\-c\\, the
criterion measured from the point between the two distributions. Under
`"treatment"` coding (0 vs. 1), the intercept is \\z(\mathrm{FAR})\\,
the criterion measured from the noise distribution.

When a task is split at the median, deviation coding sets the group
criterion to zero in balanced designs, so the model does not need to
estimate it.

## See also

[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md),
[`meyen_split()`](https://ricardoreysaez.github.io/uSDT/reference/meyen_split.md),
[`sdt_moments()`](https://ricardoreysaez.github.io/uSDT/reference/sdt_moments.md)

## Examples

``` r
# 1. Tasks in separate data frames
# Direct task: binary choices (old/new)
# Indirect task: response times (split at the median)
d_separate <- usdt_data_tasks(
  direct   = vadillo_awareness,
  indirect = vadillo_cuing,
  subject_col      = "subj",
  condition_col    = "condition",
  condition_levels = c(signal = "old", noise = "new"),
  response_col     = list(direct = "judged.old", indirect = "rt"),
  response_levels  = list(direct   = c(signal = 1, noise = 0),
                          indirect = c(signal = "faster", noise = "slower")),
  dichotomize      = list(direct = FALSE, indirect = TRUE)
)
d_separate
#> ── Data summary ──────────────────────────────────────────────────────────────── 
#> 
#>   Input:          2 data frames (usdt_data_tasks)
#>   Subjects:       104 (104 in both tasks, 0 in one only)
#>   Trials:         46,592 -> 416 aggregated rows (4 per subject)
#>   Coding:         deviation (condition coded -0.5 / +0.5; intercept estimates -c)
#>   Parameters:     7 (3 fixed effects, 4 (co)variance components)
#> 
#> ── Variable mapping ────────────────────────────────────────────────────────────
#> 
#>   Variable       Task       Column          Signal        Noise
#>   subject        Direct     subj            -             -
#>                  Indirect   subj            -             -
#>   condition      Direct     condition       old           new
#>                  Indirect   condition       old           new
#>   response       Direct     judged.old      1             0
#>                  Indirect   rt              faster        slower  [Meyen split]
#> 
#> ── Descriptives: median [min, max] across subjects ─────────────────────────────
#> 
#>   Task       Trials/cell                  HR                 FAR   d' (method-of-moments)
#>   Direct              32    .56 [ .28,  .88]    .47 [ .12,  .78]    0.24 [-0.58,  1.38]
#>   Indirect           192    .53 [ .44,  .60]    .47 [ .40,  .56]    0.13 [-0.31,  0.53]
#> 
#>   No cells at floor or ceiling.

# 2. Tasks combined in a single long data frame
long <- rbind(
  data.frame(task      = "D",
             subj      = vadillo_awareness$subj,
             condition = vadillo_awareness$condition,
             response  = vadillo_awareness$judged.old),
  data.frame(task      = "I",
             subj      = vadillo_cuing$subj,
             condition = vadillo_cuing$condition,
             response  = vadillo_cuing$rt)
)

d_long <- usdt_data_long(
  long,
  task_col         = "task",
  task_levels      = c(direct = "D", indirect = "I"),
  subject_col      = "subj",
  condition_col    = "condition",
  condition_levels = c(signal = "old", noise = "new"),
  response_col     = "response",
  response_levels  = list(direct   = c(signal = 1, noise = 0),
                          indirect = c(signal = "faster", noise = "slower")),
  dichotomize      = list(direct = FALSE, indirect = TRUE)
)
d_long
#> ── Data summary ──────────────────────────────────────────────────────────────── 
#> 
#>   Input:          1 long data frame (usdt_data_long)
#>   Subjects:       104 (104 in both tasks, 0 in one only)
#>   Trials:         46,592 -> 416 aggregated rows (4 per subject)
#>   Coding:         deviation (condition coded -0.5 / +0.5; intercept estimates -c)
#>   Parameters:     7 (3 fixed effects, 4 (co)variance components)
#> 
#> ── Variable mapping ────────────────────────────────────────────────────────────
#> 
#>   Variable       Task       Column          Signal        Noise
#>   subject        Direct     subj            -             -
#>                  Indirect   subj            -             -
#>   condition      Direct     condition       old           new
#>                  Indirect   condition       old           new
#>   response       Direct     response        1             0
#>                  Indirect   response        faster        slower  [Meyen split]
#> 
#> ── Descriptives: median [min, max] across subjects ─────────────────────────────
#> 
#>   Task       Trials/cell                  HR                 FAR   d' (method-of-moments)
#>   Direct              32    .56 [ .28,  .88]    .47 [ .12,  .78]    0.24 [-0.58,  1.38]
#>   Indirect           192    .53 [ .44,  .60]    .47 [ .40,  .56]    0.13 [-0.31,  0.53]
#> 
#>   No cells at floor or ceiling.
```
