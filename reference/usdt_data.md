# Prepare the data of a direct and an indirect task

Both functions build the same object. Use `usdt_data_tasks()` when each
task has its own data frame, and `usdt_data_long()` when a single data
frame holds both tasks together with a column that identifies them.

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

  The data frame of each task, for `usdt_data_tasks()`.
  `usdt_data_long()` ignores them and uses `task_levels` instead.

- subject_col:

  Name of the column that identifies the subject.

- condition_col:

  Name of the column that holds the signal and noise condition. It is
  not needed when the data arrive as an SDT table through `sdt_cols`.

- condition_levels:

  Which value of `condition_col` plays each role, as
  `c(signal = "old", noise = "new")`. The function guesses them and
  reports its choice when they are missing.

- response_col:

  Name of the column that holds the response. It can be a binary
  response, or a continuous measure such as response times when the task
  is dichotomized.

- response_levels:

  Which value of `response_col` counts as a signal response, as
  `c(signal = 1, noise = 0)`. A task that is dichotomized takes the side
  of the median instead, as `c(signal = "faster", noise = "slower")`.
  The function guesses them and reports its choice when they are
  missing.

- successes_col, trials_col:

  Names of the columns that hold data already summed up, the number or
  proportion of signal responses and the number of trials. Give these
  instead of `response_col`.

- successes_type:

  Format of `successes_col`. Use `"counts"` for counts and
  `"proportions"` for proportions. `"auto"` recognises clear cases and
  asks for an explicit choice when every value is zero or one.

- sdt_cols:

  Names of the columns of an SDT table, as
  `c(hit = "H", miss = "M", fa = "FA", cr = "CR")`. Give these instead
  of `condition_col` and `response_col`.

- dichotomize:

  Which tasks need the median split of
  [`meyen_split()`](https://ricardoreysaez.github.io/uSDT/reference/meyen_split.md).
  Name them with `"none"`, `"direct"`, `"indirect"` or `"both"`, or give
  one logical value per task, as
  `list(direct = FALSE, indirect = TRUE)`.

- ties:

  What to do with trials that fall exactly on the median. See
  [`meyen_split()`](https://ricardoreysaez.github.io/uSDT/reference/meyen_split.md).

- coding:

  How the condition enters the model. See Details.

- labels:

  Display names for the two tasks. They only affect printed output.

- data:

  A single data frame holding both tasks, for `usdt_data_long()`.

- task_col:

  Name of the column that identifies the task, for `usdt_data_long()`.

- task_levels:

  Which value of `task_col` belongs to each task, as
  `c(direct = "D", indirect = "I")`.

- x:

  A `usdt_data` object.

- ...:

  Ignored.

## Value

An object of class `usdt_data`. Its `agg` element is the data frame of
counts that the model uses, with one row per subject, task and
condition. Its `meta` element records how every column was read, which
tasks were split at the median, and the descriptive summaries shown when
the object is printed. Pass the object to
[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md).

## Details

The functions count the responses of every subject in each condition,
check that the two tasks describe the same subjects, and record how each
column was read. Printing the result shows all of that, so the coding
can be checked before the model runs.

## One value or one per task

The two tasks rarely come from the same experimental design, so every
argument that names a column, a level or a format accepts two forms.
Give one value and both tasks use it. Give one value per task and each
task is read on its own.

    subject_col      = "subj"                       # both tasks
    condition_col    = list(direct   = "condition",
                            indirect = "cue")       # one per task
    condition_levels = list(
      direct   = c(signal = "old",  noise = "new"),
      indirect = c(signal = "cued", noise = "uncued"))
    dichotomize      = list(direct = FALSE, indirect = TRUE)

This applies to `subject_col`, `condition_col`, `condition_levels`,
`response_col`, `response_levels`, `successes_col`, `trials_col`,
`successes_type`, `sdt_cols`, `dichotomize` and `ties`. The two tasks
may therefore use different columns, different values inside those
columns, and even different formats, with one task given trial by trial
and the other as an SDT table. Use
[`list()`](https://rdrr.io/r/base/list.html) rather than
[`c()`](https://rdrr.io/r/base/c.html) when the value of a task is
itself a vector, as happens with the `*_levels` and `sdt_cols`
arguments.

Only `coding` works differently. It describes the model itself, so it
always applies to both tasks at once.

## Condition coding

The two codings answer different questions and give different
intercepts. Under `"deviation"` the condition takes the values -0.5 and
+0.5, and the intercept is `-c`, the criterion measured from the
midpoint between the signal and noise distributions. Under `"treatment"`
the condition takes the values 0 and 1, and the intercept is `z(FAR)`,
the criterion measured from the noise distribution. The two intercepts
are related by `intercept_treatment = intercept_deviation - d'/2`.

The choice matters for a task that was split at the median. The split
leaves each subject with half signal responses, so with balanced
conditions the deviation intercept is exactly zero and needs no
estimation. The treatment intercept equals `-d'/2` instead, which is not
zero and has to be estimated. The object records this in
`criterion_zero`, and
[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md) uses
it to decide.

## See also

[`meyen_split()`](https://ricardoreysaez.github.io/uSDT/reference/meyen_split.md),
[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md),
[`sdt_moments()`](https://ricardoreysaez.github.io/uSDT/reference/sdt_moments.md)

## Examples

``` r
set.seed(1)
df <- usdt_simulate(n_subj = 30, n_trials = 80)

# Both tasks share every column name and every level here.
d  <- usdt_data_long(df, task_col = "task",
                     task_levels   = c(direct = "D", indirect = "I"),
                     subject_col   = "subj",
                     condition_col = "cond",
                     condition_levels = c(signal = 1, noise = 0),
                     response_col  = "response",
                     response_levels  = c(signal = 1, noise = 0))
d
#> ── Data summary ──────────────────────────────────────────────────────────────── 
#> 
#>   Input:          1 long data frame (usdt_data_long)
#>   Subjects:       30 (30 in both tasks, 0 in one only)
#>   Trials:         4,800 -> 120 aggregated rows (4 per subject)
#>   Coding:         deviation (condition coded -0.5 / +0.5; intercept estimates -c)
#>   Parameters:     10 (4 fixed effects, 6 (co)variance components)
#> 
#> ── Variable mapping ────────────────────────────────────────────────────────────
#> 
#>   Variable       Task       Column          Signal        Noise
#>   subject        Direct     subj            -             -
#>                  Indirect   subj            -             -
#>   condition      Direct     cond            1             0
#>                  Indirect   cond            1             0
#>   response       Direct     response        1             0
#>                  Indirect   response        1             0
#> 
#> ── Descriptives: median [min, max] across subjects ─────────────────────────────
#> 
#>   Task       Trials/cell                  HR                 FAR   d' (method-of-moments)
#>   Direct              40    .70 [ .23,  .93]    .38 [ .07,  .78]    0.90 [-0.63,  2.19]
#>   Indirect            40    .60 [ .38,  .80]    .45 [ .17,  .72]    0.28 [-0.13,  1.16]
#> 
#>   No cells at floor or ceiling.

# When they do not, each task gets its own column and its own levels.
aware <- df[df$task == "D", ]
cuing <- df[df$task == "I", ]
names(aware)[names(aware) == "cond"] <- "seen"
names(cuing)[names(cuing) == "cond"] <- "cue"
aware$seen <- ifelse(aware$seen == 1, "old",  "new")
cuing$cue  <- ifelse(cuing$cue  == 1, "cued", "uncued")

usdt_data_tasks(
  direct = aware, indirect = cuing,
  subject_col      = "subj",
  condition_col    = list(direct = "seen", indirect = "cue"),
  condition_levels = list(direct   = c(signal = "old",  noise = "new"),
                          indirect = c(signal = "cued", noise = "uncued")),
  response_col     = "response",
  response_levels  = c(signal = 1, noise = 0))
#> ── Data summary ──────────────────────────────────────────────────────────────── 
#> 
#>   Input:          2 data frames (usdt_data_tasks)
#>   Subjects:       30 (30 in both tasks, 0 in one only)
#>   Trials:         4,800 -> 120 aggregated rows (4 per subject)
#>   Coding:         deviation (condition coded -0.5 / +0.5; intercept estimates -c)
#>   Parameters:     10 (4 fixed effects, 6 (co)variance components)
#> 
#> ── Variable mapping ────────────────────────────────────────────────────────────
#> 
#>   Variable       Task       Column          Signal        Noise
#>   subject        Direct     subj            -             -
#>                  Indirect   subj            -             -
#>   condition      Direct     seen            old           new
#>                  Indirect   cue             cued          uncued
#>   response       Direct     response        1             0
#>                  Indirect   response        1             0
#> 
#> ── Descriptives: median [min, max] across subjects ─────────────────────────────
#> 
#>   Task       Trials/cell                  HR                 FAR   d' (method-of-moments)
#>   Direct              40    .70 [ .23,  .93]    .38 [ .07,  .78]    0.90 [-0.63,  2.19]
#>   Indirect            40    .60 [ .38,  .80]    .45 [ .17,  .72]    0.28 [-0.13,  1.16]
#> 
#>   No cells at floor or ceiling.
```
