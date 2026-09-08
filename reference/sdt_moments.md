# Signal detection measures for each subject

Computes the hit rate, the false-alarm rate, d' and the criterion of
every subject, without fitting a model. It can also add the sampling
variance of d'. These descriptive values are useful to inspect the data
before fitting, and to compare with the model estimates afterwards.

## Usage

``` r
sdt_moments(
  data,
  subject_col = NULL,
  condition_col = NULL,
  condition_levels = NULL,
  response_col = NULL,
  response_levels = NULL,
  coding = c("deviation", "treatment"),
  correction = c("hautus", "none"),
  variances = FALSE
)
```

## Arguments

- data:

  A `usdt_data` object, or a plain data frame with one row per trial.
  With a `usdt_data` object the function processes both tasks and adds a
  `task` column.

- subject_col, condition_col, response_col:

  Names of the columns that hold the subject, the condition and the
  response. They are needed only for a plain data frame.

- condition_levels, response_levels:

  Which value plays each role, as `c(signal = ..., noise = ...)`. The
  function guesses them and reports its choice when they are missing.

- coding:

  Which criterion to report. `"deviation"` gives the classical criterion
  `c = -(z(HR) + z(FAR)) / 2`, measured from the midpoint between the
  two distributions. `"treatment"` gives `lambda = -z(FAR)`, measured
  from the noise distribution. The two are related by
  `lambda = c + d'/2`. A `usdt_data` object supplies its own coding.

- correction:

  What to do with rates of exactly 0 or 1, which make `d'` infinite.
  `"hautus"` adds 0.5 to the four counts of the affected subject.
  `"none"` leaves the infinite values in place.

- variances:

  If `TRUE`, adds the sampling variance of `d'` from Gourevitch and
  Galanter (1967) and from Miller (1996), together with the standard
  error that follows from the second one.

## Value

A data frame with one row per subject, or one row per subject and task
when `data` is a `usdt_data` object. It holds the four response counts
(`hit`, `miss`, `fa`, `cr`), the two rates (`hr`, `far`) and their
probit values (`zhr`, `zfar`), then `dprime`, `criterion`, and
`corrected` to mark the subjects that received the edge correction. With
`variances = TRUE` it also holds `var_gg`, `var_miller`, `e_miller` and
`se_dprime`.

## Details

The edge correction applies only to the subjects that need it. Applying
it to the whole sample would change the estimates of every other subject
as well, and those estimates are already usable.

The Miller variance treats the observed hit and false-alarm rates as
binomial probabilities. Samples that reach a rate of 0 or 1 receive the
values `0.5 / n` and `(n - 0.5) / n`. The number of trials stays the
original one throughout.

`var_gg` comes from the expected information of the two probit cells,
with the criterion treated as a nuisance parameter. It equals the
standard error that a probit regression would give for `d'` if it were
fitted to that subject alone.
[`usdt_reliability()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_reliability.md)
uses the same formula, but evaluates it at the rates the hierarchical
model predicts instead of the observed ones.

## References

Gourevitch, V., & Galanter, E. (1967). A significance test for one
parameter isosensitivity functions. *Psychometrika*.

Hautus, M. J. (1995). Corrections for extreme proportions and their
biasing effects on estimated values of d'. *Behavior Research Methods*.

Miller, J. (1996). The sampling distribution of d'. *Perception &
Psychophysics*.

Suero, M., Privado, J., & Botella, J. (2017). Methods to estimate the
variance of some indices of the signal detection theory: A simulation
study. *Psicologica*.

## See also

[`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md)

## Examples

``` r
set.seed(1)
df <- usdt_simulate(n_subj = 20, n_trials = 80)
head(sdt_moments(df[df$task == "D", ],
                 subject_col   = "subj",
                 condition_col = "cond",
                 condition_levels = c(signal = 1, noise = 0),
                 response_col  = "response",
                 response_levels  = c(signal = 1, noise = 0)))
#>   subj hit miss fa cr    hr   far        zhr       zfar      dprime   criterion
#> 1    1  27   13 20 20 0.675 0.500  0.4537622  0.0000000  0.45376219 -0.22688110
#> 2   10  24   16 15 25 0.600 0.375  0.2533471 -0.3186394  0.57198647  0.03264613
#> 3   11  30   10  8 32 0.750 0.200  0.6744898 -0.8416212  1.51611098  0.08356574
#> 4   12  23   17  9 31 0.575 0.225  0.1891184 -0.7554150  0.94453345  0.28314830
#> 5   13  25   15 14 26 0.625 0.350  0.3186394 -0.3853205  0.70395983  0.03334055
#> 6   14  16   24 17 23 0.400 0.425 -0.2533471 -0.1891184 -0.06422868  0.22123276
#>   corrected
#> 1     FALSE
#> 2     FALSE
#> 3     FALSE
#> 4     FALSE
#> 5     FALSE
#> 6     FALSE
```
