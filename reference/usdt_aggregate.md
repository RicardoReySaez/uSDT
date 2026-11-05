# Count binary responses by group

Aggregates trial-level binary responses into the number of signal
responses (`y`) and the number of trials (`n`) for every combination of
the grouping columns. A binomial model fitted to these counts has the
same likelihood, and therefore the same estimates and posterior, as a
Bernoulli model fitted to the trials, and it runs much faster. Use it to
prepare custom models fitted with brms or lme4 (see
[usdt_hypotheses](https://ricardoreysaez.github.io/uSDT/reference/usdt_hypotheses.md)).

## Usage

``` r
usdt_aggregate(data, response, by)
```

## Arguments

- data:

  A data frame with one row per trial.

- response:

  Name of the binary response column, coded `0`/`1` or `FALSE`/`TRUE`.

- by:

  Names of the grouping columns. Every predictor of the model, the
  subject included, must be among them, so that all the trials of a row
  share their predictors.

## Value

A data frame with one row per combination of the `by` columns that
occurs in `data`, those columns, and the counts `y` and `n`.

## See also

[usdt_hypotheses](https://ricardoreysaez.github.io/uSDT/reference/usdt_hypotheses.md)

## Examples

``` r
# Awareness task of Vadillo et al. (2025): "old" judgements by subject,
# condition and set size
counts <- usdt_aggregate(vadillo_awareness, response = "judged.old",
                         by = c("subj", "condition", "set.size"))
head(counts)
#>   subj condition    set.size  y  n
#> 1 2001       new set size 16  9 16
#> 2 2002       new set size 16 14 16
#> 3 2003       new set size 16 14 16
#> 4 2004       new set size 16  8 16
#> 5 2005       new set size 16 11 16
#> 6 2006       new set size 16 11 16

# The counts keep every trial
sum(counts$n) == nrow(vadillo_awareness)
#> [1] TRUE
```
