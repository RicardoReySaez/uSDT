# Awareness data from a probabilistic cuing experiment

Trial-level direct-awareness data from Experiment 2 of Vadillo et al.
(2025). The object reproduces the source CSV file without filtering or
recoding and can be used as the `direct` input to
[`usdt_data_tasks()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_data.md).

## Usage

``` r
vadillo_awareness
```

## Format

A data frame with 6,656 rows, 104 participants, and 10 variables:

- `subj`:

  Participant identifier.

- `file`:

  Original participant file name.

- `trial`:

  Trial number.

- `pattId`:

  Pattern identifier.

- `judgm`:

  Awareness judgement on the original response scale.

- `judged.old`:

  Whether the pattern was judged old (1) or new (0).

- `condition`:

  Whether the pattern was old or new.

- `offset`:

  Spatial-offset condition.

- `color`:

  Colour condition.

- `set.size`:

  Set-size condition.

## Source

Vadillo, M. A., Malejka, S., & Shanks, D. R. (2025). Mapping the
reliability multiverse of contextual cuing. *Journal of Experimental
Psychology: Learning, Memory, and Cognition*.
[doi:10.1037/xlm0001410](https://doi.org/10.1037/xlm0001410) . Data
retrieved from <https://osf.io/jp3gx/>.

## Examples

``` r
data(vadillo_awareness)
str(vadillo_awareness)
#> 'data.frame':    6656 obs. of  10 variables:
#>  $ subj      : int  2001 2001 2001 2001 2001 2001 2001 2001 2001 2001 ...
#>  $ file      : chr  "y109b_subj1051.mat" "y109b_subj1051.mat" "y109b_subj1051.mat" "y109b_subj1051.mat" ...
#>  $ trial     : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ pattId    : int  1001 1002 3005 3006 3001 1006 3002 1007 3007 3008 ...
#>  $ judgm     : int  1 5 4 5 2 4 3 5 3 6 ...
#>  $ judged.old: int  0 1 1 1 0 1 0 1 0 1 ...
#>  $ condition : chr  "old" "old" "new" "new" ...
#>  $ offset    : chr  "no offset" "no offset" "offset" "offset" ...
#>  $ color     : chr  "color" "b/w" "color" "b/w" ...
#>  $ set.size  : chr  "set size 8" "set size 8" "set size 8" "set size 8" ...
```
