# Cuing data from a probabilistic cuing experiment

Trial-level indirect-task data from Experiment 2 of Vadillo et al.
(2025). The object reproduces the source CSV file without filtering or
recoding and can be used as the `indirect` input to
[`usdt_data_tasks()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_data.md).

## Usage

``` r
vadillo_cuing
```

## Format

A data frame with 39,936 rows, 104 participants, and 13 variables:

- `experiment`:

  Experiment label.

- `subj`:

  Participant identifier.

- `file`:

  Original participant file name.

- `trial`:

  Trial number.

- `block`:

  Block number.

- `epoch`:

  Epoch number.

- `pattId`:

  Pattern identifier.

- `acc`:

  Response accuracy (1 = correct, 0 = incorrect).

- `rt`:

  Response time in milliseconds.

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
data(vadillo_cuing)
str(vadillo_cuing)
#> 'data.frame':    39936 obs. of  13 variables:
#>  $ experiment: chr  "Experiment 2" "Experiment 2" "Experiment 2" "Experiment 2" ...
#>  $ subj      : int  2001 2001 2001 2001 2001 2001 2001 2001 2001 2001 ...
#>  $ file      : chr  "y109b_subj1051.mat" "y109b_subj1051.mat" "y109b_subj1051.mat" "y109b_subj1051.mat" ...
#>  $ trial     : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ block     : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ epoch     : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ pattId    : int  2001 1008 1007 2006 1001 1006 2003 1003 1005 2002 ...
#>  $ acc       : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ rt        : num  1305 4028 2484 1771 926 ...
#>  $ condition : chr  "new" "old" "old" "new" ...
#>  $ offset    : chr  "no offset" "offset" "offset" "offset" ...
#>  $ color     : chr  "color" "b/w" "color" "b/w" ...
#>  $ set.size  : chr  "set size 8" "set size 16" "set size 16" "set size 8" ...
```
