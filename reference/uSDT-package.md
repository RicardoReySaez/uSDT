# uSDT: Hierarchical Signal Detection Theory for Unconscious Processing

Fits hierarchical signal detection theory models to paired direct and
indirect measures. This is the design used to test whether a stimulus is
processed without awareness.

## The workflow

1.  [`usdt_data_long()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_data.md)
    or
    [`usdt_data_tasks()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_data.md)
    prepare the data. Printing the result reports how every column was
    read, which tasks were split at the median, and what the model will
    do with all of it.

2.  [`sdt_moments()`](https://ricardoreysaez.github.io/uSDT/reference/sdt_moments.md)
    gives descriptive estimates for each subject.

3.  [`hsdt()`](https://ricardoreysaez.github.io/uSDT/reference/hsdt.md)
    fits the model and tests the three hypotheses.

4.  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
    [`usdt_reliability()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_reliability.md)
    help to interpret the fit, and
    [`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/reference/usdt_boot.md)
    adds intervals by simulation when the model needs them.

## The three hypotheses

- H1:

  The difference between the average sensitivities of the two tasks.

- H2:

  The correlation between the two sensitivities across subjects.

- H3:

  The regression of the indirect sensitivity on the direct one. Its
  intercept is the sensitivity expected in the indirect task from a
  subject whose direct sensitivity is zero, which is the test for
  unconscious processing.

## Confidence intervals

H1 and the regression of H3 use Wald intervals. The correlation of H2
uses a Fisher-z interval, so its limits stay between -1 and 1. When the
model reaches a boundary and an interval becomes unreliable, the package
reports it as unavailable and explains why.

## See also

Useful links:

- <https://github.com/RicardoReySaez/uSDT>

- <https://ricardoreysaez.github.io/uSDT/>

- Report bugs at <https://github.com/RicardoReySaez/uSDT/issues>

## Author

**Maintainer**: Ricardo Rey-Sáez <ricardoreysaez95@gmail.com>
([ORCID](https://orcid.org/0000-0001-6739-2035))

Authors:

- Ricardo Rey-Sáez <ricardoreysaez95@gmail.com>
  ([ORCID](https://orcid.org/0000-0001-6739-2035))

- Alicia Franco-Martínez <aliciafranco96@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-9710-1240))

- Francisco Garre-Frutos <fgfrutos@gmail.com>
  ([ORCID](https://orcid.org/0000-0001-9810-186X))

- Miguel Vadillo <mgl.vadillo@gmail.com>
  ([ORCID](https://orcid.org/0000-0001-8421-816X))
