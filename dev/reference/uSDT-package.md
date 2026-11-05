# uSDT: Hierarchical Signal Detection Theory for Unconscious Processing

Fits hierarchical signal detection theory models to paired awareness and
influence measures. This is the design used to test whether a stimulus
is processed without awareness.

## The workflow

1.  [`usdt_data_long()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_data.md)
    or
    [`usdt_data_tasks()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_data.md)
    prepare the data. Printing the result reports how every column was
    read, which tasks were split at the median, and what the model will
    do with all of it.

2.  [`sdt_moments()`](https://ricardoreysaez.github.io/uSDT/dev/reference/sdt_moments.md)
    gives descriptive estimates for each subject.

3.  [`hsdt()`](https://ricardoreysaez.github.io/uSDT/dev/reference/hsdt.md)
    fits the model and tests the three hypotheses, by maximum likelihood
    with 'lme4' or, with `estimation = "bayesian"`, with 'brms' under
    the priors of
    [`usdt_priors()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_priors.md).

4.  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
    [`usdt_reliability()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_reliability.md)
    help to interpret the fit,
    [`usdt_boot()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_boot.md)
    adds intervals by simulation when the maximum likelihood model needs
    them, and
    [`usdt_bf()`](https://ricardoreysaez.github.io/uSDT/dev/reference/usdt_bf.md)
    weighs each hypothesis of a Bayesian fit with Bayes factors.

## The three hypotheses

- H1:

  The difference between the average sensitivities of the two tasks.

- H2:

  The correlation between the two sensitivities across subjects.

- H3:

  The regression of the influence sensitivity on the awareness one. Its
  intercept is the sensitivity expected in the influence task from a
  subject whose awareness sensitivity is zero, which is the test for
  unconscious processing.

## Intervals

In a maximum likelihood fit, H1 and the regression of H3 use Wald
intervals. The correlation of H2 uses a Fisher-z interval, so its limits
stay between -1 and 1. When the model reaches a boundary and an interval
becomes unreliable, the package reports it as unavailable and explains
why. A Bayesian fit reports central credible intervals from the
posterior draws, which need no such approximation.

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

- Francisco Garre-Frutos <fgfrutos@gmail.com>
  ([ORCID](https://orcid.org/0000-0001-9810-186X))

- Alicia Franco-Martínez <aliciafranco96@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-9710-1240))

- Miguel Vadillo <mgl.vadillo@gmail.com>
  ([ORCID](https://orcid.org/0000-0001-8421-816X))
