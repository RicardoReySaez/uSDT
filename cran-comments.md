# CRAN comments for uSDT 0.2.0

## Update

This updates uSDT 0.1.0 with Bayesian estimation through brms, Bayes factors,
equivalence tests and study planning. There are no known reverse dependencies
on CRAN as of 2026-11-05.

## Previous review

* The bootstrap example remains in \donttest{} because it refits 500 models.
* Examples measured below five seconds per help page are unwrapped, including
  the aggregated lme4 example with covariates.
* Only examples requiring an external C++ toolchain are in \dontrun{}.
  Their comments state the compiler and rstan header requirements.
* Package code no longer assigns or removes .Random.seed directly. Scoped
  withr seeds restore the user's RNG state and generator settings, including
  after errors. Tests check reproducibility across one and two cores.

## Bayesian requirements

brms and the default rstan backend are available from CRAN. Compiling a
Bayesian model requires a C++17 toolchain and the BH and RcppEigen headers.
No model is compiled during package installation. The optional cmdstanr
backend is in Suggests and its repository is declared in Additional_repositories.

Compiled models use the package-specific tools::R_user_dir() cache. The cache
is bounded and outdated entries are removed. Tests use a temporary cache.
The Bayesian vignette is precomputed, with reproducible source in the GitHub
repository. Rebuilding the packaged vignettes does not require Stan.

## Checks

The validation report and logs are retained with the release candidate.
The final submission comments will list the completed checks on this tarball.

## Data

The example data are unchanged from 0.1.0. Their rights holders have authorised
redistribution under GPL-3. Source attribution is retained in the help pages
and inst/COPYRIGHTS.
