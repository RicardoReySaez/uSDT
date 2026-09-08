# uSDT-package.R
# This script documents the package and records its imports.
# Author: Ricardo Rey-Sáez
# Last modified: 08-09-2026

#' uSDT: Hierarchical Signal Detection Theory for Unconscious Processing
#'
#' Fits hierarchical signal detection theory models to paired direct and
#' indirect measures. This is the design used to test whether a stimulus is
#' processed without awareness.
#'
#' @section The workflow:
#' \enumerate{
#'   \item [usdt_data_long()] or [usdt_data_tasks()] prepare the data. Printing
#'     the result reports how every column was read, which tasks were split at
#'     the median, and what the model will do with all of it.
#'   \item [sdt_moments()] gives descriptive estimates for each subject.
#'   \item [hsdt()] fits the model and tests the three hypotheses.
#'   \item [plot()] and [usdt_reliability()] help to interpret the fit, and
#'     [usdt_boot()] adds intervals by simulation when the model needs them.
#' }
#'
#' @section The three hypotheses:
#' \describe{
#'   \item{H1}{The difference between the average sensitivities of the two
#'     tasks.}
#'   \item{H2}{The correlation between the two sensitivities across subjects.}
#'   \item{H3}{The regression of the indirect sensitivity on the direct one.
#'     Its intercept is the sensitivity expected in the indirect task from a
#'     subject whose direct sensitivity is zero, which is the test for
#'     unconscious processing.}
#' }
#'
#' @section Confidence intervals:
#' H1 and the regression of H3 use Wald intervals. The correlation of H2 uses a
#' Fisher-z interval, so its limits stay between -1 and 1. When the model
#' reaches a boundary and an interval becomes unreliable, the package reports it
#' as unavailable and explains why.
#'
#' @keywords internal
"_PACKAGE"

# Imports

#' @importFrom stats qnorm pnorm rnorm rbinom median sd ave setNames
#' @importFrom stats aggregate as.formula binomial complete.cases formula
#' @importFrom stats getCall model.frame quantile reshape update vcov family
#' @importFrom utils head tail modifyList
#' @importFrom rlang .data
NULL
