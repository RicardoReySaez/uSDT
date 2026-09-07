# uSDT-package.R
# This script documents the package and records its imports.
# Author: Ricardo Rey-Sáez
# Last modified: 04-09-2026

#' uSDT: Hierarchical Signal Detection Theory for Unconscious Processing
#'
#' Fits hierarchical signal detection theory models to paired direct and
#' indirect measures, the design used to test whether a stimulus is processed
#' without awareness.
#'
#' @section The workflow:
#' \enumerate{
#'   \item [usdt_data_long()] or [usdt_data_tasks()] prepare the data. They say
#'     out loud how every column was read, what the median split did, and what
#'     all of it implies for the model.
#'   \item [usdt_reliability()] estimates the reliability of both measures.
#'   \item [hsdt()] fits the model and tests the three hypotheses.
#'   \item [usdt_boot()] replaces the intervals with a parametric bootstrap when
#'     the model is weakly identified.
#' }
#'
#' @section The three hypotheses:
#' \describe{
#'   \item{H1}{The group-level difference between the two sensitivities.}
#'   \item{H2}{Their latent correlation across subjects.}
#'   \item{H3}{The latent regression of the indirect sensitivity on the direct
#'     one. Its intercept is the indirect sensitivity expected of a subject
#'     whose direct sensitivity is zero, which is the test for unconscious
#'     processing.}
#' }
#'
#' @section Confidence intervals:
#' The group difference and latent regression use Wald intervals. The latent
#' correlation uses Fisher-z so its limits remain between -1 and 1. The package
#' reports unavailable intervals explicitly when the model reaches a boundary.
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
