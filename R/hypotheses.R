# hypotheses.R
# This script calculates the three hypothesis tests used by uSDT.
# Author: Ricardo Rey-Sáez
# Last modified: 04-09-2026

# Public functions

#' Test the three hypotheses of a hierarchical SDT model
#'
#' Computes the group-level difference between sensitivities, their latent
#' correlation, and the latent regression of the indirect sensitivity on the
#' direct one. They accept any `glmerMod` in which the two sensitivities are
#' fixed effects and share a random-effects term, not only models built by
#' [hsdt()].
#'
#' @param fit A fitted `glmerMod`, typically from [hsdt()] or from
#'   `lme4::glmer()` directly.
#' @param direct,indirect Names of the two sensitivity terms in the model.
#' @param level Confidence level.
#'
#' @return A data frame with one row per quantity and columns `term`,
#'   `estimate`, `se`, `statistic`, `p.value`, `conf.low`, `conf.high` and
#'   `ci_method`. The columns `status` and `reason` identify results that cannot
#'   support inference. `usdt_tests()` adds a `hypothesis` column.
#'
#' @details
#' H1 and the two regression terms use Wald intervals. The latent correlation
#' uses a Fisher-z interval so its limits remain between -1 and 1.
#'
#' The latent slope is zero exactly when the latent covariance is zero, which
#' is also when the latent correlation is zero, so H2 and the H3 slope are the
#' same null hypothesis. Both are therefore reported with the same Wald test on
#' the covariance.
#'
#' @seealso [hsdt()]
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' df <- usdt_simulate(n_subj = 40, n_trials = 100)
#' d  <- usdt_data_long(df, task_col = "task",
#'                      task_levels      = c(direct = "D", indirect = "I"),
#'                      subject_col      = "subj",
#'                      condition_col    = "cond",
#'                      condition_levels = c(signal = 1, noise = 0),
#'                      response_col     = "response",
#'                      response_levels  = c(signal = 1, noise = 0))
#' m <- hsdt(d)
#' usdt_tests(m$fit)
#' }
#'
#' @name usdt_hypotheses
NULL

#' @rdname usdt_hypotheses
#' @export
usdt_tests <- function(fit, direct = "d_D", indirect = "d_I", level = 0.95) {

  # The function collects the values shared by the three tests.
  .check_confidence_level(level)
  r  <- .resolve_fit(fit)
  p  <- .usdt_pars(r$fit, direct, indirect, devfun = r$devfun)

  # The function calculates the three tests from those values.
  out <- rbind(
    cbind(hypothesis = "H1", .diff_rows(p, level)),
    cbind(hypothesis = "H2", .cor_rows(p, level)),
    cbind(hypothesis = "H3", .reg_rows(p, level))
  )
  rownames(out) <- NULL
  attr(out, "boundary") <- p$boundary
  out
}

#' @rdname usdt_hypotheses
#' @export
sensitivity_diff <- function(fit, direct = "d_D", indirect = "d_I",
                             level = 0.95) {
  .check_confidence_level(level)
  .diff_rows(.pars_of(fit, direct, indirect), level)
}

#' @rdname usdt_hypotheses
#' @export
latent_cor <- function(fit, direct = "d_D", indirect = "d_I",
                       level = 0.95) {
  .check_confidence_level(level)
  .cor_rows(.pars_of(fit, direct, indirect), level)
}

#' @rdname usdt_hypotheses
#' @export
latent_regression <- function(fit, direct = "d_D", indirect = "d_I",
                              level = 0.95) {
  .check_confidence_level(level)
  .reg_rows(.pars_of(fit, direct, indirect), level)
}

# Internal calculations

# This function tests the group-level sensitivity difference.
.diff_rows <- function(p, level) {
  e <- p$est
  est  <- e[["gamma_D"]] - e[["gamma_I"]]
  V <- p$fixed_vcov
  if (is.null(V) || !all(p$fixed_names %in% rownames(V))) {
    return(.row_na("d'(direct) - d'(indirect)", est,
                   "the fixed-effects covariance is unavailable"))
  }
  V <- V[p$fixed_names, p$fixed_names, drop = FALSE]
  se <- .delta_se(c(1, -1), V)
  if (is.na(se)) {
    return(.row_na("d'(direct) - d'(indirect)", est,
                   "the fixed-effects standard error is invalid"))
  }
  w <- .wald(est, se, level)
  .row("d'(direct) - d'(indirect)", est, se, w$statistic, w$p.value,
       w$conf.low, w$conf.high, w$ci_method)
}

# This function tests the latent sensitivity correlation.
.cor_rows <- function(p, level) {
  e   <- p$est
  rho <- p$rho
  if (!p$joint_ok) {
    return(.row_na("correlation", rho, p$inference_reason))
  }

  # The gradient gives the uncertainty of the correlation.
  grad <- c(0, 0, -rho / (2 * e[["s2_D"]]), -rho / (2 * e[["s2_I"]]),
            1 / sqrt(e[["s2_D"]] * e[["s2_I"]]))
  se   <- .delta_se(grad, p$vcov)
  if (is.na(se)) {
    return(.row_na("correlation", rho,
                   "the correlation standard error is invalid"))
  }

  # Fisher-z keeps the interval inside the correlation range.
  int <- .fisher_z(rho, se, level)

  # The covariance provides the test for a null correlation.
  tst <- .cov_test(p)
  .row("correlation", rho, se, tst$statistic, tst$p.value,
       int$conf.low, int$conf.high, int$ci_method)
}

# This function tests the latent sensitivity regression.
.reg_rows <- function(p, level) {
  e <- p$est
  s2D <- e[["s2_D"]]; sDI <- e[["s_DI"]]
  gD  <- e[["gamma_D"]]; gI <- e[["gamma_I"]]
  slope <- sDI / s2D
  inter <- gI - slope * gD
  if (!p$joint_ok) {
    return(rbind(.row_na("intercept", inter, p$inference_reason),
                 .row_na("slope", slope, p$inference_reason)))
  }

  # These gradients give the uncertainty of both regression terms.
  g_slope <- c(0, 0, -sDI / s2D^2, 0, 1 / s2D)
  se_slope <- .delta_se(g_slope, p$vcov)
  se_inter <- .delta_se(.reg_grad(e), p$vcov)

  # The intercept uses a Wald interval on its own scale.
  row_i <- if (is.na(se_inter)) {
    .row_na("intercept", inter, "the intercept standard error is invalid")
  } else {
    w <- .wald(inter, se_inter, level)
    .row("intercept", inter, se_inter, w$statistic, w$p.value,
         w$conf.low, w$conf.high, "delta")
  }

  # The slope uses a Wald interval on the regression scale.
  if (is.na(se_slope)) {
    row_s <- .row_na("slope", slope, "the slope standard error is invalid")
  } else {
    int <- .wald(slope, se_slope, level)
    tst <- .cov_test(p)
    row_s <- .row("slope", slope, se_slope, tst$statistic, tst$p.value,
                  int$conf.low, int$conf.high, int$ci_method)
  }
  rbind(row_i, row_s)
}

# This function gives one latent-line gradient per direct value.
.reg_grad <- function(e, x = 0) {
  s2D <- e[["s2_D"]]; sDI <- e[["s_DI"]]; gD <- e[["gamma_D"]]
  cbind(-sDI / s2D, 1, (gD - x) * sDI / s2D^2, 0, (x - gD) / s2D)
}

# This function tests the covariance shared by H2 and H3.
.cov_test <- function(p) {
  se <- sqrt(p$vcov[["s_DI", "s_DI"]])
  if (!is.finite(se) || se <= 0) {
    return(list(statistic = NA_real_, p.value = NA_real_))
  }
  st <- p$est[["s_DI"]] / se
  list(statistic = st, p.value = 2 * stats::pnorm(-abs(st)))
}

# Model input

# This function finds the fitted model and its deviance function.
.resolve_fit <- function(x) {
  if (inherits(x, "hsdt")) return(list(fit = x$fit, devfun = x$devfun))
  list(fit = x, devfun = NULL)
}

# This function collects the five values used by the tests.
.pars_of <- function(x, direct, indirect) {
  r <- .resolve_fit(x)
  .usdt_pars(r$fit, direct, indirect, devfun = r$devfun)
}
