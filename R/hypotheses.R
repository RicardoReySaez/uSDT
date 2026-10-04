# hypotheses.R
# Hypothesis tests for hierarchical SDT models
# Author: Ricardo Rey-Sáez
# Last modified: 03-10-2026

# Public functions

#' Test the three core hypotheses of a hierarchical SDT model
#'
#' Evaluates the difference between average task sensitivities (H1), their
#' correlation across subjects (H2), and the regression of indirect sensitivity
#' on direct sensitivity (H3). `usdt_tests()` computes all three together,
#' while individual functions compute them separately.
#'
#' These tests run automatically inside [hsdt()] and appear in its summary.
#' Calling them directly is especially useful when fitting custom models with
#' [lme4::glmer()] or [brms::brm()], allowing you to test these hypotheses
#' while controlling for additional covariates (e.g., set size, experimental
#' groups).
#'
#' @param fit A fitted model: an `hsdt` object, frequentist or Bayesian, a
#'   `glmerMod` from `lme4::glmer()` or a `brmsfit` from `brms::brm()`.
#' @param direct,indirect Character strings naming the sensitivity terms in the
#'   model. Defaults match the internal naming of [hsdt()]. For custom models,
#'   both terms must be fixed effects and share a common random-effects grouping
#'   by subject.
#' @param level Confidence level for intervals (default is 0.95).
#'
#' @return A data frame with columns `term`, `estimate`, `se`, `statistic`,
#'   `p.value`, `conf.low`, `conf.high`, and `ci_method`. Columns `status` and
#'   `reason` flag unsupported estimates (e.g., singular fits). `usdt_tests()`
#'   includes an extra `hypothesis` column (`H1`, `H2`, `H3`).
#'
#'   For a Bayesian `hsdt` fit the columns are `term`, `estimate` (posterior
#'   mean), `est.error` (posterior standard deviation), `conf.low` and
#'   `conf.high` (central credible interval), `ci_method`, `p.value` (two-sided
#'   posterior p-value, below `1 - level` exactly when the interval excludes
#'   zero), `rhat` and `ess_bulk`.
#'
#' @details
#' # The three hypotheses
#'
#' * **H1 (Mean difference):** Tests whether average sensitivity differs
#'   between the direct and indirect tasks.
#' * **H2 (Correlation):** Tests the correlation between task sensitivities
#'   across participants using a Fisher-\eqn{z} transformed interval.
#' * **H3 (Latent regression):** Regresses indirect sensitivity onto direct
#'   sensitivity. The intercept represents expected indirect performance when
#'   direct sensitivity is zero (\eqn{d'_{\mathrm{Direct}} = 0}), testing for
#'   unconscious processing.
#'
#' Because both the correlation (H2) and regression slope (H3) are zero if and
#' only if the covariance between sensitivities is zero, they evaluate the same
#' association and share identical test statistics.
#'
#' # Custom models with covariates
#'
#' To adjust tests for additional factors, specify the model directly using
#' [lme4::glmer()] or [brms::brm()]. As long as the two sensitivity terms are
#' included as fixed effects and correlated across subjects via random slopes,
#' `usdt_tests()` will compute the latent tests conditional on those
#' covariates. A brms model is summarised from its posterior, like a Bayesian
#' [hsdt()] fit, and keeps the term names that lme4 would give it. Aggregating
#' the trials with [usdt_aggregate()] gives either model the same estimates
#' and makes it much faster to fit.
#'
#' @seealso [hsdt()], [usdt_boot()], [usdt_aggregate()]
#'
#' @examples
#' # 1. Standard model via hsdt()
#' d <- usdt_data_tasks(
#'   direct   = vadillo_awareness,
#'   indirect = vadillo_cuing,
#'   subject_col      = "subj",
#'   condition_col    = "condition",
#'   condition_levels = c(signal = "old", noise = "new"),
#'   response_col     = list(direct = "judged.old", indirect = "rt"),
#'   response_levels  = list(direct   = c(signal = 1, noise = 0),
#'                           indirect = c(signal = "faster", noise = "slower")),
#'   dichotomize      = list(direct = FALSE, indirect = TRUE)
#' )
#'
#' m <- hsdt(d)
#'
#' # All three tests at once
#' usdt_tests(m)
#'
#' # Or one test at a time
#' sensitivity_diff(m)   # H1
#' latent_cor(m)         # H2
#' latent_regression(m)  # H3
#'
#' \donttest{
#' # 2. Custom model with covariates via glmer()
#' # Controlling for display set size across both tasks
#' trials <- rbind(
#'   data.frame(vadillo_awareness[c("subj", "condition", "set.size")],
#'              task = "D", resp = vadillo_awareness$judged.old),
#'   data.frame(vadillo_cuing[c("subj", "condition", "set.size")],
#'              task = "I", resp = meyen_split(vadillo_cuing$rt,
#'                                             by = vadillo_cuing$subj))
#' )
#'
#' # Deviation contrasts (-0.5 vs 0.5); `direct` flags the direct task
#' trials <- within(trials, {
#'   cond   <- ifelse(condition == "old", 0.5, -0.5)
#'   size   <- ifelse(set.size == "set size 16", 0.5, -0.5)
#'   direct <- as.integer(task == "D")
#' })
#'
#' # Counting responses per subject, task and condition keeps the likelihood
#' # and makes the fit much faster than one row per trial
#' counts <- usdt_aggregate(trials, response = "resp",
#'                          by = c("subj", "task", "cond", "size", "direct"))
#'
#' # Standard glmer formula: indirect criterion is omitted (fixed at 0
#' # by the median split). Random effects estimate the direct criterion
#' # and correlated task sensitivities across subjects.
#' fit <- lme4::glmer(
#'   cbind(y, n - y) ~ 0 + direct + task:size + task:cond +
#'     (0 + direct | subj) + (0 + task:cond | subj),
#'   data = counts, family = binomial("probit"),
#'   control = lme4::glmerControl(optimizer = "bobyqa")
#' )
#'
#' # Check the names lme4 assigned to the sensitivity terms
#' names(lme4::fixef(fit))
#'
#' # Evaluate hypotheses conditional on set size
#' usdt_tests(fit, direct = "taskD:cond", indirect = "taskI:cond")
#' }
#'
#' \donttest{
#' # 3. The same model fitted with brms. Compiling it needs a C++ toolchain.
#' if (requireNamespace("brms", quietly = TRUE)) {
#'   fit_b <- brms::brm(
#'     y | trials(n) ~ 0 + direct + task:size + task:cond +
#'       (0 + direct | subj) + (0 + task:cond | subj),
#'     data = counts, family = binomial("probit"),
#'     prior = brms::prior(normal(0, 1), class = b), seed = 1, refresh = 0
#'   )
#'
#'   # The same term names; the tests now summarise the posterior
#'   usdt_tests(fit_b, direct = "taskD:cond", indirect = "taskI:cond")
#' }
#' }
#'
#' @name usdt_hypotheses
NULL

#' @rdname usdt_hypotheses
#' @export
usdt_tests <- function(fit, direct = "d_D", indirect = "d_I", level = 0.95) {

  # The function collects the values shared by the three tests.
  .check_confidence_level(level)
  if (.is_posterior(fit)) {
    return(.bayes_tests(.bivariate_draws(fit, direct, indirect), level))
  }
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
  if (.is_posterior(fit)) {
    return(.bayes_rows(fit, direct, indirect, "H1", level))
  }
  .diff_rows(.pars_of(fit, direct, indirect), level)
}

#' @rdname usdt_hypotheses
#' @export
latent_cor <- function(fit, direct = "d_D", indirect = "d_I",
                       level = 0.95) {
  .check_confidence_level(level)
  if (.is_posterior(fit)) {
    return(.bayes_rows(fit, direct, indirect, "H2", level))
  }
  .cor_rows(.pars_of(fit, direct, indirect), level)
}

#' @rdname usdt_hypotheses
#' @export
latent_regression <- function(fit, direct = "d_D", indirect = "d_I",
                              level = 0.95) {
  .check_confidence_level(level)
  if (.is_posterior(fit)) {
    return(.bayes_rows(fit, direct, indirect, "H3", level))
  }
  .reg_rows(.pars_of(fit, direct, indirect), level)
}

# Internal calculations

# This function tests the group-level sensitivity difference.
.diff_rows <- function(p, level) {
  est <- .usdt_quantities(p$est)[, "diff"]
  V <- p$fixed_vcov
  if (is.null(V) || !all(p$fixed_names %in% rownames(V))) {
    return(.row_na("d'(indirect) - d'(direct)", est,
                   "the fixed-effects covariance is unavailable"))
  }
  V <- V[p$fixed_names, p$fixed_names, drop = FALSE]
  se <- .delta_se(c(-1, 1), V)
  if (is.na(se)) {
    return(.row_na("d'(indirect) - d'(direct)", est,
                   "the fixed-effects standard error is invalid"))
  }
  w <- .wald(est, se, level)
  .row("d'(indirect) - d'(direct)", est, se, w$statistic, w$p.value,
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
  q <- .usdt_quantities(e)
  slope <- q[, "slope"]
  inter <- q[, "intercept"]
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

# This function tells a Bayesian hsdt fit from every frequentist input.
.is_bayes <- function(x) {
  inherits(x, "hsdt") && identical(x$estimation, "bayesian")
}

# This function tells a posterior, from a Bayesian hsdt fit or brms, from a
# maximum-likelihood fit.
.is_posterior <- function(x) {
  .is_bayes(x) || inherits(x, "brmsfit")
}

# This function returns the posterior rows of one hypothesis.
.bayes_rows <- function(fit, direct, indirect, hypothesis, level) {
  tests <- .bayes_tests(.bivariate_draws(fit, direct, indirect), level)
  out <- tests[tests$hypothesis == hypothesis, -1L, drop = FALSE]
  rownames(out) <- NULL
  out
}

# This function draws the bivariate normal of the two sensitivities from a
# posterior. A Bayesian hsdt fit names it itself; a brms model names it after
# the sensitivity terms.
.bivariate_draws <- function(fit, direct, indirect) {
  if (.is_bayes(fit)) return(.stan_bivariate(fit$draws))
  rlang::check_installed(c("brms", "posterior"), reason = "to test a brms model.")
  variables <- .brms_variables(posterior::variables(fit), direct, indirect)
  B <- posterior::as_draws_array(fit, variable = variables)
  .select_bivariate(B, variables)
}

# This function finds, among the variable names of a brms model, the means,
# SDs and correlation of the two sensitivity terms: `b_<term>`,
# `sd_<group>__<term>` and `cor_<group>__<term>__<term>` in either order. The
# names are matched as plain text, because terms such as `taskD:cond` hold
# characters that regular expressions reserve.
.brms_variables <- function(variables, direct, indirect) {
  fixed <- paste0("b_", c(direct, indirect))
  missing <- c(direct, indirect)[!fixed %in% variables]
  if (length(missing)) {
    available <- substring(variables[startsWith(variables, "b_")], 3L)
    .usdt_stop("population-level effect", if (length(missing) > 1L) "s" else "",
               " ", paste0("`", missing, "`", collapse = " and "),
               " not found in the model.\n  Available: ",
               paste(available, collapse = ", "))
  }

  # The groups whose SDs include a term.
  sds <- variables[startsWith(variables, "sd_")]
  groups_of <- function(term) {
    suffix <- paste0("__", term)
    hit <- sds[endsWith(sds, suffix)]
    substring(hit, 4L, nchar(hit) - nchar(suffix))
  }
  for (group in intersect(groups_of(direct), groups_of(indirect))) {
    pairs <- paste0("cor_", group, "__",
                    c(paste0(direct, "__", indirect),
                      paste0(indirect, "__", direct)))
    correlation <- pairs[pairs %in% variables]
    if (length(correlation)) {
      return(c(fixed, paste0("sd_", group, "__", c(direct, indirect)),
               correlation[1L]))
    }
  }
  .usdt_stop("no group-level term correlates `", direct, "` and `", indirect,
             "`.\n  The latent correlation requires the two sensitivities to ",
             "vary together in the same term, e.g. (0 + ", direct, " + ",
             indirect, " | subject).")
}

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
