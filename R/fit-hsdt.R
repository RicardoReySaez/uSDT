# fit-hsdt.R
# Fit hierarchical Signal Detection Theory models
# Author: Ricardo Rey-Sáez
# Last modified: 03-10-2026

# Public functions

#' Fit a hierarchical signal detection theory model
#'
#' Fits a bivariate hierarchical SDT model and evaluates the core unconscious
#' processing hypotheses. The model estimates task-specific sensitivities
#' (\eqn{d'}) and response criteria (\eqn{c}) for the population, and their
#' variation and correlation across participants. The frequentist estimation
#' uses [lme4::glmer()]; the Bayesian estimation uses Stan.
#'
#' @param data A `usdt_data` object from [usdt_data_tasks()] or
#'   [usdt_data_long()].
#' @param estimation Estimation framework: `"frequentist"` (maximum likelihood
#'   via Laplace approximation with lme4) or `"bayesian"` (Hamiltonian Monte
#'   Carlo with Stan).
#' @param fix_criteria How to handle response criteria. `"auto"` fixes to zero
#'   any criterion that is zero by design (such as a task split at the median
#'   under deviation coding). `"none"` estimates all criteria.
#' @param level Confidence level for Wald intervals, or credible level for the
#'   posterior intervals of a Bayesian fit (default is 0.95).
#' @param optimizer Primary optimizer passed to [lme4::glmerControl()].
#'   Alternative optimizers are automatically evaluated if the default fails to
#'   converge or produces a singular fit. Frequentist estimation only.
#' @param priors A `usdt_priors` object from [usdt_priors()]. Bayesian
#'   estimation only.
#' @param unequal_variances Logical. Estimate the standard deviation of the
#'   signal distribution of each task instead of fixing it to that of the noise
#'   distribution. Bayesian estimation only.
#' @param backend Stan interface used for Bayesian estimation: `"rstan"`
#'   (default, available from CRAN) or `"cmdstanr"` (installed from the Stan
#'   R-universe together with CmdStan). Bayesian estimation only.
#' @param ... For frequentist estimation, additional arguments passed to
#'   [lme4::glmer()]; the model formula, family and data remain managed by the
#'   package. For Bayesian estimation, sampling options: `chains` (default 4),
#'   `iter` (iterations per chain including warmup, default 6000), `warmup`
#'   (default 1000, so each chain keeps 5000 draws), `cores` (chains run at
#'   once; default `min(parallel::detectCores() - 2, 4)`, at least 1), `seed`,
#'   `control = list(adapt_delta, max_treedepth)` (default 0.95 and 10) and
#'   `refresh` (iterations between progress updates of each chain; default 1%
#'   of `iter`, and `0` hides the progress bar).
#'
#' @return An object of class `hsdt` containing:
#' * `$fit`: The underlying fit: a `glmerMod` object from lme4, or the Stan fit
#'   of the chosen backend.
#' * `$tests`: Summary table for hypotheses H1, H2, and H3.
#' * `$design`: Summary of the model specification.
#' * `$diagnostics`: Convergence diagnostics.
#'
#' A frequentist fit also contains `$pars`, the parameter estimates and their
#' joint covariance. A Bayesian fit also contains `$draws`, the posterior draws
#' of the population and subject parameters, and `$priors`.
#'
#' @details
#' The model fits trial counts with a binomial probit link, directly mapping
#' coefficients to standard Signal Detection Theory parameters. Fixed effects
#' capture population sensitivities and criteria, while random effects estimate
#' participant variation and the latent correlation between direct and indirect
#' sensitivity.
#'
#' Hypotheses evaluated by default:
#' * **H1:** Mean sensitivity difference between tasks.
#' * **H2:** Latent correlation of sensitivities across participants.
#' * **H3:** Latent regression of indirect on direct sensitivity. Its intercept
#'   reflects expected indirect performance when direct awareness is zero
#'   (\eqn{d'_{\mathrm{Direct}} = 0}).
#'
#' When sample sizes or trial counts are low, variance components can reach
#' singular boundaries. In these cases, the function issues a warning, and
#' parametric bootstrap intervals can be calculated using [usdt_boot()], or the
#' model can be fitted with `estimation = "bayesian"`.
#'
#' # Bayesian estimation
#'
#' The Bayesian model places its priors on the mean, standard deviation and
#' correlation of the two sensitivities and on the criteria; see
#' [usdt_priors()]. The slope and intercept of the latent regression are
#' derived from them draw by draw. The hypothesis table reports the posterior
#' mean and standard deviation, a central credible interval and a two-sided
#' posterior p-value, \eqn{2 \min\{P(\theta > 0 \mid y), P(\theta < 0 \mid y)\}}.
#' The p-value falls below `1 - level` exactly when zero lies outside the
#' credible interval, so with the default `level = 0.95` a p-value below .05
#' means that the 95% credible interval excludes zero.
#'
#' With `unequal_variances = TRUE`, the noise distribution of each task is
#' \eqn{N(-d/2, 1)} and the signal distribution \eqn{N(d/2, \sigma_s^2)}, so
#' \eqn{d} is the distance between their means in units of the noise standard
#' deviation, and H1-H3 compare these \eqn{d}. When \eqn{\sigma_s} differs
#' between the tasks, \eqn{d} no longer measures the same discriminability in
#' both: the measure comparable with an equal-variance \eqn{d'} is
#' \eqn{d_a = d \sqrt{2 / (1 + \sigma_s^2)}}, and the mean \eqn{d_a} of a task
#' is its mean \eqn{d} times that factor, which can be computed from the draws
#' of `mu_d` and `sigma_s` in `$draws`. The correlation (H2) is the same on
#' either scale.
#'
#' The Stan model is compiled the first time it is used, which needs a C++
#' toolchain, and stays in a per-user cache directory for later sessions. The
#' default sampling settings keep 20,000 posterior draws, from four chains
#' that run in parallel behind a progress bar. A summary after the bar, and a
#' warning when something fails, report divergent transitions, transitions at
#' the maximum tree depth, an R-hat above 1.01 or an effective sample size
#' below 400.
#'
#' @seealso [usdt_data_tasks()], [usdt_priors()], [usdt_tests()],
#'   [usdt_boot()], [plot.hsdt()]
#'
#' @examples
#' # Contextual cuing data from Vadillo et al. (2025)
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
#' # Full summary table with SDT parameters and hypothesis tests
#' summary(m)
#'
#' # Inspect the model formula (indirect criterion omitted by default)
#' m$design$formula
#'
#' \donttest{
#' # Bayesian estimation. The first call compiles the Stan model.
#' if (requireNamespace("rstan", quietly = TRUE)) {
#'   mb <- hsdt(d, estimation = "bayesian", seed = 1)
#'   summary(mb)
#' }
#' }
#'
#' @export
hsdt <- function(data,
                 estimation   = c("frequentist", "bayesian"),
                 fix_criteria = c("auto", "none"),
                 level        = 0.95,
                 optimizer    = "bobyqa",
                 priors       = usdt_priors(),
                 unequal_variances = FALSE,
                 backend      = c("rstan", "cmdstanr"),
                 ...) {

  # The function checks the data and the requested options.
  if (!inherits(data, "usdt_data")) {
    .usdt_stop("`data` must come from usdt_data_long() or usdt_data_tasks(), ",
               "not a plain ", class(data)[1L], ".")
  }
  estimation <- tryCatch(match.arg(estimation), error = function(e)
    .usdt_stop("`estimation` must be `\"frequentist\"` or `\"bayesian\"`."))
  fix_criteria <- match.arg(fix_criteria)
  .check_confidence_level(level)
  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))))) {
    .usdt_stop("every argument in `...` must have a name.")
  }

  # Each estimation takes its own options.
  if (estimation == "bayesian") {
    if (!missing(optimizer)) {
      .usdt_stop("`optimizer` applies only to `estimation = \"frequentist\"`.")
    }
    if (!inherits(priors, "usdt_priors")) {
      .usdt_stop("`priors` must come from usdt_priors(), not a plain ",
                 class(priors)[1L], ".")
    }
    if (!is.logical(unequal_variances) || length(unequal_variances) != 1L ||
        is.na(unequal_variances)) {
      .usdt_stop("`unequal_variances` must be `TRUE` or `FALSE`.")
    }
    backend <- tryCatch(match.arg(backend), error = function(e)
      .usdt_stop("`backend` must be `\"rstan\"` or `\"cmdstanr\"`."))
    return(.hsdt_bayes(data, fix_criteria, level, priors, unequal_variances,
                       backend, dots, match.call()))
  }
  if (!missing(priors) || !missing(unequal_variances) || !missing(backend)) {
    .usdt_stop("`priors`, `unequal_variances` and `backend` apply only to ",
               "`estimation = \"bayesian\"`.")
  }
  if (!is.character(optimizer) || length(optimizer) != 1L ||
      is.na(optimizer) || !nzchar(optimizer)) {
    .usdt_stop("`optimizer` must be one optimizer name.")
  }
  removed <- intersect(names(dots), c("re", "ci"))
  if (length(removed)) {
    .usdt_stop("`", removed[1L], "` is not a uSDT option. The model always ",
               "correlates the two sensitivities and uses its fixed interval ",
               "methods.")
  }
  reserved <- intersect(names(dots), c("formula", "data", "family", "control",
                                       "devFunOnly"))
  if (length(reserved)) {
    .usdt_stop("uSDT controls `", paste(reserved, collapse = "`, `"),
               "`, so it cannot be supplied through `...`.")
  }
  if ("nAGQ" %in% names(dots)) {
    nAGQ <- dots$nAGQ
    if (!is.numeric(nAGQ) || length(nAGQ) != 1L || is.na(nAGQ) || nAGQ != 1) {
      .usdt_stop("`nAGQ` must be 1. uSDT requires the Laplace fit for its ",
                 "joint inference.")
    }
    dots$nAGQ <- NULL
  }

  # The selected settings define the model formula.
  design <- .usdt_formula(data, fix_criteria = fix_criteria)
  agg    <- data$agg

  # The formula keeps the aggregated data available.
  fenv <- new.env(parent = parent.frame())
  assign("agg", agg, envir = fenv)
  environment(design$formula) <- fenv

  # The code fits the model and may try another optimizer.
  fitted <- .fit_retry(design$formula, agg, optimizer, dots)
  fit    <- fitted$fit

  # The fitted settings recreate the deviance function.
  devfun <- tryCatch(
    .glmer_fit(design$formula, agg, fitted$optimizer, dots, devfun = TRUE),
    error = function(e) NULL)

  # The fitted parameters produce the three hypothesis tests.
  pars  <- .usdt_pars(fit, "d_D", "d_I", devfun = devfun)
  tests <- rbind(
    cbind(hypothesis = "H1", .diff_rows(pars, level)),
    cbind(hypothesis = "H2", .cor_rows(pars, level)),
    cbind(hypothesis = "H3", .reg_rows(pars, level))
  )
  rownames(tests) <- NULL

  # A warning explains when the fit cannot support these intervals.
  diag <- .diagnose(fit, pars, fitted)
  if (!pars$joint_ok) {
    .usdt_warn("H2 and H3 are not estimable with reliable standard errors.\n",
               if (length(pars$at_bound))
                 paste0("  At the bound: ",
                        paste(pars$at_bound, collapse = ", "), ".\n") else "",
               "  Reason: ", pars$inference_reason, "\n",
               "  Use usdt_boot() for a parametric bootstrap, or ",
               "estimation = \"bayesian\".")
  }

  structure(list(fit = fit, tests = tests, pars = pars, design = design,
                 data = data, diagnostics = diag, devfun = devfun,
                 call = match.call(), level = level,
                 estimation = "frequentist"),
            class = "hsdt")
}

# Internal functions

# This function sends the same arguments to every model fit.
.glmer_fit <- function(formula, agg, optimizer, dots, devfun = FALSE) {
  args <- c(list(formula = formula, data = agg,
                 family = stats::binomial("probit"),
                 control = lme4::glmerControl(optimizer = optimizer),
                 nAGQ = 1L), dots)
  if (devfun) args$devFunOnly <- TRUE
  suppressMessages(do.call(lme4::glmer, args))
}

# This function reads the convergence result stored by lme4.
.fit_convergence <- function(fit) {
  opt <- fit@optinfo$conv$opt
  opt_ok <- !length(opt) || (is.numeric(opt) && all(opt == 0))
  messages <- fit@optinfo$conv$lme4$messages
  if (length(messages)) {
    messages <- messages[!grepl("boundary.*singular", messages,
                                ignore.case = TRUE)]
  }
  list(ok = opt_ok && !length(messages),
       message = paste(c(if (!opt_ok) paste("optimizer code", paste(opt, collapse = ", ")),
                         messages), collapse = "; "))
}

# This function tries several optimizers until one converges.
.fit_retry <- function(formula, agg, optimizer, dots) {
  tried <- character(0)
  reasons <- character(0)
  queue <- unique(c(optimizer, "Nelder_Mead", "nlminbwrap"))
  singular_fit <- NULL

  for (opt in queue) {
    failure <- NULL
    fit <- tryCatch(.glmer_fit(formula, agg, opt, dots),
                    error = function(e) {
                      failure <<- conditionMessage(e)
                      NULL
                    })
    tried <- c(tried, opt)
    if (is.null(fit)) {
      reasons <- c(reasons, paste0(opt, ": ", failure))
      next
    }
    conv <- .fit_convergence(fit)
    if (!conv$ok) {
      reasons <- c(reasons, paste0(opt, ": ", conv$message))
      next
    }
    candidate <- list(fit = fit, optimizer = opt, tried = tried,
                      retried = length(tried) > 1L)
    if (!isTRUE(lme4::isSingular(fit, tol = 1e-4))) {
      return(candidate)
    }
    if (is.null(singular_fit)) {
      singular_fit <- candidate
    }
  }

  # A converged boundary fit can still provide point estimates.
  if (!is.null(singular_fit)) {
    singular_fit$tried <- tried
    singular_fit$retried <- length(tried) > 1L
    return(singular_fit)
  }
  detail <- unique(reasons[nzchar(reasons)])
  .usdt_stop("the model did not converge with any available optimizer.",
             if (length(detail)) paste0("\n  ", paste(detail, collapse = "\n  ")) else "",
             "\n  Check the data and model identification, or fit the model ",
             "with estimation = \"bayesian\".")
}

# This function collects the main fitting diagnostics.
.diagnose <- function(fit, pars, fitted) {
  derivs <- fit@optinfo$derivs
  grad <- tryCatch(
    max(abs(solve(derivs$Hessian, derivs$gradient))),
    error = function(e) NA_real_)
  sds <- sqrt(c(pars$est[["s2_D"]], pars$est[["s2_I"]]))
  list(optimizer = fitted$optimizer, retried = fitted$retried,
       tried = fitted$tried, max_grad = grad,
       singular = pars$singular, boundary = !pars$joint_ok,
       min_sd = min(sds), rho = pars$rho,
       converged = .fit_convergence(fit)$ok,
       inference_reason = pars$inference_reason)
}
