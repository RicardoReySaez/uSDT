# fit-freq.R
# This script fits the frequentist hierarchical uSDT model.
# Author: Ricardo Rey-Sáez
# Last modified: 04-09-2026

# Public functions

#' Fit a hierarchical SDT model by maximum likelihood
#'
#' Fits the binomial probit mixed model in which the sensitivities of the
#' direct and indirect tasks are correlated random effects, and tests the three
#' hypotheses of the unconscious-processing design.
#'
#' @param data A `usdt_data` object from [usdt_data_long()] or
#'   [usdt_data_tasks()].
#' @param fix_criteria `"auto"` fixes to zero every criterion the data show to
#'   be zero by construction, which a Meyen median split under deviation coding
#'   guarantees. `"none"` estimates them all.
#' @param level Confidence level.
#' @param optimizer Optimizer passed to `lme4::glmerControl()`. Other optimizers
#'   are tried when the model does not converge or reaches a singular fit.
#' @param ... Named arguments passed unchanged to every `lme4::glmer()` call.
#'   The formula, data, family, optimizer and `nAGQ = 1` remain fixed by uSDT.
#'   Print and summary methods ignore this argument.
#'
#' @return An object of class `usdt_freq`: a list with the fitted model
#'   (`fit`), the hypothesis table (`tests`), the extracted parameters
#'   (`pars`), the formula information (`design`) and the diagnostics
#'   (`diagnostics`).
#'
#' @seealso [usdt_data_long()], [usdt_tests()], [plot.usdt_freq()]
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
#' m <- usdt_freq(d)
#' summary(m)
#' }
#'
#' @export
usdt_freq <- function(data,
                      fix_criteria = c("auto", "none"),
                      level        = 0.95,
                      optimizer    = "bobyqa",
                      ...) {

  # The function checks the data and the requested options.
  if (!inherits(data, "usdt_data")) {
    .usdt_stop("`data` must come from usdt_data_long() or usdt_data_tasks(), ",
               "not a plain ", class(data)[1L], ".")
  }
  fix_criteria <- match.arg(fix_criteria)
  .check_confidence_level(level)
  if (!is.character(optimizer) || length(optimizer) != 1L ||
      is.na(optimizer) || !nzchar(optimizer)) {
    .usdt_stop("`optimizer` must be one optimizer name.")
  }
  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))))) {
    .usdt_stop("every argument in `...` must have a name.")
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
               "  Use usdt_boot() for a parametric bootstrap. A future ",
               "Bayesian version will provide regularizing priors.")
  }

  structure(list(fit = fit, tests = tests, pars = pars, design = design,
                 data = data, diagnostics = diag, devfun = devfun,
                 call = match.call(), level = level),
            class = "usdt_freq")
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
             "\n  Check the data and model identification. A future Bayesian ",
             "version will provide regularizing priors.")
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
