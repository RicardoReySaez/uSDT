# population.R
# Population values for planning hierarchical SDT studies
# Author: Ricardo Rey-Sáez
# Last modified: 28-10-2026

# Public functions

#' Population values for planning a uSDT study
#'
#' Describes the population a planned study samples from: how direct and
#' indirect sensitivity (\eqn{d'}) vary and covary across subjects, and how
#' the response criteria vary. The arguments take the names of
#' [usdt_priors()], with numbers instead of priors, and the population can be
#' written in either parameterization of [hsdt()].
#'
#' @param dprime Mean \eqn{d'} of each task. Under the latent regression
#'   parameterization (when `intercept` is given), the mean of the direct task
#'   only.
#' @param sd_dprime Between-subject standard deviation of \eqn{d'} in each
#'   task.
#' @param cor_dprime Correlation between direct and indirect \eqn{d'} across
#'   subjects.
#' @param criterion Mean criterion \eqn{c} of each task, on the classical scale
#'   where positive values are conservative.
#' @param sd_criterion Between-subject standard deviation of the criterion in
#'   each task. Zero gives every subject the same criterion.
#' @param cor_criterion Correlation between the two criteria across subjects.
#' @param sd_signal Standard deviation of the signal distribution of each task,
#'   relative to the noise distribution. The default, 1, is the equal-variance
#'   model that [hsdt()] fits.
#' @param intercept Intercept of the latent regression of indirect on direct
#'   \eqn{d'}. Giving it selects the latent regression parameterization, in
#'   which it replaces the indirect mean.
#' @param from A frequentist [hsdt()] fit, such as a pilot study, whose
#'   estimates give every value not set in the call.
#'
#' @details
#' # Values for each task
#'
#' `dprime`, `sd_dprime`, `criterion`, `sd_criterion` and `sd_signal` take one
#' value for both tasks or one for each, as `c(direct = , indirect = )` or
#' `list(direct = , indirect = )`. The correlations and the intercept take a
#' single value.
#'
#' # The two parameterizations
#'
#' The correlation parameterization gives the bivariate normal distribution of
#' \eqn{d'}: its means \eqn{\mu_D} and \eqn{\mu_I}, standard deviations
#' \eqn{\sigma_D} and \eqn{\sigma_I} and correlation \eqn{\rho}. The latent
#' regression parameterization gives the intercept \eqn{\beta_0} instead of
#' \eqn{\mu_I}, which then follows as
#' \eqn{\mu_I = \beta_0 + \rho \, (\sigma_I / \sigma_D) \, \mu_D}, as in
#' `hsdt(parameterization = "regression")`. Both describe the same
#' distribution, and the printed population reports the tested quantities of
#' either: the difference \eqn{\mu_I - \mu_D} (H1), the correlation (H2) and
#' the intercept and slope of the latent regression (H3).
#'
#' The parameterization decides what stays fixed when a value changes, in a
#' scenario or in [update()]. Under the correlation parameterization a new
#' correlation keeps both means and moves the intercept (unless
#' \eqn{\mu_D = 0}); under the latent regression it keeps the intercept and
#' moves the indirect mean, and with it the difference.
#'
#' # Criteria
#'
#' The criteria matter only through the precision with which each subject's
#' trials measure \eqn{d'}. A task split at the median has its criterion fixed
#' at zero by the split, whatever its value here. The criterion correlation
#' applies when both criteria are estimated.
#'
#' # Scenarios
#'
#' Several values in an argument give one scenario for each, and several
#' arguments with several values give one scenario for each combination. A
#' value shared by both tasks keeps them equal in every scenario:
#' `sd_dprime = c(0.2, 0.3)` gives two scenarios, while
#' `sd_dprime = list(direct = c(0.2, 0.3), indirect = c(0.2, 0.3))` gives four.
#'
#' # Values from a fitted model
#'
#' With `from`, the estimates of a frequentist [hsdt()] fit give the means,
#' standard deviations and correlations, and the arguments set in the call
#' replace them, so `usdt_population(from = m, cor_dprime = 0.3)` keeps every
#' estimate of `m` except the correlation. Setting `intercept` switches to the
#' latent regression parameterization and drops the estimated indirect mean. A
#' criterion fixed at zero by the median split is zero with no variation.
#'
#' @return An object of class `usdt_population`: a data frame with one row per
#'   scenario and the columns `scenario`, `mu_D`, `mu_I`, `sigma_D`,
#'   `sigma_I`, `rho`, the tested quantities `diff`, `intercept` and `slope`,
#'   and `criterion_D`, `criterion_I`, `sd_criterion_D`, `sd_criterion_I`,
#'   `cor_criterion`, `sd_signal_D` and `sd_signal_I`. Its `parameterization`
#'   attribute is `"correlation"` or `"regression"`.
#'
#' @seealso [usdt_priors()], [hsdt()]
#'
#' @examples
#' # The correlation parameterization: means, SDs and correlation of d'
#' pop <- usdt_population(dprime    = c(direct = 0.15, indirect = 0.30),
#'                        sd_dprime = c(direct = 0.35, indirect = 0.30),
#'                        cor_dprime = 0.3)
#' pop
#'
#' # The latent regression parameterization: the intercept replaces the
#' # indirect mean
#' usdt_population(dprime = 0.15, intercept = 0.25,
#'                 sd_dprime = c(direct = 0.35, indirect = 0.30),
#'                 cor_dprime = 0.3)
#'
#' # Three scenarios for the correlation
#' update(pop, cor_dprime = c(0.1, 0.3, 0.5))
#'
#' # The estimates of a pilot study, with a weaker correlation
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
#' usdt_population(from = hsdt(d), cor_dprime = 0.3)
#'
#' @export
usdt_population <- function(dprime        = NULL,
                            sd_dprime     = NULL,
                            cor_dprime    = NULL,
                            criterion     = 0,
                            sd_criterion  = 0.3,
                            cor_criterion = 0,
                            sd_signal     = 1,
                            intercept     = NULL,
                            from          = NULL) {

  # The arguments set in the call replace the defaults or the fitted values.
  set <- setdiff(names(match.call())[-1L], "from")
  if (is.null(from)) {
    values <- list(criterion = criterion, sd_criterion = sd_criterion,
                   cor_criterion = cor_criterion, sd_signal = sd_signal)
    labels <- c(direct = "Direct", indirect = "Indirect")
  } else {
    values <- .population_from_fit(from)
    labels <- from$data$meta$labels

    # An estimate on the boundary cannot describe a population, so the call
    # has to replace it.
    boundary <- setdiff(attr(values, "boundary"), set)
    if (length(boundary)) {
      .usdt_stop("the fit in `from` has ", paste0("`", boundary, "`",
                 collapse = " and "), " on the boundary (a standard deviation ",
                 "of zero or a correlation of -1 or 1), which no population ",
                 "can have. Set ", if (length(boundary) > 1L) "them" else "it",
                 " in the call.")
    }
    attr(values, "boundary") <- NULL
  }
  for (a in set) values[a] <- list(get(a))

  # A fitted indirect mean gives way to an intercept set in the call.
  if (!is.null(from) && !is.null(values$intercept) && !"dprime" %in% set) {
    values$dprime <- .direct_value(values$dprime)
  }
  .population_build(values, labels)
}

#' @param x A `usdt_population` object.
#' @param ... Ignored.
#' @rdname usdt_population
#' @export
print.usdt_population <- function(x, ...) {
  regression <- identical(attr(x, "parameterization"), "regression")
  labels <- attr(x, "labels")
  chars <- .usdt_chars()
  n <- nrow(x)
  cat(.rule("Population", paste0(if (regression) "latent regression" else
    "correlation", if (n > 1L) paste0(" ", chars$dot, " ", n, " scenarios"))),
    "\n\n")

  # Each parameter prints its value, or `varies` when the scenarios differ.
  number <- function(v) formatC(v, format = "f", digits = 3L)
  cell <- function(column) {
    v <- x[[column]]
    if (length(unique(v)) > 1L) "varies" else number(v[1L])
  }
  rows <- list(
    c("mean d'", cell("mu_D"), if (regression) paste0(cell("mu_I"), "*") else
      cell("mu_I")),
    c("sd(d')", cell("sigma_D"), cell("sigma_I")),
    c("cor(d')", cell("rho"), ""),
    if (regression) c("intercept", cell("intercept"), ""),
    c("mean c", cell("criterion_D"), cell("criterion_I")),
    c("sd(c)", cell("sd_criterion_D"), cell("sd_criterion_I")),
    if (any(x$cor_criterion != 0)) c("cor(c)", cell("cor_criterion"), ""),
    if (any(c(x$sd_signal_D, x$sd_signal_I) != 1))
      c("sd(signal)", cell("sd_signal_D"), cell("sd_signal_I")))
  rows <- rbind(c("Parameter", labels[["direct"]], labels[["indirect"]]),
                do.call(rbind, rows))

  # formatC() pads by characters, so labels such as the delta stay aligned.
  for (i in seq_len(nrow(rows))) {
    cat("  ", formatC(rows[i, 1L], width = -12L), " ",
        formatC(rows[i, 2L], width = 10L), " ",
        formatC(rows[i, 3L], width = 10L), "\n", sep = "")
  }

  # The tested quantities follow, for one scenario or for each.
  delta <- paste0(chars$delta, "d'")
  if (n == 1L) {
    tested <- list(c("", "Tested quantity", "Value"),
                   c("H1", paste0(delta, " (indirect - direct)"),
                     number(x$diff)),
                   c("H2", "correlation", number(x$rho)),
                   c("H3", "intercept", number(x$intercept)),
                   c("", "slope", number(x$slope)))
    cat("\n")
    for (t in tested) {
      cat("  ", formatC(t[1L], width = -4L), " ", formatC(t[2L], width = -28L),
          " ", formatC(t[3L], width = 8L), "\n", sep = "")
    }
  } else {
    .print_population_scenarios(x, regression, delta, number)
  }

  cat("\n")
  if (regression) {
    cat("  * The indirect mean follows from the intercept:",
        "mu_I = b0 + rho (sd_I / sd_D) mu_D.\n")
  }
  cat("  A task split at the median has its criterion fixed at 0.\n")
  invisible(x)
}

#' @param object A `usdt_population` object.
#' @param ... For `update()`, the arguments of `usdt_population()` to change.
#'   A task argument may name one task only, which keeps the value of the
#'   other.
#' @rdname usdt_population
#' @export
update.usdt_population <- function(object, ...) {
  changes <- list(...)
  if (length(changes) && (is.null(names(changes)) ||
                          any(!nzchar(names(changes))))) {
    .usdt_stop("every argument of `update()` must have a name.")
  }
  allowed <- setdiff(names(formals(usdt_population)), "from")
  unknown <- setdiff(names(changes), allowed)
  if (length(unknown)) {
    .usdt_stop("`", unknown[1L], "` is not an argument of usdt_population(). ",
               "It takes ", paste0("`", allowed, "`", collapse = ", "), ".")
  }
  values <- attr(object, "values")

  # A task argument that names one task keeps the stored value of the other.
  for (a in names(changes)) {
    new <- changes[[a]]
    task_names <- names(new)
    if (a %in% .population_task_args && length(task_names) == 1L &&
        task_names %in% c("direct", "indirect") &&
        !(a == "dprime" && !is.null(values$intercept))) {
      old <- values[[a]]
      if (!is.list(old)) {
        old <- if (is.null(names(old))) list(direct = old, indirect = old) else
          as.list(old)
      }
      old[[task_names]] <- new[[task_names]]
      values[[a]] <- old
    } else {
      values[a] <- list(new)
    }
  }

  # An intercept added to a correlation population replaces the indirect mean.
  if (!is.null(values$intercept) && is.null(attr(object, "values")$intercept) &&
      !"dprime" %in% names(changes)) {
    values$dprime <- .direct_value(values$dprime)
  }
  .population_build(values, attr(object, "labels"))
}

# A subset of the scenarios is a plain data frame: the printed population and
# update() need every scenario the arguments describe.
#' @export
`[.usdt_population` <- function(x, ...) {
  out <- NextMethod()
  if (is.data.frame(out)) {
    class(out) <- setdiff(class(out), "usdt_population")
    attr(out, "parameterization") <- attr(out, "labels") <-
      attr(out, "values") <- NULL
  }
  out
}

# Internal functions

# The arguments of usdt_population() that take one value for each task.
.population_task_args <- c("dprime", "sd_dprime", "criterion", "sd_criterion",
                           "sd_signal")

# This function keeps the direct value of a task argument: the `direct`
# element when the tasks are named, and the shared value otherwise.
.direct_value <- function(x) {
  if (is.null(names(x))) return(x)
  x[["direct"]]
}

# This function reads the values of a population from a frequentist hsdt()
# fit, on the scales of usdt_population(). The model writes each criterion as
# the intercept of its task, the negative of the classical criterion; under
# treatment coding that intercept is measured from the noise distribution, so
# the classical criterion also subtracts half of d', whose variance it adds.
.population_from_fit <- function(fit) {
  if (!inherits(fit, "hsdt")) {
    .usdt_stop("`from` must come from hsdt(), not a plain ", class(fit)[1L],
               ".")
  }
  if (.is_bayes(fit)) {
    .usdt_stop("`from` takes a frequentist hsdt() fit.")
  }
  e <- fit$pars$est
  sd_d <- sqrt(c(direct = e[["s2_D"]], indirect = e[["s2_I"]]))
  mu_d <- c(direct = e[["gamma_D"]], indirect = e[["gamma_I"]])
  beta <- lme4::fixef(fit$fit)
  re <- fit$pars$re
  treatment <- identical(fit$data$meta$coding, "treatment")

  criterion <- sd_criterion <- c(direct = 0, indirect = 0)
  for (k in c("direct", "indirect")) {
    term <- if (k == "direct") "c_D" else "c_I"
    if (!term %in% fit$design$criteria) next
    sd_c <- re$sd[re$term == term]
    criterion[[k]] <- -beta[[term]] - if (treatment) mu_d[[k]] / 2 else 0
    sd_criterion[[k]] <- sqrt(sd_c^2 + if (treatment) sd_d[[k]]^2 / 4 else 0)
  }
  cor_criterion <- if (is.null(fit$pars$criterion_cor)) 0 else
    fit$pars$criterion_cor$estimate

  # The values a population cannot take are named, so the call can set them.
  boundary <- c(
    if (any(sd_d < 1e-4)) "sd_dprime",
    if (!is.finite(fit$pars$rho) || abs(fit$pars$rho) > 1 - 1e-6) "cor_dprime",
    if (!is.finite(cor_criterion) || abs(cor_criterion) > 1 - 1e-6)
      "cor_criterion")

  structure(list(dprime = as.list(mu_d), sd_dprime = as.list(sd_d),
                 cor_dprime = fit$pars$rho, criterion = as.list(criterion),
                 sd_criterion = as.list(sd_criterion),
                 cor_criterion = cor_criterion, sd_signal = 1),
            boundary = boundary)
}

# This function checks the values of a population, expands its scenarios and
# derives the tested quantities. `values` holds the arguments of
# usdt_population(); the parameterization is the latent regression when it
# holds an intercept.
.population_build <- function(values, labels) {
  for (a in c("dprime", "sd_dprime", "cor_dprime")) {
    if (is.null(values[[a]])) {
      .usdt_stop("`", a, "` is required, unless `from` gives it.")
    }
  }
  regression <- !is.null(values$intercept)

  # The direct mean stands alone under the latent regression.
  if (regression) {
    d <- values$dprime
    task_names <- names(d)
    if ("indirect" %in% task_names) {
      .usdt_stop("`intercept` replaces the indirect mean d', so `dprime` ",
                 "takes the direct mean only.")
    }
    if (!is.null(task_names) && !identical(task_names, "direct")) {
      .usdt_stop("`dprime` names `", task_names[1L], "`; under the latent ",
                 "regression it takes the direct mean only.")
    }
    values$dprime <- unname(if (is.list(d)) d[[1L]] else d)
  } else if (identical(names(values$dprime), "direct")) {
    .usdt_stop("`dprime` needs the indirect mean too, unless `intercept` ",
               "gives the latent regression.")
  }

  # Each argument becomes a table of its scenarios.
  frames <- list(
    dprime = if (regression) .population_values(values$dprime, "dprime",
                                                kind = "mean", tasks = FALSE)
             else .population_values(values$dprime, "dprime", "mean"),
    sd_dprime = .population_values(values$sd_dprime, "sd_dprime", "sd"),
    cor_dprime = .population_values(values$cor_dprime, "cor_dprime", "cor",
                                    tasks = FALSE),
    intercept = if (regression) .population_values(values$intercept,
                                                   "intercept", "mean",
                                                   tasks = FALSE),
    criterion = .population_values(values$criterion, "criterion", "mean"),
    sd_criterion = .population_values(values$sd_criterion, "sd_criterion",
                                      "sd0"),
    cor_criterion = .population_values(values$cor_criterion, "cor_criterion",
                                       "cor", tasks = FALSE),
    sd_signal = .population_values(values$sd_signal, "sd_signal", "sd"))
  frames <- frames[!vapply(frames, is.null, TRUE)]

  # Every combination of the scenarios of each argument is one scenario.
  grid <- expand.grid(lapply(frames, function(f) seq_len(nrow(f))),
                      KEEP.OUT.ATTRS = FALSE)
  pick <- function(arg, column) frames[[arg]][[column]][grid[[arg]]]

  mu_D <- if (regression) pick("dprime", "value") else pick("dprime", "direct")
  sigma_D <- pick("sd_dprime", "direct")
  sigma_I <- pick("sd_dprime", "indirect")
  rho <- pick("cor_dprime", "value")
  mu_I <- if (regression) {
    pick("intercept", "value") + rho * sigma_I / sigma_D * mu_D
  } else {
    pick("dprime", "indirect")
  }
  B <- cbind(mu_D = mu_D, mu_I = mu_I, sigma_D = sigma_D, sigma_I = sigma_I,
             rho = rho)
  q <- .usdt_quantities(.bivariate_primitives(B))

  out <- data.frame(
    scenario = seq_len(nrow(grid)), B,
    diff = unname(q[, "diff"]),
    intercept = if (regression) pick("intercept", "value") else
      unname(q[, "intercept"]),
    slope = unname(q[, "slope"]),
    criterion_D = pick("criterion", "direct"),
    criterion_I = pick("criterion", "indirect"),
    sd_criterion_D = pick("sd_criterion", "direct"),
    sd_criterion_I = pick("sd_criterion", "indirect"),
    cor_criterion = pick("cor_criterion", "value"),
    sd_signal_D = pick("sd_signal", "direct"),
    sd_signal_I = pick("sd_signal", "indirect"),
    row.names = NULL)

  if (!regression) values$intercept <- NULL
  if (regression) values$dprime <- list(direct = values$dprime)
  structure(out, class = c("usdt_population", "data.frame"),
            parameterization = if (regression) "regression" else "correlation",
            labels = labels, values = values)
}

# This function checks one argument of usdt_population() and returns its
# scenarios, one row each: a `value` column for an argument without tasks,
# and `direct` and `indirect` columns for one with tasks. A value shared by
# both tasks keeps them equal; values given for each task are crossed. `kind`
# is the range of the argument: any real value (`mean`), positive (`sd`),
# non-negative (`sd0`) or a correlation (`cor`).
.population_values <- function(x, arg, kind, tasks = TRUE) {
  check <- function(v, what = arg) {
    if (is.list(v)) v <- unlist(v, use.names = FALSE)
    if (!is.numeric(v) || !length(v) || anyNA(v) || any(!is.finite(v))) {
      .usdt_stop("`", what, "` must be one or more finite numbers.")
    }
    bad <- switch(kind, mean = FALSE, sd = any(v <= 0), sd0 = any(v < 0),
                  cor = any(abs(v) >= 1))
    if (bad) {
      .usdt_stop("`", what, "` ", switch(kind,
        sd = "is a standard deviation, so it must be positive.",
        sd0 = "is a standard deviation, so it cannot be negative.",
        cor = "is a correlation, so it must lie inside (-1, 1)."))
    }
    unname(as.numeric(v))
  }

  if (!tasks) {
    if (!is.null(names(x))) {
      .usdt_stop("`", arg, "` takes a single value, not one per task.")
    }
    return(data.frame(value = unique(check(x))))
  }
  task_names <- names(x)
  if (is.null(task_names)) {
    v <- unique(check(x))
    return(data.frame(direct = v, indirect = v))
  }
  per_task <- .per_task(x, arg)
  expand.grid(direct = unique(check(per_task$direct,
                                    paste0(arg, "$direct"))),
              indirect = unique(check(per_task$indirect,
                                      paste0(arg, "$indirect"))),
              KEEP.OUT.ATTRS = FALSE)
}

# This function prints one row per scenario with the parameters that vary and
# the tested quantities.
.print_population_scenarios <- function(x, regression, delta, number) {
  inputs <- c(mu_D = "mean d' D", mu_I = "mean d' I", sigma_D = "sd(d') D",
              sigma_I = "sd(d') I", rho = "cor(d')", intercept = "intercept",
              criterion_D = "mean c D", criterion_I = "mean c I",
              sd_criterion_D = "sd(c) D", sd_criterion_I = "sd(c) I",
              cor_criterion = "cor(c)", sd_signal_D = "sd(sig) D",
              sd_signal_I = "sd(sig) I")
  if (!regression) inputs <- inputs[names(inputs) != "intercept"]
  varies <- names(inputs)[vapply(names(inputs), function(column) {
    length(unique(x[[column]])) > 1L
  }, TRUE)]
  derived <- c(diff = delta, intercept = "intercept", slope = "slope")
  derived <- derived[setdiff(names(derived), varies)]
  columns <- c(inputs[varies], derived)

  cat("\n  ", formatC("Scenario", width = -9L),
      paste(formatC(columns, width = 10L), collapse = ""), "\n", sep = "")
  for (i in seq_len(nrow(x))) {
    cat("  ", formatC(x$scenario[i], width = -9L),
        paste(formatC(number(unlist(x[i, names(columns)])), width = 10L),
              collapse = ""), "\n", sep = "")
  }
  invisible(NULL)
}
