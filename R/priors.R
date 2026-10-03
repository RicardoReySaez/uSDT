# priors.R
# Specify the priors of the Bayesian hierarchical SDT model
# Author: Ricardo Rey-Sáez
# Last modified: 03-10-2026

# Public functions

#' Priors for the Bayesian hierarchical SDT model
#'
#' Specifies the prior distributions of the Bayesian model. The family of every
#' prior is fixed by the model; you choose its parameters. Each prior is
#' written as a string in Stan notation:
#' * means of \eqn{d'} and of the criterion: `"normal(location, scale)"`;
#' * between-subject standard deviations: `"student_t(df, 0, scale)"`, a
#'   Student-\eqn{t} distribution truncated at zero (half-Student-\eqn{t});
#' * correlations: `"scaled_beta(alpha, beta)"`, a beta distribution stretched
#'   from \eqn{(0, 1)} to \eqn{(-1, 1)};
#' * signal standard deviation under unequal variances:
#'   `"lognormal(location, scale)"`.
#'
#' @param dprime Prior of the mean \eqn{d'} of each task.
#' @param sd_dprime Prior of the between-subject standard deviation of
#'   \eqn{d'} in each task.
#' @param cor_dprime Prior of the correlation between direct and indirect
#'   \eqn{d'} across subjects.
#' @param criterion Prior of the mean criterion \eqn{c} of each task, on the
#'   classical scale where positive values are conservative.
#' @param sd_criterion Prior of the between-subject standard deviation of the
#'   criterion in each task.
#' @param cor_criterion Prior of the correlation between the two criteria. It is
#'   used only when both criteria are estimated.
#' @param sd_signal Prior of the standard deviation of the signal distribution,
#'   relative to the noise distribution. It is used only under unequal
#'   variances.
#'
#' @details
#' `dprime`, `sd_dprime`, `criterion`, `sd_criterion` and `sd_signal` take one
#' prior for both tasks or one per task as `list(direct = ..., indirect = ...)`.
#' The correlations take a single prior.
#'
#' The defaults are centred at zero and weakly informative on the probit
#' scale: `normal(0, 1)` for every mean, `student_t(4, 0, 1)` for every
#' between-subject standard deviation and a uniform `scaled_beta(1, 1)` for
#' every correlation. With equal locations for both means, the prior of the
#' difference between them (H1) is centred at zero. A symmetric scaled beta,
#' `scaled_beta(a, a)`, equals the marginal of an LKJ(a) prior on a 2 x 2
#' correlation matrix and is centred at zero (H2).
#'
#' The slope and intercept of the latent regression (H3) have no priors of their
#' own: their priors follow from those of the means, standard deviations and
#' correlation. With the defaults, the intercept prior is symmetric around zero
#' but has heavy tails. It does not change when both standard deviation priors
#' are rescaled by the same factor, because the slope depends on their ratio.
#'
#' @return An object of class `usdt_priors`: a data frame with one row per
#'   prior, giving the parameter, the task, the family and its arguments.
#'
#' @examples
#' # Default priors
#' usdt_priors()
#'
#' # A narrower prior for the mean indirect sensitivity
#' usdt_priors(dprime = list(direct = "normal(0, 1)",
#'                           indirect = "normal(0, 0.5)"))
#'
#' # A correlation prior that concentrates near zero
#' usdt_priors(cor_dprime = "scaled_beta(2, 2)")
#'
#' @export
usdt_priors <- function(dprime        = "normal(0, 1)",
                        sd_dprime     = "student_t(4, 0, 1)",
                        cor_dprime    = "scaled_beta(1, 1)",
                        criterion     = "normal(0, 1)",
                        sd_criterion  = "student_t(4, 0, 1)",
                        cor_criterion = "scaled_beta(1, 1)",
                        sd_signal     = "lognormal(0, 0.5)") {

  # Each argument has one admissible family; task arguments may vary by task.
  spec <- list(
    dprime        = list(value = dprime,        family = "normal",      tasks = TRUE),
    sd_dprime     = list(value = sd_dprime,     family = "student_t",   tasks = TRUE),
    cor_dprime    = list(value = cor_dprime,    family = "scaled_beta", tasks = FALSE),
    criterion     = list(value = criterion,     family = "normal",      tasks = TRUE),
    sd_criterion  = list(value = sd_criterion,  family = "student_t",   tasks = TRUE),
    cor_criterion = list(value = cor_criterion, family = "scaled_beta", tasks = FALSE),
    sd_signal     = list(value = sd_signal,     family = "lognormal",   tasks = TRUE)
  )

  rows <- lapply(names(spec), function(parameter) {
    s <- spec[[parameter]]
    values <- if (s$tasks) .per_task(s$value, parameter) else list(both = s$value)
    do.call(rbind, lapply(names(values), function(task) {
      cbind(data.frame(parameter = parameter, task = task,
                       stringsAsFactors = FALSE),
            .parse_prior(values[[task]], parameter, s$family))
    }))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  structure(out, class = c("usdt_priors", "data.frame"))
}

#' @param x A `usdt_priors` object.
#' @param ... Ignored.
#' @rdname usdt_priors
#' @export
print.usdt_priors <- function(x, ...) {
  labels <- c(dprime = "mean d'", sd_dprime = "sd(d')", cor_dprime = "cor(d')",
              criterion = "mean c", sd_criterion = "sd(c)",
              cor_criterion = "cor(c)", sd_signal = "sd(signal)")
  cat(.rule("Priors"), "\n\n")
  cat(sprintf("  %-12s %-10s %s\n", "Parameter", "Task", "Prior"))
  for (i in seq_len(nrow(x))) {
    cat(sprintf("  %-12s %-10s %s\n", labels[[x$parameter[i]]], x$task[i],
                x$prior[i]))
  }
  cat("\n  Standard deviation priors are truncated at zero. cor(c) applies when\n",
      "  both criteria are estimated, and sd(signal) under unequal variances.\n",
      sep = "")
  invisible(x)
}

# Internal functions

# The arguments of each admissible family, in Stan order.
.prior_families <- list(
  normal      = c("location", "scale"),
  student_t   = c("df", "location", "scale"),
  scaled_beta = c("alpha", "beta"),
  lognormal   = c("location", "scale")
)

# This function reads one prior string and checks its arguments.
.parse_prior <- function(x, arg, family) {
  example <- switch(family, normal = "normal(0, 1)",
                    student_t = "student_t(4, 0, 1)",
                    scaled_beta = "scaled_beta(1, 1)",
                    lognormal = "lognormal(0, 0.5)")
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    .usdt_stop("`", arg, "` must be one prior written as a string, e.g. \"",
               example, "\".")
  }
  parts <- regmatches(x, regexec("^\\s*([A-Za-z_]+)\\s*\\((.*)\\)\\s*$", x))[[1L]]
  if (!length(parts) || parts[2L] != family) {
    .usdt_stop("`", arg, "` must be a ", family, " prior, e.g. \"", example,
               "\", not \"", x, "\".")
  }

  arg_names <- .prior_families[[family]]
  values <- suppressWarnings(as.numeric(strsplit(parts[3L], ",", fixed = TRUE)[[1L]]))
  if (length(values) != length(arg_names) || any(!is.finite(values))) {
    .usdt_stop("`", arg, "` needs ", length(arg_names), " numbers (",
               paste(arg_names, collapse = ", "), "), e.g. \"", example, "\".")
  }
  values <- stats::setNames(values, arg_names)
  positive <- intersect(c("df", "scale", "alpha", "beta"), arg_names)
  if (any(values[positive] <= 0)) {
    .usdt_stop("`", arg, "` needs positive ",
               paste(positive, collapse = " and "), ".")
  }
  if (family == "student_t" && values[["location"]] != 0) {
    .usdt_stop("`", arg, "` is a standard deviation prior, so its location ",
               "must be 0: the Student-t is truncated at zero.")
  }

  out <- data.frame(family = family, df = NA_real_, location = NA_real_,
                    scale = NA_real_, alpha = NA_real_, beta = NA_real_,
                    prior = paste0(family, "(", paste(values, collapse = ", "), ")"),
                    stringsAsFactors = FALSE)
  out[arg_names] <- as.list(values)
  out
}

# This function returns the rows of one prior parameter in the given task order.
.prior_rows <- function(priors, parameter, tasks = c("direct", "indirect")) {
  rows <- priors[priors$parameter == parameter, , drop = FALSE]
  rows[match(tasks, rows$task), , drop = FALSE]
}

# This function turns the priors into the hyperparameters of the Stan model.
# Criterion priors are kept only for the estimated criteria, and the signal SD
# priors only under unequal variances.
.stan_priors <- function(priors, free_c, UV) {
  vec <- function(x) array(x, dim = length(x))
  d <- .prior_rows(priors, "dprime")
  sd_d <- .prior_rows(priors, "sd_dprime")
  cor_d <- .prior_rows(priors, "cor_dprime", "both")
  c_ <- .prior_rows(priors, "criterion", c("direct", "indirect")[free_c == 1L])
  sd_c <- .prior_rows(priors, "sd_criterion", c("direct", "indirect")[free_c == 1L])
  cor_c <- .prior_rows(priors, "cor_criterion", "both")
  s <- .prior_rows(priors, "sd_signal", if (UV) c("direct", "indirect") else character(0))
  list(mu_d_loc = vec(d$location), mu_d_scale = vec(d$scale),
       sd_d_df = vec(sd_d$df), sd_d_scale = vec(sd_d$scale),
       rho_d_a = cor_d$alpha, rho_d_b = cor_d$beta,
       mu_c_loc = vec(c_$location), mu_c_scale = vec(c_$scale),
       sd_c_df = vec(sd_c$df), sd_c_scale = vec(sd_c$scale),
       rho_c_a = cor_c$alpha, rho_c_b = cor_c$beta,
       sd_s_loc = vec(s$location), sd_s_scale = vec(s$scale))
}

# This function simulates the five values that describe the two sensitivities
# from their priors. Passed through .usdt_quantities(), the draws give the
# priors that the model induces on the difference, correlation, slope and
# intercept.
.prior_draws <- function(priors, n) {
  d <- .prior_rows(priors, "dprime")
  s <- .prior_rows(priors, "sd_dprime")
  r <- .prior_rows(priors, "cor_dprime", "both")
  mu <- vapply(1:2, function(j) stats::rnorm(n, d$location[j], d$scale[j]),
               numeric(n))
  sigma <- vapply(1:2, function(j) abs(stats::rt(n, s$df[j])) * s$scale[j],
                  numeric(n))
  rho <- 2 * stats::rbeta(n, r$alpha, r$beta) - 1
  cbind(gamma_D = mu[, 1L], gamma_I = mu[, 2L],
        s2_D = sigma[, 1L]^2, s2_I = sigma[, 2L]^2,
        s_DI = rho * sigma[, 1L] * sigma[, 2L])
}
