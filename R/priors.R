# priors.R
# Specify the priors of the Bayesian hierarchical SDT model
# Author: Ricardo Rey-Sáez
# Last modified: 15-10-2026

# Public functions

#' Priors for the Bayesian hierarchical SDT model
#'
#' Specifies the prior distributions of the Bayesian model. Each prior is
#' written as a string in Stan notation, `"family(arguments)"`, with the
#' families and parameterisations of Stan, e.g. `"normal(0, 1)"`,
#' `"student_t(4, 0, 0.5)"` or `"gamma(2, 4)"`. Each argument takes the
#' families that cover the values its parameter can take:
#' * means of \eqn{d'} and of the criterion can take any real value, so their
#'   prior must be defined on the whole real line: `normal`, `std_normal`,
#'   `student_t`, `cauchy`, `double_exponential`, `logistic`, `gumbel`,
#'   `skew_normal`, `exp_mod_normal` or `skew_double_exponential`;
#' * standard deviations take any of those, truncated at zero (with location
#'   0, a `student_t` is a half-Student-\eqn{t}), or a family on the positive
#'   values: `lognormal`, `exponential`, `gamma`, `inv_gamma`, `weibull`,
#'   `frechet`, `rayleigh`, `chi_square`, `inv_chi_square`,
#'   `scaled_inv_chi_square`, `loglogistic`, or `pareto_type_2` with
#'   `mu <= 0`;
#' * correlations take `"lkj(eta)"`, the LKJ prior on a correlation matrix and
#'   the only correlation prior brms accepts. On a 2 x 2 matrix the
#'   correlation follows a beta(eta, eta) distribution stretched from
#'   \eqn{(0, 1)} to \eqn{(-1, 1)}: `lkj(1)` is uniform, and larger values of
#'   eta concentrate it around zero.
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
#' @param intercept,slope Priors of the intercept and slope of the latent
#'   regression of indirect on direct \eqn{d'}. They are used only by the
#'   latent regression parameterization.
#' @param sd_residual Prior of the residual standard deviation of the
#'   indirect \eqn{d'} around the latent regression line. It is used only by
#'   the latent regression parameterization.
#'
#' @details
#' `dprime`, `sd_dprime`, `criterion`, `sd_criterion` and `sd_signal` take one
#' prior for both tasks or one per task as `list(direct = ..., indirect = ...)`.
#' The correlations take a single prior.
#'
#' The defaults are centred at zero and weakly informative on the probit
#' scale: `normal(0, 1)` for the mean \eqn{d'}, `normal(0, 0.5)` for the mean
#' criterion, `student_t(4, 0, 0.5)` for every between-subject standard
#' deviation and a uniform `lkj(1)` for every correlation. With them, the hit
#' and false-alarm rates of a new subject spread over \eqn{(0, 1)} without
#' piling up at either end; see [plot.usdt_priors()]. With the same prior for
#' both means, the prior of the difference between them (H1) is symmetric
#' around zero, and every `lkj` prior is centred at zero (H2).
#'
#' The model describes the two sensitivities by their means, standard
#' deviations and correlation, so the slope and intercept of the latent
#' regression (H3) have no priors of their own: their priors follow from
#' those of the means, standard deviations and correlation. With the defaults,
#' the intercept prior is symmetric around zero but has heavy tails. It does
#' not change when both standard deviation priors are rescaled by the same
#' factor, because the slope depends on their ratio.
#'
#' The latent regression parameterization describes the indirect sensitivity
#' of each subject instead as `intercept + slope * d'(direct) + residual`, with
#' priors on the intercept, the slope and the residual standard deviation, and
#' the indirect mean, standard deviation and correlation follow from them.
#' Then `dprime` and `sd_dprime` apply to the direct task only, and
#' `cor_dprime` is not used.
#'
#' @return An object of class `usdt_priors`: a data frame with one row per
#'   prior, giving the parameter, the task, the family, its arguments (a list
#'   column) and the prior as a string.
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
#' usdt_priors(cor_dprime = "lkj(2)")
#'
#' # Standard deviations of d' expected around 0.3 rather than near zero
#' usdt_priors(sd_dprime = "student_t(4, 0.3, 0.2)")
#'
#' # Any family that covers the positive values, e.g. a gamma distribution
#' usdt_priors(sd_dprime = "gamma(2, 4)")
#'
#' @export
usdt_priors <- function(dprime        = "normal(0, 1)",
                        sd_dprime     = "student_t(4, 0, 0.5)",
                        cor_dprime    = "lkj(1)",
                        criterion     = "normal(0, 0.5)",
                        sd_criterion  = "student_t(4, 0, 0.5)",
                        cor_criterion = "lkj(1)",
                        sd_signal     = "lognormal(0, 0.5)",
                        intercept     = "normal(0, 0.5)",
                        slope         = "normal(0, 1)",
                        sd_residual   = "student_t(4, 0, 0.5)") {
  values <- as.list(environment())

  # Task arguments take one prior for both tasks or one per task.
  rows <- lapply(names(.prior_parameters), function(parameter) {
    value <- values[[parameter]]
    tasks <- if (.prior_parameters[[parameter]]$tasks) {
      .per_task(value, parameter)
    } else {
      list(both = value)
    }
    lapply(names(tasks), function(task) {
      c(list(parameter = parameter, task = task),
        .parse_prior(tasks[[task]], parameter))
    })
  })
  .priors_frame(unlist(rows, recursive = FALSE))
}

#' @param x A `usdt_priors` object.
#' @param ... Ignored.
#' @rdname usdt_priors
#' @export
print.usdt_priors <- function(x, ...) {
  cat(.rule("Priors"), "\n\n")
  .print_prior_rows(x)
  cat("\n  Standard deviation priors are truncated at zero. cor(c) applies when\n",
      "  both criteria are estimated, and sd(signal) under unequal variances.\n",
      "  Under the latent regression, intercept, slope and sd(residual) replace\n",
      "  the indirect mean d', sd(d') and cor(d').\n",
      sep = "")
  invisible(x)
}

# Internal functions

# The parameters that take a prior: whether each task has its own, the kind
# of value they are, and the label they print with.
.prior_parameters <- list(
  dprime        = list(tasks = TRUE,  kind = "mean", label = "mean d'"),
  sd_dprime     = list(tasks = TRUE,  kind = "sd",   label = "sd(d')"),
  cor_dprime    = list(tasks = FALSE, kind = "cor",  label = "cor(d')"),
  criterion     = list(tasks = TRUE,  kind = "mean", label = "mean c"),
  sd_criterion  = list(tasks = TRUE,  kind = "sd",   label = "sd(c)"),
  cor_criterion = list(tasks = FALSE, kind = "cor",  label = "cor(c)"),
  sd_signal     = list(tasks = TRUE,  kind = "sd",   label = "sd(signal)"),
  intercept     = list(tasks = FALSE, kind = "mean", label = "intercept"),
  slope         = list(tasks = FALSE, kind = "mean", label = "slope"),
  sd_residual   = list(tasks = FALSE, kind = "sd",   label = "sd(residual)")
)

# The values each kind of parameter can take.
.prior_domains <- list(mean = c(-Inf, Inf), sd = c(0, Inf), cor = c(-1, 1))

# This function prints one line per prior.
.print_prior_rows <- function(x) {
  cat(sprintf("  %-12s %-10s %s\n", "Parameter", "Task", "Prior"))
  for (i in seq_len(nrow(x))) {
    cat(sprintf("  %-12s %-10s %s\n", .prior_parameters[[x$parameter[i]]]$label,
                x$task[i], x$prior[i]))
  }
  invisible(NULL)
}

# This function keeps the priors a fit with the given design uses: those of
# its parameterization of the two sensitivities, of the estimated criteria, of
# the criterion correlation when both are estimated, and of the signal SDs
# under unequal variances.
.used_priors <- function(priors, design) {
  free <- c("direct", "indirect")[design$free_c == 1L]
  sensitivities <- if (design$parameterization == "regression") {
    (priors$parameter %in% c("dprime", "sd_dprime") &
       priors$task == "direct") |
      priors$parameter %in% c("intercept", "slope", "sd_residual")
  } else {
    priors$parameter %in% c("dprime", "sd_dprime", "cor_dprime")
  }
  keep <- sensitivities |
    (priors$parameter %in% c("criterion", "sd_criterion") &
       priors$task %in% free) |
    (priors$parameter == "cor_criterion" & length(free) == 2L) |
    (priors$parameter == "sd_signal" & design$unequal_variances)
  priors[keep, , drop = FALSE]
}

# This function reads one prior of usdt_priors() and checks that its family is
# known, that its arguments are valid, and that it covers every value its
# parameter can take.
.parse_prior <- function(x, parameter) {
  kind <- .prior_parameters[[parameter]]$kind
  example <- formals(usdt_priors)[[parameter]]
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    .usdt_stop("`", parameter, "` must be one prior written as a string, ",
               "e.g. \"", example, "\".")
  }
  read <- .read_prior(x)
  if (is.null(read)) {
    .usdt_stop("cannot read the prior \"", x, "\" of `", parameter, "`. ",
               "Write it as family(arguments), e.g. \"", example, "\".")
  }
  family <- .distributions[[read$family]]
  if (is.null(family)) {
    .usdt_stop("`", parameter, "` uses `", read$family, "`, which is not a ",
               "Stan distribution uSDT can use as a prior. The families that ",
               "suit `", parameter, "` are ",
               paste0("`", .families_for(kind), "`", collapse = ", "), ".")
  }
  n <- length(family$args)
  if (length(read$args) != n || !all(is.finite(read$args))) {
    .usdt_stop("`", parameter, "`: ", read$family, " takes ",
               if (n) paste0(n, " number", if (n > 1L) "s", " (",
                             paste(family$args, collapse = ", "), ")")
               else "no arguments",
               ", not \"", x, "\".")
  }
  if (!family$valid(read$args)) {
    .usdt_stop("`", parameter, "`: ", read$family, " needs ", family$rule, ".")
  }

  prior <- paste0(read$family, "(", paste(read$args, collapse = ", "), ")")
  if (kind == "cor" && read$family != "lkj") {
    .usdt_stop("`", parameter, "` is a correlation, so its prior must be ",
               "lkj(eta), the only correlation prior brms accepts, e.g. ",
               "\"lkj(1)\".")
  }
  support <- family$support(read$args)
  domain <- .prior_domains[[kind]]
  if (support[1L] > domain[1L] || support[2L] < domain[2L]) {
    .usdt_stop("`", parameter, "` ",
               switch(kind, mean = "can take any real value",
                      sd = "is a standard deviation, which can take any positive value"),
               ", but ", prior, " is defined only on [", support[1L], ", ",
               support[2L], "].")
  }
  list(family = read$family, args = read$args, prior = prior)
}

# This function names the families whose support covers every value of a kind
# of parameter whatever their arguments.
.families_for <- function(kind) {
  domain <- .prior_domains[[kind]]
  covers <- vapply(.distributions, function(f) {
    support <- f$support(rep(NA_real_, length(f$args)))
    isTRUE(support[1L] <= domain[1L] && support[2L] >= domain[2L])
  }, TRUE)
  names(.distributions)[covers]
}

# This function builds a usdt_priors table from a list of rows, each a list
# with the parameter, the task, the family, its arguments and the prior text.
.priors_frame <- function(rows) {
  field <- function(name) vapply(rows, `[[`, "", name)
  out <- data.frame(parameter = field("parameter"), task = field("task"),
                    family = field("family"), stringsAsFactors = FALSE)
  out$args <- lapply(rows, `[[`, "args")
  out$prior <- field("prior")
  structure(out, class = c("usdt_priors", "data.frame"))
}

# This function returns the rows of one prior parameter in the given task order.
.prior_rows <- function(priors, parameter, tasks = c("direct", "indirect")) {
  rows <- priors[priors$parameter == parameter, , drop = FALSE]
  rows[match(tasks, rows$task), , drop = FALSE]
}

# This function turns the priors into the hyperparameters of the Stan model,
# which has the correlation parameterization. The model takes normal means,
# Student-t SDs, LKJ correlations and lognormal signal SDs; criterion priors
# are kept only for the estimated criteria, and the signal SD priors only
# under unequal variances.
.stan_priors <- function(priors, design) {
  used <- .used_priors(priors, design)
  free_c <- design$free_c
  UV <- design$unequal_variances
  families <- c(dprime = "normal", sd_dprime = "student_t", cor_dprime = "lkj",
                criterion = "normal", sd_criterion = "student_t",
                cor_criterion = "lkj", sd_signal = "lognormal")
  wrong <- which(used$family != families[used$parameter])
  if (length(wrong)) {
    i <- wrong[1L]
    .usdt_stop("the Stan model takes ", families[[used$parameter[i]]],
               " priors for `", used$parameter[i], "`, not ", used$prior[i],
               ".")
  }

  # Stan reads every hyperparameter of a task as a vector, even of length 1.
  free <- c("direct", "indirect")[free_c == 1L]
  signal <- if (UV) c("direct", "indirect") else character(0)
  arg <- function(parameter, k, tasks = c("direct", "indirect")) {
    x <- vapply(.prior_rows(priors, parameter, tasks)$args, `[`, 0, k)
    array(x, dim = length(x))
  }
  eta <- function(parameter) .prior_rows(priors, parameter, "both")$args[[1L]]
  list(mu_d_loc = arg("dprime", 1L), mu_d_scale = arg("dprime", 2L),
       sd_d_df = arg("sd_dprime", 1L), sd_d_loc = arg("sd_dprime", 2L),
       sd_d_scale = arg("sd_dprime", 3L),
       rho_d_a = eta("cor_dprime"), rho_d_b = eta("cor_dprime"),
       mu_c_loc = arg("criterion", 1L, free),
       mu_c_scale = arg("criterion", 2L, free),
       sd_c_df = arg("sd_criterion", 1L, free),
       sd_c_loc = arg("sd_criterion", 2L, free),
       sd_c_scale = arg("sd_criterion", 3L, free),
       rho_c_a = eta("cor_criterion"), rho_c_b = eta("cor_criterion"),
       sd_s_loc = arg("sd_signal", 1L, signal),
       sd_s_scale = arg("sd_signal", 2L, signal))
}

# Number of draws that approximate a prior the model induces rather than
# states, such as the prior of the latent intercept, and the smaller number
# that suffices to draw it.
.prior_samples <- 1e5L
.curve_samples <- 2e4L

# This function returns the log density, the distribution function and a
# generator of one prior, restricted to the values its parameter can take: a
# standard deviation prior defined on the whole real line is truncated at zero
# and renormalised. Truncated draws are kept by rejection.
.prior_fun <- function(row) {
  a <- row$args[[1L]]
  spec <- .prior_parameters[[row$parameter]]
  if (!.known_prior(list(family = row$family, args = a))) {
    .usdt_stop("uSDT cannot simulate the ", spec$label, " prior \"",
               row$prior, "\". The families it can use are ",
               paste0("`", names(.distributions), "`", collapse = ", "), ".")
  }
  family <- .distributions[[row$family]]
  support <- family$support(a)
  domain <- .prior_domains[[spec$kind]]
  lower <- max(domain[1L], support[1L])
  upper <- min(domain[2L], support[2L])

  # The distribution function of the family, at values clamped to the domain.
  p <- function(q) {
    q <- pmin(pmax(q, lower), upper)
    inside <- q > support[1L] & q < support[2L]
    out <- as.numeric(q >= support[2L])
    out[inside] <- family$p(q[inside], a)
    out
  }
  base <- p(lower)
  mass <- p(upper) - base
  truncated <- lower > support[1L] || upper < support[2L]

  list(
    log_density = function(x) {
      if (lower == -Inf && upper == Inf) return(family$logd(x, a))
      out <- rep(-Inf, length(x))
      inside <- x >= lower & x <= upper
      out[inside] <- family$logd(x[inside], a) - log(mass)
      out
    },
    cdf = function(q) (p(q) - base) / mass,
    r = function(n) {
      if (!truncated) return(family$r(n, a))
      if (mass < 1e-4) {
        .usdt_stop("the ", spec$label, " prior \"", row$prior, "\" puts ",
                   "almost no mass on the values its parameter can take.")
      }
      draws <- numeric(0)
      while (length(draws) < n) {
        x <- family$r(ceiling((n - length(draws)) * 1.1 / mass) + 10L, a)
        draws <- c(draws, x[x >= lower & x <= upper])
      }
      draws[seq_len(n)]
    })
}

# This function simulates the population distribution of the two sensitivities
# from its priors, one draw per row.
.prior_draws <- function(priors, n) {
  d <- .prior_rows(priors, "dprime")
  s <- .prior_rows(priors, "sd_dprime")
  draw <- function(row) .prior_fun(row)$r(n)
  cbind(mu_D = draw(d[1L, ]), mu_I = draw(d[2L, ]),
        sigma_D = draw(s[1L, ]), sigma_I = draw(s[2L, ]),
        rho = draw(.prior_rows(priors, "cor_dprime", "both")))
}

# This function adds one new subject to every population draw and returns that
# subject's hit and false-alarm rates in each task, following the subject level
# of the model. A criterion fixed by the Meyen split takes the value that makes
# HR + FAR = 1.
.prior_predictive <- function(priors, n, free_c, UV) {
  draw <- function(row) .prior_fun(row)$r(n)
  B <- .prior_draws(priors, n)
  z <- matrix(stats::rnorm(2L * n), n)
  d <- cbind(B[, "mu_D"] + B[, "sigma_D"] * z[, 1L],
             B[, "mu_I"] + B[, "sigma_I"] *
               (B[, "rho"] * z[, 1L] + sqrt(1 - B[, "rho"]^2) * z[, 2L]))
  signal <- .prior_rows(priors, "sd_signal")
  sigma_s <- if (UV) {
    vapply(1:2, function(j) draw(signal[j, ]), numeric(n))
  } else {
    matrix(1, n, 2L)
  }

  criterion <- d * (1 - sigma_s) / (2 * (1 + sigma_s))
  free <- which(free_c == 1L)
  if (length(free)) {
    mean_c <- .prior_rows(priors, "criterion")
    sd_c <- .prior_rows(priors, "sd_criterion")
    w <- matrix(stats::rnorm(2L * n), n)
    if (length(free) == 2L) {
      r <- draw(.prior_rows(priors, "cor_criterion", "both"))
      w[, 2L] <- r * w[, 1L] + sqrt(1 - r^2) * w[, 2L]
    }
    for (k in seq_along(free)) {
      j <- free[k]
      criterion[, j] <- draw(mean_c[j, ]) + draw(sd_c[j, ]) * w[, k]
    }
  }

  list(draws = B,
       hit = stats::pnorm((d / 2 - criterion) / sigma_s),
       fa = stats::pnorm(-d / 2 - criterion))
}

# This function returns the log density and the distribution function of the
# prior of each requested tested quantity. Each one is a + b * theta, where
# theta has a prior of its own and a and b > 0 depend on the other parameters,
# so its prior averages the transformed prior of theta over prior draws of a
# and b (a Rao-Blackwell estimate over n draws), whatever the families of the
# priors. The correlation, with a and b fixed, is exact.
.quantity_priors <- function(priors, quantities, n = .prior_samples) {
  mu <- .prior_rows(priors, "dprime")
  mu_I <- .prior_fun(mu[2L, ])
  rho <- .prior_fun(.prior_rows(priors, "cor_dprime", "both"))
  out <- list(rho = .affine_prior(rho, 0, 1, bounds = c(-1, 1)))
  if (any(c("diff", "intercept") %in% quantities)) {
    mu_D <- .prior_fun(mu[1L, ])$r(n)
  }
  if ("diff" %in% quantities) out$diff <- .affine_prior(mu_I, -mu_D, 1)
  if (any(c("intercept", "slope") %in% quantities)) {
    s <- .prior_rows(priors, "sd_dprime")
    ratio <- .prior_fun(s[2L, ])$r(n) / .prior_fun(s[1L, ])$r(n)
    slope <- rho$r(n) * ratio
    if ("intercept" %in% quantities) {
      out$intercept <- .affine_prior(mu_I, -slope * mu_D, 1)
    }
    out$slope <- .affine_prior(rho, 0, ratio)
  }
  out
}

# This function gives the prior of a + b * theta from that of theta (from
# .prior_fun()) and from draws, or fixed values, of a and b > 0.
.affine_prior <- function(theta, a, b, bounds = NULL) {
  list(
    log_density = function(x) vapply(x, function(v) {
      .log_mean_exp(theta$log_density((v - a) / b) - log(b))
    }, 0),
    cdf = function(x) vapply(x, function(v) mean(theta$cdf((v - a) / b)), 0),
    bounds = bounds)
}

# This function averages probabilities given on the log scale.
.log_mean_exp <- function(l) {
  top <- max(l)
  if (!is.finite(top)) return(top)
  top + log(mean(exp(l - top)))
}

# Prior predictive check

#' Prior predictive check for the Bayesian hierarchical SDT model
#'
#' Draws what the priors imply before seeing any data: the hit and
#' false-alarm rates of a new subject in each task, and the prior of the
#' quantity tested by each hypothesis. Use it to check that the priors allow
#' the data you expect and to see the prior that the model induces on the
#' latent intercept.
#'
#' @param x A `usdt_priors` object from [usdt_priors()].
#' @param data Optional `usdt_data` object. When given, a criterion fixed by
#'   the median split stays fixed, as in [hsdt()], and the tasks take the data
#'   labels. Without it, both criteria are estimated. The caption says which
#'   criteria are fixed.
#' @param unequal_variances Logical. Simulate the signal standard deviation
#'   from its prior instead of fixing it to 1.
#' @param ... Ignored.
#'
#' @return A `ggplot` object. Its underlying data frame is stored in `$data`.
#'
#' @details
#' The plot simulates 100,000 draws from the priors in R; no model is fitted
#' and Stan is not needed.
#'
#' The rates panel shows the hit and false-alarm rates of one new subject drawn
#' from each population draw. Rates piled up at 0 or 1 mean that the priors
#' allow extreme performance.
#'
#' The H2 panel shows the exact prior of the correlation. The difference
#' between the mean sensitivities (H1) and the intercept of the latent
#' regression (H3) have no priors of their own: the model induces them through
#' the priors of the means, standard deviations and correlation. Each is drawn
#' as the density that the prior of the indirect mean gives it, averaged over
#' 20,000 prior draws of the other parameters (a Rao-Blackwell estimate), the
#' same priors that [usdt_bf()] averages over 100,000 draws. The intercept has
#' heavy tails, so only its central 95% is drawn. Dotted lines mark the value
#' of each hypothesis under the null.
#'
#' @seealso [usdt_priors()]
#'
#' @examples
#' plot(usdt_priors())
#'
#' # A narrower prior for the mean indirect sensitivity
#' plot(usdt_priors(dprime = list(direct = "normal(0, 1)",
#'                                indirect = "normal(0, 0.5)")))
#'
#' @export
plot.usdt_priors <- function(x, data = NULL, unequal_variances = FALSE, ...) {
  if (!is.null(data) && !inherits(data, "usdt_data")) {
    .usdt_stop("`data` must come from usdt_data_long() or usdt_data_tasks(), ",
               "not a plain ", class(data)[1L], ".")
  }
  if (!is.logical(unequal_variances) || length(unequal_variances) != 1L ||
      is.na(unequal_variances)) {
    .usdt_stop("`unequal_variances` must be `TRUE` or `FALSE`.")
  }
  if (length(list(...))) {
    .usdt_stop("`...` is not used by the prior predictive plot.")
  }
  .plot_priors(.prior_check_data(x, data, unequal_variances))
}

# This function collects the curves drawn by the prior predictive check.
.prior_check_data <- function(priors, data, unequal_variances) {
  labels <- if (is.null(data)) c(direct = "Direct", indirect = "Indirect") else
    data$meta$labels
  free_c <- if (is.null(data)) c(1L, 1L) else
    as.integer(!unlist(data$meta$criterion_zero)[c("direct", "indirect")])
  simulated <- .prior_predictive(priors, .prior_samples, free_c,
                                 unequal_variances)

  # Panel titles are drawn on a graphics device, so they stay ASCII.
  panels <- c("Hit and false-alarm rates of a new subject",
              "H1: difference in mean d'",
              "H2: correlation between d'",
              "H3: intercept of the latent regression")
  curve <- function(panel, x, density, task = NA_character_,
                    rate = NA_character_) {
    data.frame(panel = panel, x = x, density = density, task = task,
               rate = rate, stringsAsFactors = FALSE)
  }

  # Rates are bounded, so their densities reflect the draws at 0 and 1.
  rates <- do.call(rbind, lapply(1:2, function(j) {
    rbind(curve(panels[1L], .unit_grid, .unit_density(simulated$hit[, j]),
                labels[[j]], "Hit rate"),
          curve(panels[1L], .unit_grid, .unit_density(simulated$fa[, j]),
                labels[[j]], "False-alarm rate"))
  }))

  # The hypotheses take the priors of the Bayes factors. The difference spans
  # four robust SDs either side of its median; the induced intercept prior has
  # heavy tails, so only its central 95% is drawn. Each point of an induced
  # density averages over prior draws, so these smooth curves take fewer
  # points and draws than a Bayes factor, which reads one point precisely.
  priors_of <- .quantity_priors(priors, c("diff", "intercept"),
                                n = .curve_samples)
  Q <- .usdt_quantities(.bivariate_primitives(simulated$draws))
  centre <- stats::median(Q[, "diff"])
  spread <- stats::IQR(Q[, "diff"]) / (2 * stats::qnorm(0.75))
  grid <- seq(centre - 4 * spread, centre + 4 * spread, length.out = 128L)
  h1 <- curve(panels[2L], grid, exp(priors_of$diff$log_density(grid)))

  grid <- seq(-0.999, 0.999, length.out = 512L)
  h2 <- curve(panels[3L], grid, exp(priors_of$rho$log_density(grid)))

  limits <- stats::quantile(Q[, "intercept"], c(0.025, 0.975), names = FALSE)
  grid <- seq(limits[1L], limits[2L], length.out = 128L)
  h3 <- curve(panels[4L], grid, exp(priors_of$intercept$log_density(grid)))

  curves <- rbind(rates, h1, h2, h3)
  curves$panel <- factor(curves$panel, levels = panels)
  curves$task <- factor(curves$task, levels = unname(labels))
  curves$rate <- factor(curves$rate,
                        levels = c("Hit rate", "False-alarm rate"))
  list(curves = curves, labels = labels,
       fixed = unname(labels[c("direct", "indirect")][free_c == 0L]),
       unequal_variances = unequal_variances,
       nulls = data.frame(panel = factor(panels[-1L], levels = panels), x = 0))
}

# The grid on which rate densities are evaluated.
.unit_grid <- seq(0, 1, length.out = 512L)

# This function estimates a density on (0, 1), reflecting the draws at both
# ends so that the mass near a bound is not lost.
.unit_density <- function(v) {
  k <- stats::density(c(v, -v, 2 - v), bw = stats::bw.nrd0(v),
                      from = 0, to = 1, n = length(.unit_grid))
  3 * k$y
}

# This function draws the prior predictive check.
.plot_priors <- function(values) {
  curves <- values$curves
  hypotheses <- curves[is.na(curves$task), , drop = FALSE]
  rates <- curves[!is.na(curves$task), , drop = FALSE]

  # The tasks keep the colours of the ROC plot, and the hypotheses those of the
  # model estimates in the shrinkage plot.
  task_colours <- stats::setNames(c("#176B60", "#B6543A"),
                                  unname(values$labels[c("direct", "indirect")]))
  model_colour <- "#116B60"
  model_fill <- "#3FA88E"

  # The rates depend on which criteria the median split fixes, so the caption
  # names them. Under equal variances a fixed criterion is exactly zero.
  fixed <- values$fixed
  where <- if (values$unequal_variances) "by the median split (HR + FAR = 1)" else
    "at 0 by the median split"
  criteria <- if (!length(fixed)) {
    "Both criteria are estimated."
  } else if (length(fixed) == 2L) {
    paste0("Both criteria are fixed ", where, ".")
  } else {
    paste0("The ", fixed, " criterion is fixed ", where,
           "; the other is estimated.")
  }

  ggplot2::ggplot(curves) +
    ggplot2::geom_vline(
      data = values$nulls, ggplot2::aes(xintercept = .data[["x"]]),
      colour = "#697076", linetype = "dotted", linewidth = 0.5
    ) +
    ggplot2::geom_area(
      data = hypotheses,
      ggplot2::aes(x = .data[["x"]], y = .data[["density"]]),
      fill = model_fill, alpha = 0.17, colour = NA
    ) +
    ggplot2::geom_line(
      data = hypotheses,
      ggplot2::aes(x = .data[["x"]], y = .data[["density"]]),
      colour = model_colour, linewidth = 0.9
    ) +
    ggplot2::geom_line(
      data = rates,
      ggplot2::aes(x = .data[["x"]], y = .data[["density"]],
                   colour = .data[["task"]], linetype = .data[["rate"]]),
      linewidth = 0.9
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data[["panel"]]), ncol = 2,
                        scales = "free") +
    # The task legend always comes first, so the legend reads the same way
    # whatever the labels.
    ggplot2::scale_colour_manual(
      values = task_colours, name = NULL,
      guide = ggplot2::guide_legend(order = 1L)
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Hit rate" = "solid", "False-alarm rate" = "dashed"),
      name = NULL, guide = ggplot2::guide_legend(order = 2L)
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.03)) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.06))
    ) +
    ggplot2::labs(
      x = NULL, y = NULL,
      caption = .wrap_caption(
        "Prior densities from 100,000 draws of the priors; no responses are ",
        "used. Rates belong to one new subject per draw. ", criteria, "\n",
        "H2 is the exact prior of the correlation. H1 and the H3 intercept are ",
        "induced by the other priors and averaged over 20,000 of the draws; ",
        "the intercept is shown over its central 95%. Dotted lines mark each ",
        "null value."
      )
    ) +
    .panel_theme() +
    # Wider keys keep the dashed false-alarm line recognisable in the legend.
    ggplot2::theme(legend.key.width = grid::unit(30, "pt"))
}
