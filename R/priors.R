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
#' * between-subject standard deviations: `"student_t(df, location, scale)"`,
#'   a Student-\eqn{t} distribution truncated at zero. With location 0 it is a
#'   half-Student-\eqn{t};
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
#' # Standard deviations of d' expected around 0.3 rather than near zero
#' usdt_priors(sd_dprime = "student_t(4, 0.3, 0.2)")
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
  free <- c("direct", "indirect")[free_c == 1L]
  d <- .prior_rows(priors, "dprime")
  sd_d <- .prior_rows(priors, "sd_dprime")
  cor_d <- .prior_rows(priors, "cor_dprime", "both")
  c_ <- .prior_rows(priors, "criterion", free)
  sd_c <- .prior_rows(priors, "sd_criterion", free)
  cor_c <- .prior_rows(priors, "cor_criterion", "both")
  s <- .prior_rows(priors, "sd_signal",
                   if (UV) c("direct", "indirect") else character(0))
  list(mu_d_loc = vec(d$location), mu_d_scale = vec(d$scale),
       sd_d_df = vec(sd_d$df), sd_d_loc = vec(sd_d$location),
       sd_d_scale = vec(sd_d$scale),
       rho_d_a = cor_d$alpha, rho_d_b = cor_d$beta,
       mu_c_loc = vec(c_$location), mu_c_scale = vec(c_$scale),
       sd_c_df = vec(sd_c$df), sd_c_loc = vec(sd_c$location),
       sd_c_scale = vec(sd_c$scale),
       rho_c_a = cor_c$alpha, rho_c_b = cor_c$beta,
       sd_s_loc = vec(s$location), sd_s_scale = vec(s$scale))
}

# Number of draws that approximate a prior the model induces rather than
# states, such as the prior of the latent intercept.
.prior_samples <- 1e5L

# This function draws n values from one prior. The Student-t is truncated at
# zero, so its draws invert the distribution function above that point.
.r_prior <- function(row, n) {
  switch(row$family,
         normal = stats::rnorm(n, row$location, row$scale),
         student_t = row$location + row$scale * stats::qt(
           stats::runif(n, stats::pt(-row$location / row$scale, row$df), 1),
           row$df),
         scaled_beta = 2 * stats::rbeta(n, row$alpha, row$beta) - 1,
         lognormal = stats::rlnorm(n, row$location, row$scale))
}

# This function simulates the population distribution of the two sensitivities
# from its priors, one draw per row.
.prior_draws <- function(priors, n) {
  d <- .prior_rows(priors, "dprime")
  s <- .prior_rows(priors, "sd_dprime")
  cbind(mu_D = .r_prior(d[1L, ], n), mu_I = .r_prior(d[2L, ], n),
        sigma_D = .r_prior(s[1L, ], n), sigma_I = .r_prior(s[2L, ], n),
        rho = .r_prior(.prior_rows(priors, "cor_dprime", "both"), n))
}

# This function adds one new subject to every population draw and returns that
# subject's hit and false-alarm rates in each task, following the subject level
# of the Stan model. A criterion fixed by the Meyen split takes the value that
# makes HR + FAR = 1.
.prior_predictive <- function(priors, n, free_c, UV) {
  B <- .prior_draws(priors, n)
  z <- matrix(stats::rnorm(2L * n), n)
  d <- cbind(B[, "mu_D"] + B[, "sigma_D"] * z[, 1L],
             B[, "mu_I"] + B[, "sigma_I"] *
               (B[, "rho"] * z[, 1L] + sqrt(1 - B[, "rho"]^2) * z[, 2L]))
  signal <- .prior_rows(priors, "sd_signal")
  sigma_s <- if (UV) {
    vapply(1:2, function(j) .r_prior(signal[j, ], n), numeric(n))
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
      r <- .r_prior(.prior_rows(priors, "cor_criterion", "both"), n)
      w[, 2L] <- r * w[, 1L] + sqrt(1 - r^2) * w[, 2L]
    }
    for (k in seq_along(free)) {
      j <- free[k]
      criterion[, j] <- .r_prior(mean_c[j, ], n) +
        .r_prior(sd_c[j, ], n) * w[, k]
    }
  }

  list(draws = B,
       hit = stats::pnorm((d / 2 - criterion) / sigma_s),
       fa = stats::pnorm(-d / 2 - criterion))
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
#'   labels. Without it, both criteria are estimated.
#' @param unequal_variances Logical. Simulate the signal standard deviation
#'   from its prior instead of fixing it to 1.
#' @param ... Ignored.
#'
#' @return A `ggplot` object. Its underlying data frame is stored in `$data`.
#'
#' @details
#' The plot simulates 100,000 draws from the priors in R; no model is fitted
#' and Stan is not needed. These draws match those of the Stan model sampling
#' from its priors alone.
#'
#' The rates panel shows the hit and false-alarm rates of one new subject drawn
#' from each population draw. Rates piled up at 0 or 1 mean that the priors
#' allow extreme performance.
#'
#' The H1 and H2 panels show exact prior densities: the difference between two
#' normal means is normal, and the correlation follows its scaled beta. The H3
#' panel shows the intercept of the latent regression, whose prior the model
#' induces through the means, standard deviations and correlation. It is
#' estimated from the draws and drawn over its central 95%, because its tails
#' are heavy. Dotted lines mark the value of each hypothesis under the null.
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

  # The difference of two normal means and the scaled beta are exact.
  d <- .prior_rows(priors, "dprime")
  centre <- d$location[2L] - d$location[1L]
  spread <- sqrt(sum(d$scale^2))
  grid <- seq(centre - 4 * spread, centre + 4 * spread, length.out = 512L)
  h1 <- curve(panels[2L], grid, stats::dnorm(grid, centre, spread))

  r <- .prior_rows(priors, "cor_dprime", "both")
  grid <- seq(-0.999, 0.999, length.out = 512L)
  h2 <- curve(panels[3L], grid,
              stats::dbeta((grid + 1) / 2, r$alpha, r$beta) / 2)

  # The induced intercept prior has heavy tails, so only its central 95% is
  # drawn.
  intercept <- .usdt_quantities(
    .bivariate_primitives(simulated$draws))[, "intercept"]
  limits <- stats::quantile(intercept, c(0.025, 0.975), names = FALSE)
  kernel <- stats::density(intercept, from = limits[1L], to = limits[2L],
                           n = 512L)
  h3 <- curve(panels[4L], kernel$x, kernel$y)

  curves <- rbind(rates, h1, h2, h3)
  curves$panel <- factor(curves$panel, levels = panels)
  curves$task <- factor(curves$task, levels = unname(labels))
  curves$rate <- factor(curves$rate,
                        levels = c("Hit rate", "False-alarm rate"))
  list(curves = curves, labels = labels,
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
    ggplot2::scale_colour_manual(values = task_colours, name = NULL) +
    ggplot2::scale_linetype_manual(
      values = c("Hit rate" = "solid", "False-alarm rate" = "dashed"),
      name = NULL
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.03)) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.06))
    ) +
    ggplot2::labs(
      x = NULL, y = NULL,
      caption = .wrap_caption(
        "Prior densities from 100,000 draws of the priors; no data are used. ",
        "Rates belong to one new subject per draw. H1 and H2 are exact ",
        "densities; the H3 intercept is induced by the other priors and is ",
        "shown over its central 95%. Dotted lines mark each null value."
      )
    ) +
    .panel_theme() +
    # Wider keys keep the dashed false-alarm line recognisable in the legend.
    ggplot2::theme(legend.key.width = grid::unit(30, "pt"))
}
