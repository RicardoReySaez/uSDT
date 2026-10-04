# bayes-factor.R
# Bayes factors for the three hypotheses of a Bayesian hierarchical SDT model
# Author: Ricardo Rey-Sáez
# Last modified: 04-10-2026

# Public functions

#' Bayes factors for the hypotheses of a Bayesian hierarchical SDT model
#'
#' Tests point, directional and interval hypotheses about the quantities of
#' the three uSDT hypotheses in a Bayesian [hsdt()] fit: the difference between
#' mean sensitivities (`diff`, H1), their correlation (`rho`, H2) and the slope
#' and intercept of the latent regression (`slope`, `intercept`, H3).
#'
#' @param fit A Bayesian `hsdt` object, fitted with `estimation = "bayesian"`.
#' @param hypothesis Character vector of hypotheses. Each one is
#'   `"<quantity> <op> <value>"`, where the quantity is `diff`, `rho`, `slope`
#'   or `intercept` and the operator is `=`, `<` or `>` (`<=` and `>=` are read
#'   as `<` and `>`), or `"<quantity> in [a, b]"` and `"<quantity> out [a, b]"`
#'   for a region and its complement. The default tests the three point nulls.
#' @param level Credible level of the reported intervals. Defaults to the
#'   level of `fit`.
#'
#' @return An object of class `usdt_bf`: a data frame with one row per
#'   hypothesis and the columns `hypothesis`, `test`, `H1`, `H0`, `estimate`
#'   and `est.error` (posterior mean and SD), `conf.low` and `conf.high`
#'   (central credible interval), `post.prob` (posterior probability of H1 for
#'   a directional or interval test), `log_BF10`, `BF10`, `BF01`, `evidence`
#'   and `favours`. The prior and posterior curves are stored for
#'   [plot.usdt_bf()].
#'
#' @details
#' # Point hypotheses: the Savage-Dickey density ratio
#'
#' A point hypothesis \eqn{\theta = c} compares a model in which
#' \eqn{\theta} is fixed at \eqn{c} with the fitted model, in which it is free.
#' When the first is nested in the second and keeps the conditional prior of
#' the other parameters, the Bayes factor is the ratio of the posterior and
#' prior densities at \eqn{c} (Dickey & Lientz, 1970; Wagenmakers et al.,
#' 2010):
#' \deqn{\mathrm{BF}_{01} = p(\theta = c \mid y) / p(\theta = c).}
#' The posterior density comes from a logspline fit to the posterior draws.
#' Density estimates are least precise far from the posterior mass, so a test
#' value deep in a tail is read on the log scale and with caution.
#'
#' # Directional and interval hypotheses
#'
#' An inequality or a region is evidence about where the quantity lies, so its
#' Bayes factor is the ratio of posterior to prior odds of H1:
#' \deqn{\mathrm{BF}_{10} = \frac{P(H_1 \mid y) / P(H_0 \mid y)}{P(H_1) / P(H_0)}.}
#' `post.prob` reports \eqn{P(H_1 \mid y)}, read from the same logspline fit.
#'
#' # The prior of each quantity
#'
#' * `diff`: the difference of two normal means is normal, so its prior is
#'   exact.
#' * `rho`: the scaled beta of [usdt_priors()], exact.
#' * `intercept`: induced by the other priors. Given the slope
#'   \eqn{\beta_1}, the intercept is normal,
#'   \eqn{N(m_I - \beta_1 m_D, s_I^2 + \beta_1^2 s_D^2)}, so its prior density
#'   and probabilities are averages of normal ones over 100,000 prior draws of
#'   the slope (a Rao-Blackwell estimate). The prior has no mean or variance,
#'   because the slope involves \eqn{1/\sigma_D}.
#' * `slope`: given the ratio \eqn{\sigma_I / \sigma_D}, the slope is a
#'   rescaled correlation, so its prior is likewise averaged over 100,000 prior
#'   draws of that ratio. Its density is infinite at zero, and the slope always
#'   has the sign of the correlation, so `slope = 0`, `slope > 0` and
#'   `slope < 0` are tested as the same statements about `rho`. The row still
#'   summarises the posterior of the slope.
#'
#' The prior draws use the random number generator; call [set.seed()] first
#' for reproducible results.
#'
#' The Bayes factor of a point hypothesis depends on the prior of the tested
#' quantity: a wider prior puts less density at the tested value and favours
#' the null. The intercept prior is induced by the means, SDs and correlation,
#' so its Bayes factor depends on all of them. Report results with the priors
#' used, and check them with other reasonable priors.
#'
#' `evidence` follows the bands of Jeffreys (1961) with the labels of Lee and
#' Wagenmakers (2013): anecdotal (\eqn{\mathrm{BF} < 3}), moderate (< 10),
#' strong (< 30), very strong (< 100) and extreme.
#'
#' @references
#' Dickey, J. M., & Lientz, B. P. (1970). The weighted likelihood ratio, sharp
#' hypotheses about chances, the order of a Markov chain. *The Annals of
#' Mathematical Statistics*, 41(1), 214-226.
#' \doi{10.1214/aoms/1177697203}
#'
#' Jeffreys, H. (1961). *Theory of probability* (3rd ed.). Oxford University
#' Press.
#'
#' Lee, M. D., & Wagenmakers, E.-J. (2013). *Bayesian cognitive modeling: A
#' practical course*. Cambridge University Press.
#'
#' Wagenmakers, E.-J., Lodewyckx, T., Kuriyal, H., & Grasman, R. (2010).
#' Bayesian hypothesis testing for psychologists: A tutorial on the
#' Savage-Dickey method. *Cognitive Psychology*, 60(3), 158-189.
#' \doi{10.1016/j.cogpsych.2009.12.001}
#'
#' @seealso [hsdt()], [usdt_priors()], [plot.usdt_bf()]
#'
#' @examples
#' \donttest{
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
#' if (requireNamespace("rstan", quietly = TRUE) &&
#'     requireNamespace("logspline", quietly = TRUE)) {
#'   m <- hsdt(d, estimation = "bayesian", seed = 1)
#'
#'   # The three point nulls
#'   set.seed(1)
#'   usdt_bf(m)
#'
#'   # Unconscious processing as a directional hypothesis, and a region of
#'   # practical equivalence for the difference
#'   b <- usdt_bf(m, c("intercept > 0", "diff in [-0.1, 0.1]"))
#'   b
#'   plot(b)
#' }
#' }
#'
#' @export
usdt_bf <- function(fit, hypothesis = c("diff = 0", "rho = 0", "intercept = 0"),
                    level = fit$level) {
  if (!.is_bayes(fit)) {
    .usdt_stop("`fit` must come from hsdt(estimation = \"bayesian\"), whose ",
               "priors uSDT knows.")
  }
  if (!is.character(hypothesis) || !length(hypothesis) || anyNA(hypothesis)) {
    .usdt_stop("`hypothesis` must be one or more strings such as \"rho = 0\".")
  }
  .check_confidence_level(level)
  rlang::check_installed("logspline", reason = "to estimate posterior densities.")

  tests <- lapply(hypothesis, .bf_parse)
  draws <- .draw_quantities(.stan_bivariate(fit$draws))
  priors <- .bf_priors(fit$priors,
                       unique(vapply(tests, `[[`, "", "quantity")))
  results <- lapply(tests, function(h) {
    .bf_one(h, draws[, h$quantity], draws[, h$summary], priors[[h$quantity]],
            level)
  })

  out <- do.call(rbind, lapply(results, `[[`, "row"))
  rownames(out) <- NULL
  structure(out, curves = lapply(results, `[[`, "curve"), level = level,
            class = c("usdt_bf", "data.frame"))
}

#' @param x A `usdt_bf` object.
#' @param ... Ignored.
#' @rdname usdt_bf
#' @export
print.usdt_bf <- function(x, ...) {
  cat(.rule("Bayes factors"), "\n")
  for (i in seq_len(nrow(x))) {
    r <- x[i, ]
    cat(sprintf("\n  %s (%s)\n", r$hypothesis, r$test))
    cat(sprintf("    H1: %s   vs   H0: %s\n", r$H1, r$H0))
    cat(sprintf("    Mean %s, SD %s, %.0f%% CrI %s\n",
                .fmt_n(r$estimate, 4L, 0L), .fmt_n(r$est.error, 4L, 0L),
                100 * attr(x, "level"), .fmt_ci(r$conf.low, r$conf.high)))
    if (!is.na(r$post.prob)) {
      cat(sprintf("    P(H1 | data) = %s\n",
                  formatC(r$post.prob, format = "f", digits = 3)))
    }
    cat(sprintf("    log BF10 = %s (BF10 = %s, BF01 = %s): %s evidence for %s\n",
                formatC(r$log_BF10, format = "f", digits = 2),
                formatC(r$BF10, format = "g", digits = 3),
                formatC(r$BF01, format = "g", digits = 3),
                r$evidence, r$favours))
  }
  invisible(x)
}

# Internal functions

# This function reads one hypothesis string.
.bf_parse <- function(text) {
  quantities <- c("diff", "rho", "slope", "intercept")
  number <- "([-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"
  region <- regmatches(text, regexec(paste0(
    "^\\s*([a-z]+)\\s+(in|out)\\s*[\\[(]\\s*", number, "\\s*,\\s*", number,
    "\\s*[\\])]\\s*$"), text, perl = TRUE))[[1L]]
  compare <- regmatches(text, regexec(paste0(
    "^\\s*([a-z]+)\\s*(==|=|<=|>=|<|>)\\s*", number, "\\s*$"), text,
    perl = TRUE))[[1L]]

  if (length(region)) {
    h <- list(quantity = region[2L], op = region[3L],
              value = sort(as.numeric(region[4:5])))
    if (h$value[1L] == h$value[2L]) {
      .usdt_stop("the region of \"", text, "\" has no width.")
    }
  } else if (length(compare)) {
    op <- switch(compare[3L], "==" = "=", "<=" = "<", ">=" = ">", compare[3L])
    h <- list(quantity = compare[2L], op = op, value = as.numeric(compare[4L]))
  } else {
    .usdt_stop("cannot read the hypothesis \"", text, "\". Write it as, e.g., ",
               "\"rho = 0\", \"intercept > 0\" or \"diff in [-0.1, 0.1]\".")
  }
  if (!h$quantity %in% quantities) {
    .usdt_stop("the hypothesis \"", text, "\" tests `", h$quantity, "`; the ",
               "quantities are ", paste0("`", quantities, "`", collapse = ", "),
               ".")
  }
  if (h$quantity == "rho" && any(abs(h$value) > 1) ||
      h$quantity == "rho" && h$op == "=" && abs(h$value) == 1) {
    .usdt_stop("a correlation can only be tested inside (-1, 1).")
  }
  h$text <- text

  # The slope has the sign of the correlation, and its prior is infinite at
  # zero, so a hypothesis that compares the slope with zero is tested on the
  # correlation, while the row still summarises the slope.
  h$summary <- h$quantity
  h$as_rho <- h$quantity == "slope" && h$op %in% c("=", "<", ">") &&
    h$value == 0
  if (h$as_rho) h$quantity <- "rho"
  h
}

# This function returns the log density and the distribution function of the
# prior of each requested quantity.
.bf_priors <- function(priors, quantities) {
  d <- .prior_rows(priors, "dprime")
  r <- .prior_rows(priors, "cor_dprime", "both")
  rho_log_density <- function(x) {
    stats::dbeta((x + 1) / 2, r$alpha, r$beta, log = TRUE) - log(2)
  }
  rho_cdf <- function(x) stats::pbeta((x + 1) / 2, r$alpha, r$beta)
  centre <- d$location[2L] - d$location[1L]
  spread <- sqrt(sum(d$scale^2))
  out <- list(
    diff = list(log_density = function(x) stats::dnorm(x, centre, spread, log = TRUE),
                cdf = function(x) stats::pnorm(x, centre, spread), bounds = NULL),
    rho = list(log_density = rho_log_density, cdf = rho_cdf, bounds = c(-1, 1))
  )

  # The slope and intercept priors are averaged over prior draws of the
  # quantities they depend on (Rao-Blackwell).
  if (any(c("intercept", "slope") %in% quantities)) {
    B <- .prior_draws(priors, .prior_samples)
    ratio <- B[, "sigma_I"] / B[, "sigma_D"]
    slope <- B[, "rho"] * ratio
    location <- d$location[2L] - slope * d$location[1L]
    scale <- sqrt(d$scale[2L]^2 + slope^2 * d$scale[1L]^2)
    out$intercept <- list(
      log_density = function(x) vapply(x, function(v) {
        .log_mean_exp(stats::dnorm(v, location, scale, log = TRUE))
      }, 0),
      cdf = function(x) vapply(x, function(v) {
        mean(stats::pnorm(v, location, scale))
      }, 0),
      bounds = NULL)
    out$slope <- list(
      log_density = function(x) vapply(x, function(v) {
        .log_mean_exp(rho_log_density(v / ratio) - log(ratio))
      }, 0),
      cdf = function(x) vapply(x, function(v) mean(rho_cdf(v / ratio)), 0),
      bounds = NULL)
  }
  out
}

# This function averages probabilities given on the log scale.
.log_mean_exp <- function(l) {
  top <- max(l)
  if (!is.finite(top)) return(top)
  top + log(mean(exp(l - top)))
}

# This function fits a logspline density to posterior draws, within the bounds
# of the quantity, and falls back to the older algorithm when it fails.
.bf_posterior <- function(x, bounds = NULL) {
  args <- list(x)
  if (!is.null(bounds)) args <- c(args, lbound = bounds[1L], ubound = bounds[2L])
  quiet <- function(f) {
    out <- NULL
    utils::capture.output(out <- suppressWarnings(
      tryCatch(do.call(f, args), error = function(e) NULL)))
    out
  }
  fit <- quiet(logspline::logspline)
  if (!is.null(fit)) {
    return(list(log_density = function(q) logspline::dlogspline(q, fit, log = TRUE),
                cdf = function(q) logspline::plogspline(q, fit)))
  }
  fit <- quiet(logspline::oldlogspline)
  if (is.null(fit)) {
    .usdt_stop("the posterior density could not be estimated from the draws.")
  }
  list(log_density = function(q) log(logspline::doldlogspline(q, fit)),
       cdf = function(q) logspline::poldlogspline(q, fit))
}

# This function evaluates one hypothesis on the draws of the tested quantity,
# and summarises the draws of the quantity it names.
.bf_one <- function(h, draws, summary, prior, level) {
  post <- .bf_posterior(draws, prior$bounds)
  a <- (1 - level) / 2
  limits <- stats::quantile(summary, c(a, 1 - a), names = FALSE)
  name <- h$summary
  as_rho <- if (h$as_rho) sprintf(", as rho %s 0", h$op) else ""

  if (h$op == "=") {
    test <- paste0("Savage-Dickey", as_rho)
    h1 <- sprintf("%s != %s", name, format(h$value))
    h0 <- sprintf("%s = %s", name, format(h$value))
    log_bf10 <- prior$log_density(h$value) - post$log_density(h$value)
    prob <- NA_real_
  } else {
    region <- function(cdf) switch(h$op,
      ">" = 1 - cdf(h$value),
      "<" = cdf(h$value),
      "in" = cdf(h$value[2L]) - cdf(h$value[1L]),
      out = 1 - (cdf(h$value[2L]) - cdf(h$value[1L])))
    interval <- sprintf("[%s, %s]", format(h$value[1L]), format(h$value[2L]))
    test <- paste0(if (h$op %in% c("<", ">")) "directional" else "interval",
                   as_rho)
    h1 <- switch(h$op, ">" = sprintf("%s > %s", name, format(h$value)),
                 "<" = sprintf("%s < %s", name, format(h$value)),
                 sprintf("%s %s %s", name, h$op, interval))
    h0 <- switch(h$op, ">" = sprintf("%s <= %s", name, format(h$value)),
                 "<" = sprintf("%s >= %s", name, format(h$value)),
                 sprintf("%s %s %s", name, if (h$op == "in") "out" else "in",
                         interval))
    prob <- region(post$cdf)
    log_bf10 <- stats::qlogis(prob) - stats::qlogis(region(prior$cdf))
  }

  evidence <- as.character(cut(abs(log_bf10), log(c(0, 3, 10, 30, 100, Inf)),
                               c("anecdotal", "moderate", "strong",
                                 "very strong", "extreme"),
                               right = FALSE, include.lowest = TRUE))
  row <- data.frame(
    hypothesis = h$text, test = test, H1 = h1, H0 = h0,
    estimate = mean(summary), est.error = stats::sd(summary),
    conf.low = limits[1L], conf.high = limits[2L], post.prob = prob,
    log_BF10 = log_bf10, BF10 = exp(log_bf10), BF01 = exp(-log_bf10),
    evidence = evidence, favours = if (log_bf10 >= 0) "H1" else "H0",
    stringsAsFactors = FALSE)

  # The plotted range covers the central 95% of the posterior and every tested
  # value; heavy posterior tails would otherwise squeeze the curves.
  range <- range(c(stats::quantile(draws, c(0.025, 0.975), names = FALSE),
                   h$value))
  pad <- 0.25 * diff(range)
  grid <- seq(range[1L] - pad, range[2L] + pad, length.out = 256L)
  if (!is.null(prior$bounds)) {
    grid <- grid[grid > prior$bounds[1L] & grid < prior$bounds[2L]]
  }
  curve <- list(hypothesis = paste0(h$text, if (h$as_rho) " (as rho)" else ""),
                op = h$op, value = h$value, grid = grid,
                prior = exp(prior$log_density(grid)),
                posterior = exp(post$log_density(grid)),
                at = if (h$op == "=") {
                  exp(c(prior = prior$log_density(h$value),
                        posterior = post$log_density(h$value)))
                },
                log_bf10 = log_bf10, evidence = evidence,
                favours = row$favours, prob = prob)
  list(row = row, curve = curve)
}

# Prior and posterior plot

#' Plot the prior and posterior of tested quantities
#'
#' Draws, for each hypothesis of a [usdt_bf()] result, the posterior density
#' (solid) and the prior density (dashed) of the tested quantity. For a point
#' hypothesis, a dotted line marks the tested value and two points mark the
#' heights whose ratio is the Savage-Dickey Bayes factor. For a directional
#' or interval hypothesis, the shaded area is the posterior probability of H1.
#'
#' @param x A `usdt_bf` object from [usdt_bf()].
#' @param ... Ignored.
#'
#' @return A `ggplot` object. Its underlying data frame is stored in `$data`.
#'
#' @seealso [usdt_bf()]
#'
#' @export
plot.usdt_bf <- function(x, ...) {
  if (length(list(...))) {
    .usdt_stop("`...` is not used by the Bayes factor plot.")
  }
  curves <- attr(x, "curves")
  panels <- vapply(curves, `[[`, "", "hypothesis")
  panel <- function(i, n) factor(rep(panels[i], n), levels = panels)

  densities <- do.call(rbind, lapply(seq_along(curves), function(i) {
    cv <- curves[[i]]
    n <- length(cv$grid)
    data.frame(panel = panel(i, 2L * n), x = rep(cv$grid, 2L),
               density = c(cv$posterior, cv$prior),
               distribution = factor(rep(c("Posterior", "Prior"), each = n),
                                     levels = c("Posterior", "Prior")))
  }))

  # The shaded pieces are the posterior mass of H1: one for an inequality or
  # a region, two for the outside of a region.
  areas <- do.call(rbind, lapply(seq_along(curves), function(i) {
    cv <- curves[[i]]
    pieces <- switch(cv$op,
      ">" = list(c(cv$value, Inf)), "<" = list(c(-Inf, cv$value)),
      "in" = list(cv$value),
      out = list(c(-Inf, cv$value[1L]), c(cv$value[2L], Inf)),
      list())
    do.call(rbind, lapply(seq_along(pieces), function(k) {
      inside <- cv$grid >= pieces[[k]][1L] & cv$grid <= pieces[[k]][2L]
      if (sum(inside) < 2L) return(NULL)
      data.frame(panel = panel(i, sum(inside)), x = cv$grid[inside],
                 density = cv$posterior[inside],
                 piece = paste(i, k))
    }))
  }))
  values <- do.call(rbind, lapply(seq_along(curves), function(i) {
    data.frame(panel = panel(i, length(curves[[i]]$value)),
               x = curves[[i]]$value)
  }))
  points <- do.call(rbind, lapply(seq_along(curves), function(i) {
    cv <- curves[[i]]
    if (is.null(cv$at)) return(NULL)
    data.frame(panel = panel(i, 2L), x = cv$value, density = cv$at[c(2L, 1L)],
               distribution = factor(c("Posterior", "Prior"),
                                     levels = c("Posterior", "Prior")))
  }))
  labels <- do.call(rbind, lapply(seq_along(curves), function(i) {
    cv <- curves[[i]]
    data.frame(panel = panel(i, 1L), x = Inf, y = Inf, label = paste0(
      sprintf("log BF10 = %.2f\n%s, favours %s", cv$log_bf10, cv$evidence,
              cv$favours),
      if (!is.na(cv$prob)) sprintf("\nP(H1 | data) = %.3f", cv$prob) else ""))
  }))

  # The posterior takes the colours of the model estimates and the prior the
  # grey of what does not come from the fit.
  colours <- c(Posterior = "#116B60", Prior = "#A8AEB3")
  plot <- ggplot2::ggplot(densities)
  if (!is.null(areas)) {
    plot <- plot + ggplot2::geom_ribbon(
      data = areas,
      ggplot2::aes(x = .data[["x"]], ymin = 0, ymax = .data[["density"]],
                   group = .data[["piece"]]),
      inherit.aes = FALSE, fill = "#3FA88E", alpha = 0.17
    )
  }
  plot <- plot +
    ggplot2::geom_vline(
      data = values, ggplot2::aes(xintercept = .data[["x"]]),
      colour = "#697076", linetype = "dotted", linewidth = 0.5
    ) +
    ggplot2::geom_line(
      ggplot2::aes(x = .data[["x"]], y = .data[["density"]],
                   colour = .data[["distribution"]],
                   linetype = .data[["distribution"]]),
      linewidth = 0.9
    )
  if (!is.null(points)) {
    plot <- plot + ggplot2::geom_point(
      data = points,
      ggplot2::aes(x = .data[["x"]], y = .data[["density"]],
                   fill = .data[["distribution"]]),
      inherit.aes = FALSE, shape = 21, size = 2.8, stroke = 0.8,
      colour = "white", show.legend = FALSE
    )
  }
  plot +
    ggplot2::geom_text(
      data = labels,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                   label = .data[["label"]]),
      inherit.aes = FALSE, hjust = 1.05, vjust = 1.2, lineheight = 1.05,
      size = 3.35, colour = "#3E4347"
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data[["panel"]]), scales = "free",
                        ncol = min(length(curves), 2L)) +
    ggplot2::scale_colour_manual(values = colours, name = NULL) +
    ggplot2::scale_fill_manual(values = colours, name = NULL) +
    ggplot2::scale_linetype_manual(
      values = c(Posterior = "solid", Prior = "dashed"), name = NULL
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.03)) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.25))
    ) +
    ggplot2::labs(
      x = NULL, y = NULL,
      caption = .wrap_caption(
        "Posterior (solid) and prior (dashed) densities of each tested ",
        "quantity. For a point hypothesis the Bayes factor is the ratio of ",
        "the two marked heights at the dotted line (Savage-Dickey); for a ",
        "directional or interval hypothesis the shaded area is the posterior ",
        "probability of H1. The intercept and slope priors are induced by the ",
        "other priors."
      )
    ) +
    .panel_theme() +
    ggplot2::theme(legend.key.width = grid::unit(30, "pt"))
}
