# bayes-factor.R
# Bayes factors for the three hypotheses of a Bayesian hierarchical SDT model
# Author: Ricardo Rey-Sáez
# Last modified: 27-10-2026

# Public functions

#' Bayes factors for the hypotheses of a Bayesian hierarchical SDT model
#'
#' Tests point, directional and interval hypotheses about the quantities of
#' the three uSDT hypotheses in a Bayesian [hsdt()] fit or a custom model
#' fitted with brms: the difference between mean sensitivities (`diff`, H1,
#' reported as \eqn{\Delta d'}, indirect minus direct), their correlation
#' (`rho`, H2) and the slope and intercept of the latent regression (`slope`,
#' `intercept`, H3).
#'
#' @param fit A Bayesian `hsdt` object, fitted with `estimation = "bayesian"`,
#'   or a `brmsfit` whose priors meet the conditions in Details.
#' @param hypothesis Character vector of hypotheses. Each one is
#'   `"<quantity> <op> <value>"`, where the quantity is `diff`, `rho`, `slope`
#'   or `intercept` and the operator is `=`, `<` or `>` (`<=` and `>=` are read
#'   as `<` and `>`), or `"<quantity> in [a, b]"` and `"<quantity> out [a, b]"`
#'   for a region and its complement. The default tests the three point nulls.
#' @param direct,indirect For a brms model, the population-level terms of the
#'   direct and indirect sensitivities, as in [usdt_tests()]. Ignored for an
#'   `hsdt` fit.
#' @param level Credible level of the reported intervals. Defaults to the
#'   level of an `hsdt` fit, and to 0.95 for a brms model.
#' @param plot Logical. Draw, for each hypothesis, a panel with the prior and
#'   posterior of the tested quantity (default `TRUE`). Set `plot = FALSE`
#'   inside loops, reports and other functions, and draw the figure later with
#'   [plot.usdt_bf()].
#'
#' @return An object of class `usdt_bf`: a data frame with one row per
#'   hypothesis and the columns `hypothesis`, `test` (`"Savage-Dickey density
#'   ratio test"`, `"Directional test"` or `"Interval test"`), `H1`, `H0`,
#'   `estimate` and `est.error` (posterior mean and SD), `conf.low` and
#'   `conf.high` (central credible interval), `post.prob` (posterior
#'   probability of H1 for a directional or interval test), `log_BF10`,
#'   `BF10`, `BF01`, `evidence` and `favours`. It prints one report per
#'   hypothesis. The prior and posterior curves are stored, so
#'   [plot.usdt_bf()] can draw the figure again.
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
#' A fit is used only if its distribution function stays within 0.05 of that
#' of the draws at every percentile. When neither logspline algorithm gives
#' such a fit, which happens with the very heavy tails of the intercept and
#' slope when the direct sensitivities barely vary across subjects, the draws
#' are fitted on the scale
#' \eqn{\mathrm{asinh}((x - m)/s)}, with \eqn{m} their median and \eqn{s}
#' their median absolute deviation, and the density is mapped back.
#' Density estimates are least precise far from the posterior mass, so a test
#' value deep in a tail is read on the log scale and with caution.
#'
#' # Directional and interval hypotheses
#'
#' An inequality or a region is evidence about where the quantity lies, so its
#' Bayes factor is the ratio of posterior to prior odds of H1:
#' \deqn{\mathrm{BF}_{10} = \frac{P(H_1 \mid y) / P(H_0 \mid y)}{P(H_1) / P(H_0)}.}
#' `post.prob` reports \eqn{P(H_1 \mid y)}, the share of posterior draws in
#' H1, and the prior probabilities come from the priors described below. When
#' every draw falls on one side, the draws cannot bound the Bayes factor and
#' `log_BF10` is infinite.
#'
#' # The prior of each quantity
#'
#' * `rho`: the marginal of the `lkj` prior of [usdt_priors()], exact.
#' * `diff`: given the mean direct sensitivity \eqn{\mu_D}, the difference is
#'   the mean indirect sensitivity shifted by \eqn{\mu_D}, so its prior density
#'   and probabilities are those of the prior of \eqn{\mu_I} at
#'   \eqn{x + \mu_D}, averaged over 100,000 prior draws of \eqn{\mu_D} (a
#'   Rao-Blackwell estimate). With normal priors on both means this is the
#'   normal prior of their difference.
#' * `intercept`: induced by the other priors. Given the slope
#'   \eqn{\beta_1} and \eqn{\mu_D}, the intercept is \eqn{\mu_I} shifted by
#'   \eqn{\beta_1 \mu_D}, so its prior is likewise averaged over 100,000 prior
#'   draws of the slope and of \eqn{\mu_D}. The prior has no mean or variance,
#'   because the slope involves \eqn{1/\sigma_D}.
#' * `slope`: given the ratio \eqn{\sigma_I / \sigma_D}, the slope is a
#'   rescaled correlation, so its prior is averaged over 100,000 prior draws of
#'   that ratio. Its density is infinite at zero, and the slope always has the
#'   sign of the correlation, so `slope = 0`, `slope > 0` and `slope < 0` are
#'   tested as the same statements about `rho`. The row still summarises the
#'   posterior of the slope.
#'
#' Under the latent regression parameterization the intercept has a prior of
#' its own instead of \eqn{\mu_I}, so it is exact like `rho`, and `slope` is
#' as above. The difference is induced:
#' * `diff`: \eqn{\beta_0 + (\beta_1 - 1) \mu_D}, so its prior is that of the
#'   intercept, shifted and averaged over 100,000 prior draws of the slope and
#'   of \eqn{\mu_D}. Like the induced intercept above, it has heavy tails and
#'   no mean or variance.
#'
#' The data fix \eqn{\mu_I = \beta_0 + \beta_1 \mu_D}, so under the latent
#' regression a narrow intercept prior also moves the posterior of the slope
#' and of `rho`, even though their priors are independent of it: the Bayes
#' factors of `rho` and `slope` then depend on the intercept prior as well.
#'
#' These averages work for any family of [usdt_priors()]. The prior draws use
#' the random number generator; call [set.seed()] first for reproducible
#' results.
#'
#' # Custom brms models
#'
#' For a `brmsfit`, the priors are read from [brms::prior_summary()], each
#' parameter taking its most specific prior (coefficient, then group, then
#' class), and the quantities get the same priors as above. The model must
#' have:
#' * a proper prior on each sensitivity term, from the families of
#'   [usdt_priors()], e.g. `prior(normal(0, 1), class = b)`. The flat default
#'   of brms is improper, and a Bayes factor needs a proper prior;
#' * an LKJ(\eqn{\eta}) prior on the correlations of the group-level term that
#'   holds both sensitivities. One correlation of a \eqn{K \times K} LKJ matrix
#'   has the marginal of a \eqn{2 \times 2} LKJ(\eqn{\eta - 1 + K/2}), so a
#'   block with other terms besides the two sensitivities still has an exact
#'   prior for `rho`;
#' * for `intercept` and `slope`, priors on the two SDs from the families of
#'   [usdt_priors()], truncated at zero as brms does. The brms default,
#'   `student_t(3, 0, 2.5)`, qualifies.
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
#' m <- hsdt(d, estimation = "bayesian", seed = 1)
#'
#' # The three point nulls, with their figure
#' set.seed(1)
#' usdt_bf(m)
#'
#' # Unconscious processing as a directional hypothesis, and a region of
#' # practical equivalence for the difference, without the figure
#' b <- usdt_bf(m, c("intercept > 0", "diff in [-0.1, 0.1]"), plot = FALSE)
#' b
#' plot(b[1, ])
#' }
#'
#' @export
usdt_bf <- function(fit, hypothesis = c("diff = 0", "rho = 0", "intercept = 0"),
                    direct = "d_D", indirect = "d_I", level = NULL,
                    plot = TRUE) {
  if (!.is_posterior(fit)) {
    .usdt_stop("`fit` must come from hsdt(estimation = \"bayesian\") or from ",
               "brms.")
  }
  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    .usdt_stop("`plot` must be TRUE or FALSE.")
  }
  if (!is.character(hypothesis) || !length(hypothesis) || anyNA(hypothesis)) {
    .usdt_stop("`hypothesis` must be one or more strings such as \"rho = 0\".")
  }
  level <- level %||% if (.is_bayes(fit)) fit$level else 0.95
  .check_confidence_level(level)

  parameterization <- if (.is_bayes(fit)) fit$design$parameterization else
    "correlation"
  tests <- lapply(hypothesis, .bf_parse)
  quantities <- unique(vapply(tests, `[[`, "", "quantity"))
  draws <- .draw_quantities(.bivariate_draws(fit, direct, indirect))
  priors <- if (.is_bayes(fit)) fit$priors else {
    .brms_priors(brms::prior_summary(fit), posterior::variables(fit), direct,
                 indirect)
  }
  priors <- .quantity_priors(priors, quantities, parameterization)
  results <- lapply(tests, function(h) {
    .bf_one(h, draws[, h$quantity], draws[, h$summary], priors[[h$quantity]],
            level)
  })

  out <- do.call(rbind, lapply(results, `[[`, "row"))
  rownames(out) <- NULL
  out <- structure(out, curves = lapply(results, `[[`, "curve"),
                   level = level, class = c("usdt_bf", "data.frame"))

  # A hypothesis is read from its picture as much as from its row, so the
  # figure is drawn unless it is turned off.
  if (plot) print(plot.usdt_bf(out))
  out
}

#' @param x A `usdt_bf` object.
#' @param digits Number of decimals printed.
#' @param width Width of the printed report, in characters.
#' @param ... Ignored.
#' @rdname usdt_bf
#' @export
print.usdt_bf <- function(x, digits = 3L, width = 80L, ...) {
  # A subset that lost columns is printed as the data frame it now is.
  needed <- c("test", "H1", "H0", "estimate", "est.error", "conf.low",
              "conf.high", "post.prob", "log_BF10", "evidence", "favours")
  if (!all(needed %in% names(x))) return(NextMethod())

  width <- max(60L, as.integer(width))
  chars <- .usdt_chars()
  level <- attr(x, "level")
  interval <- if (is.null(level)) "CrI" else sprintf("%g%% CrI", 100 * level)
  number <- function(v) {
    if (is.na(v)) return("NA")
    if (is.infinite(v)) return(if (v > 0) "Inf" else "-Inf")
    formatC(v, format = "f", digits = digits)
  }
  row <- function(label, value) cat(sprintf("  %-27s : %s\n", label, value))

  # Every hypothesis gets a complete card, so mixed test types stay explicit.
  for (i in seq_len(nrow(x))) {
    if (i > 1L) cat("\n")
    cat(" uSDT hypothesis ", chars$h, " ", x$test[[i]], "\n", sep = "")
    cat(strrep(chars$hh, width), "\n", sep = "")
    row("Hypothesis (H1)", x$H1[[i]])
    row("Hypothesis (H0)", x$H0[[i]])
    row("Posterior mean", number(x$estimate[[i]]))
    row("Posterior SD", number(x$est.error[[i]]))
    row(interval, sprintf("[%s, %s]", number(x$conf.low[[i]]),
                          number(x$conf.high[[i]])))
    if (!is.na(x$post.prob[[i]])) {
      row("Posterior probability (H1)", number(x$post.prob[[i]]))
    }
    row("log BF10", number(x$log_BF10[[i]]))
    row("Evidence", paste(x$evidence[[i]], "for", x$favours[[i]]))
    cat(strrep(chars$h, width), "\n", sep = "")
  }
  invisible(x)
}

# Internal functions

# This function reads one hypothesis string for a Bayes factor.
.bf_parse <- function(text) {
  h <- .parse_hypothesis(text)

  # The slope has the sign of the correlation and an infinite prior density at
  # zero, so a comparison of the slope with zero is tested on the correlation,
  # which has a prior of its own. The row still summarises the slope.
  h$summary <- h$quantity
  h$swapped <- h$quantity == "slope" && h$op %in% c("=", "<", ">") &&
    h$value == 0
  if (h$swapped) h$quantity <- "rho"
  h
}

# This function translates the priors of a brms model into those of uSDT, so
# that the tested quantities get the same priors as a Bayesian hsdt fit. Each
# parameter takes its most specific prior, as brms does: that of its
# coefficient, then of its group, then of its class. `table` is
# brms::prior_summary() and `variables` the names of the model's draws.
.brms_priors <- function(table, variables, direct, indirect) {
  found <- .brms_variables(variables, direct, indirect)
  terms <- c(direct = direct, indirect = indirect)
  group <- substring(found[3L], 4L, nchar(found[3L]) - nchar(direct) - 2L)
  table <- as.data.frame(table)[c("prior", "class", "coef", "group")]
  table[is.na(table)] <- ""

  effective <- function(class, coef, group) {
    rows <- table[table$class %in% class, , drop = FALSE]
    for (level in list(c(coef, group), c("", group), c("", ""))) {
      prior <- rows$prior[rows$coef == level[1L] & rows$group == level[2L]]
      prior <- prior[nzchar(prior)]
      if (length(prior)) return(prior[1L])
    }
    ""
  }
  row <- function(parameter, task, prior, read = .read_prior(prior)) {
    list(parameter = parameter, task = task, family = read$family %||% "",
         args = read$args %||% numeric(0), prior = prior)
  }

  # Each mean needs a proper prior: brms leaves population-level effects flat
  # unless told otherwise, and a flat prior has no density.
  means <- lapply(names(terms), function(task) {
    prior <- effective("b", terms[[task]], "")
    out <- row("dprime", task, prior)
    if (!.known_prior(out)) {
      .usdt_stop("Bayes factors need a proper prior on `b_", terms[[task]],
                 "` from a family uSDT knows, and it has ",
                 if (nzchar(prior)) paste0("\"", prior, "\"")
                 else "the flat default", ".\n  Set one, e.g. ",
                 "prior(normal(0, 1), class = b).")
    }
    out
  })

  # The SDs enter only the intercept and slope priors, which read them when
  # they are needed, so an SD prior uSDT cannot simulate stops only those.
  sds <- lapply(names(terms), function(task) {
    row("sd_dprime", task, effective("sd", terms[[task]], group))
  })

  # A fitted model stores the LKJ prior on the Cholesky factor (class L), a
  # specified one on the correlations (class cor). One correlation of a K x K
  # LKJ(eta) block has the marginal of a 2 x 2 LKJ(eta - 1 + K/2), and the
  # block holds every term of the group that correlates with the direct
  # sensitivity.
  prior <- effective(c("cor", "L"), "", group)
  read <- .read_prior(prior)
  if (is.null(read) || !read$family %in% c("lkj", "lkj_corr_cholesky") ||
      !.known_prior(list(family = "lkj", args = read$args))) {
    .usdt_stop("Bayes factors need an LKJ prior on the correlations of `",
               group, "`, which have \"", prior, "\".")
  }
  sd_names <- variables[startsWith(variables, paste0("sd_", group, "__"))]
  members <- substring(sd_names, nchar(group) + 6L)
  pairs <- c(paste0("cor_", group, "__", direct, "__", members),
             paste0("cor_", group, "__", members, "__", direct))
  K <- 1 + sum(pairs %in% variables)
  cor <- row("cor_dprime", "both", prior,
             list(family = "lkj", args = read$args - 1 + K / 2))

  .priors_frame(c(means, sds, list(cor)))
}

# This function fits a density to posterior draws, within the bounds of the
# quantity: logspline, then oldlogspline, then, for an unbounded quantity,
# logspline on an asinh scale. A fit counts only if its distribution function
# stays within 0.05 of that of the draws at every percentile: on very
# heavy-tailed draws either logspline algorithm can report success and still
# return a density far from them. Sound fits stay within 0.01.
.bf_posterior <- function(x, bounds = NULL) {
  at <- stats::quantile(x, seq(0.01, 0.99, by = 0.01), names = FALSE)
  empirical <- stats::ecdf(x)(at)
  fits <- list(
    function() {
      fit <- .quiet_logspline(logspline::logspline, x, bounds)
      if (!is.null(fit)) {
        list(log_density = function(q) logspline::dlogspline(q, fit, log = TRUE),
             cdf = function(q) logspline::plogspline(q, fit))
      }
    },
    function() {
      fit <- .quiet_logspline(logspline::oldlogspline, x, bounds)
      if (!is.null(fit)) {
        list(log_density = function(q) log(logspline::doldlogspline(q, fit)),
             cdf = function(q) logspline::poldlogspline(q, fit))
      }
    },
    function() if (is.null(bounds)) .bf_posterior_asinh(x))
  for (fit in fits) {
    out <- fit()
    if (!is.null(out) && max(abs(out$cdf(at) - empirical)) < 0.05) return(out)
  }
  .usdt_stop("the posterior density could not be estimated from the draws.")
}

# This function runs one logspline algorithm without its messages and
# returns NULL when it fails.
.quiet_logspline <- function(f, x, bounds = NULL) {
  args <- c(list(x), if (!is.null(bounds)) {
    list(lbound = bounds[1L], ubound = bounds[2L])
  })
  out <- NULL
  utils::capture.output(out <- suppressWarnings(
    tryCatch(do.call(f, args), error = function(e) NULL)))
  out
}

# This function fits logspline to y = asinh((x - median) / MAD) and returns
# the density of x through the Jacobian 1 / (MAD * sqrt(1 + z^2)). The
# transform is close to the identity near the centre and logarithmic in the
# tails, so the very heavy tails that the intercept and slope have when
# sigma_D is near 0, which defeat both logspline algorithms, reach logspline
# as light ones. Against bridge sampling it is less precise than logspline
# when logspline converges, so it is only the last resort.
.bf_posterior_asinh <- function(x) {
  centre <- stats::median(x)
  scale <- stats::mad(x)
  z <- function(q) (q - centre) / scale
  fit <- .quiet_logspline(logspline::logspline, asinh(z(x)))
  if (is.null(fit)) return(NULL)
  list(log_density = function(q) {
         logspline::dlogspline(asinh(z(q)), fit, log = TRUE) - log(scale) -
           0.5 * log1p(z(q)^2)
       },
       cdf = function(q) logspline::plogspline(asinh(z(q)), fit))
}

# This function writes the part of the caption of the prior and posterior
# plot that depends on the panels: which priors are induced by the others,
# what the difference is, and which comparisons with zero are tested on the
# quantity with a prior of its own.
.bf_caption <- function(curves) {
  field <- function(name) vapply(curves, `[[`, curves[[1L]][[name]], name)
  labels <- c(diff = "difference", rho = "correlation", slope = "slope",
              intercept = "intercept")
  tested <- field("tested")
  induced <- unname(labels[unique(tested[field("induced")])])
  swaps <- unique(paste(field("quantity"), tested)[field("swapped")])
  paste0(c(
    if (length(induced)) {
      sprintf(" The %s prior%s induced by the other priors.",
              paste(induced, collapse = " and "),
              if (length(induced) > 1L) "s are" else " is")
    },
    if ("diff" %in% field("quantity")) {
      " The difference is the indirect minus the direct mean sensitivity."
    },
    vapply(strsplit(swaps, " ", fixed = TRUE), function(s) {
      sprintf(paste(" A %s compared with zero is tested on %s, which has its",
                    "sign and a prior of its own, so that panel shows %s."),
              s[1L], s[2L], s[2L])
    }, "")), collapse = "")
}

# This function returns the printed name of a quantity. The difference is
# named as in summary(), Delta d' (indirect minus direct).
.bf_name <- function(quantity) {
  if (quantity == "diff") paste0(.usdt_chars()$delta, "d'") else quantity
}

# This function writes a hypothesis with the printed name of its quantity.
.bf_label <- function(name, h) {
  if (h$op %in% c("in", "out")) {
    sprintf("%s %s [%s, %s]", name, h$op, format(h$value[1L]),
            format(h$value[2L]))
  } else {
    sprintf("%s %s %s", name, h$op, format(h$value))
  }
}

# This function writes the panel title as a plotmath expression, so the Greek
# letters and the set symbols are drawn on every graphics device.
.bf_math <- function(h, short) {
  symbol <- switch(h$summary, diff = "Delta*d*\"'\"", h$summary)
  statement <- switch(h$op,
    "=" = sprintf("%s == \"%s\"", symbol, format(h$value)),
    "<" = , ">" = sprintf("%s %s \"%s\"", symbol, h$op, format(h$value)),
    sprintf("%s %s \"[%s, %s]\"", symbol,
            if (h$op == "in") "%in%" else "%notin%",
            format(h$value[1L]), format(h$value[2L])))
  sprintf("bold(%s ~~ \"(%s)\")", statement, short)
}

# This function evaluates one hypothesis on the draws of the tested quantity,
# and summarises the draws of the quantity it names.
.bf_one <- function(h, draws, summary, prior, level) {
  post <- .bf_posterior(draws, prior$bounds)
  a <- (1 - level) / 2
  limits <- stats::quantile(summary, c(a, 1 - a), names = FALSE)
  name <- .bf_name(h$summary)

  # The printed card names the test in full, the panel strip in short.
  test <- switch(h$op, "=" = "Savage-Dickey density ratio test",
                 "<" = , ">" = "Directional test", "Interval test")
  short <- switch(h$op, "=" = "Savage-Dickey", "<" = , ">" = "directional",
                  "interval")
  if (h$swapped) {
    test <- sprintf("%s, on %s %s 0", test, h$quantity, h$op)
    short <- paste0(short, ", on ", h$quantity)
  }

  statements <- .hypothesis_statements(name, h)
  if (h$op == "=") {
    log_bf10 <- prior$log_density(h$value) - post$log_density(h$value)
    prob <- NA_real_
  } else {
    region <- function(cdf) switch(h$op,
      ">" = 1 - cdf(h$value),
      "<" = cdf(h$value),
      "in" = cdf(h$value[2L]) - cdf(h$value[1L]),
      out = 1 - (cdf(h$value[2L]) - cdf(h$value[1L])))
    # The posterior probability is the share of draws in H1: unlike the
    # logspline tail, which overstated it 3 SDs out, it is unbiased.
    prob <- region(stats::ecdf(draws))
    log_bf10 <- stats::qlogis(prob) - stats::qlogis(region(prior$cdf))
  }

  evidence <- as.character(cut(abs(log_bf10), log(c(0, 3, 10, 30, 100, Inf)),
                               c("anecdotal", "moderate", "strong",
                                 "very strong", "extreme"),
                               right = FALSE, include.lowest = TRUE))
  row <- data.frame(
    hypothesis = .bf_label(name, h), test = test,
    H1 = statements[["H1"]], H0 = statements[["H0"]],
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
  curve <- list(hypothesis = .bf_math(h, short), quantity = h$summary,
                tested = h$quantity, swapped = h$swapped,
                induced = prior$induced, op = h$op, value = h$value, grid = grid,
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
#' [usdt_bf()] draws this figure by default. Call `plot()` to draw it again,
#' or to draw some of its rows, e.g. `plot(b[3, ])`.
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
  # A row subset keeps every curve, so its row names pick the curves it shows.
  curves <- attr(x, "curves")
  rows <- suppressWarnings(as.integer(rownames(x)))
  if (is.null(curves) || !nrow(x) || anyNA(rows) ||
      any(rows < 1L | rows > length(curves))) {
    .usdt_stop("`x` has lost its prior and posterior curves; plot the ",
               "object returned by usdt_bf(), or a subset of its rows.")
  }
  curves <- curves[rows]
  panels <- vapply(curves, `[[`, "", "hypothesis")
  panel <- function(i, n) factor(rep(panels[i], n), levels = unique(panels))

  # Up to three hypotheses share one row; more fill a near-square grid.
  columns <- if (length(curves) <= 3L) length(curves) else
    ceiling(sqrt(length(curves)))

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

  # The posterior is a light purple and the prior the grey of what does not
  # come from the fit; the dashed prior line also tells them apart.
  colours <- c(Posterior = "#9B7FD4", Prior = "#A8AEB3")
  plot <- ggplot2::ggplot(densities)
  if (!is.null(areas)) {
    plot <- plot + ggplot2::geom_ribbon(
      data = areas,
      ggplot2::aes(x = .data[["x"]], ymin = 0, ymax = .data[["density"]],
                   group = .data[["piece"]]),
      inherit.aes = FALSE, fill = colours[["Posterior"]], alpha = 0.2
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
    ) +
      ggplot2::scale_fill_manual(values = colours, name = NULL)
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
                        ncol = columns, labeller = ggplot2::label_parsed) +
    ggplot2::scale_colour_manual(values = colours, name = NULL) +
    ggplot2::scale_linetype_manual(
      values = c(Posterior = "solid", Prior = "dashed"), name = NULL
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.03)) +
    # A third line of text, the posterior probability, needs more headroom.
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, if (any(!is.na(
        vapply(curves, `[[`, 0, "prob")))) 0.4 else 0.3))
    ) +
    ggplot2::labs(
      x = NULL, y = NULL,
      caption = .wrap_caption(
        "Posterior (solid) and prior (dashed) densities of each tested ",
        "quantity. For a point hypothesis the Bayes factor is the ratio of ",
        "the two marked heights at the dotted line (Savage-Dickey); for a ",
        "directional or interval hypothesis the shaded area is the posterior ",
        "probability of H1.", .bf_caption(curves),
        width = max(88L, 50L * columns)
      )
    ) +
    .panel_theme() +
    # The caption spans every column, and wider keys keep the dashed prior
    # recognisable in the legend.
    ggplot2::theme(legend.key.width = grid::unit(30, "pt"),
                   plot.caption.position = "plot")
}
