# equivalence.R
# Equivalence, one-sided and minimal-effect tests for the three hypotheses of
# a frequentist hierarchical SDT model
# Author: Ricardo Rey-Sáez
# Last modified: 27-10-2026

# Public functions

#' Equivalence, non-inferiority and minimal-effect tests
#'
#' Tests interval, one-sided and point hypotheses about the quantities of the
#' three uSDT hypotheses in a frequentist [hsdt()] fit or a custom model fitted
#' with [lme4::glmer()]: the difference between mean sensitivities (`diff`, H1,
#' reported as \eqn{\Delta d'}, indirect minus direct), their correlation
#' (`rho`, H2) and the slope and intercept of the latent regression (`slope`,
#' `intercept`, H3). The hypotheses are written as in [usdt_bf()], so a
#' frequentist and a Bayesian analysis can state the same bounds.
#'
#' @param fit A frequentist `hsdt` object, with or without the replicates of
#'   [usdt_boot()], or a `glmerMod` from [lme4::glmer()].
#' @param hypothesis Character vector of hypotheses. It has no default: the
#'   bounds are the smallest effects of interest, which only the analysis can
#'   set. Each one is `"<quantity> <op> <value>"`, where the quantity is
#'   `diff`, `rho`, `slope` or `intercept` and the operator is `=`, `<` or `>`
#'   (`<=` and `>=` are read as `<` and `>`), or `"<quantity> in [a, b]"` and
#'   `"<quantity> out [a, b]"` for a region and its complement.
#' @param direct,indirect For a `glmerMod`, the fixed effects of the direct
#'   and indirect sensitivities, as in [usdt_tests()]. Defaults match the
#'   internal naming of [hsdt()].
#' @param alpha Significance level of every test, between 0 and 0.5 (default
#'   0.05).
#'
#' @return An object of class `usdt_equivalence`: a data frame with one row per
#'   hypothesis and the columns `hypothesis`, `test` (`"Two-sided test"`,
#'   `"One-sided test"`, `"Equivalence test (TOST)"` or `"Minimal-effect
#'   test"`), `H1`, `H0`, `estimate`, `se`, `conf.level`, `conf.low` and
#'   `conf.high` (the interval that decides the test), `ci_method`, `p.lower`
#'   and `p.upper` (the one-sided p-values at the two bounds of a region),
#'   `p.value`, `decision`, `status` and `reason`. Rows that cannot be
#'   estimated have `status = "not estimable"` and give the `reason`. It prints
#'   one report per hypothesis.
#'
#' @details
#' # The tests
#'
#' A point statement names the null hypothesis H0; an inequality or a region
#' names the alternative H1, as in [usdt_bf()]. Each test rejects H0 at level
#' `alpha`:
#'
#' * `q = c`: the two-sided test of H0: \eqn{q = c}.
#' * `q > c` and `q < c`: one-sided tests. With a margin on the unfavourable
#'   side of the reference value, `diff > -0.2` is a non-inferiority test and
#'   `diff < 0.2` a non-superiority test.
#' * `q in [a, b]`: an equivalence test by two one-sided tests (TOST;
#'   Schuirmann, 1987; Lakens et al., 2018). H0, \eqn{q \le a} or
#'   \eqn{q \ge b}, is rejected when both one-sided tests reject their part of
#'   it, so `p.value` is the larger of `p.lower` and `p.upper`.
#' * `q out [a, b]`: a minimal-effect test. H0, \eqn{a \le q \le b}, is
#'   rejected when either one-sided test rejects it at `alpha / 2`, so
#'   `p.value` is twice the smaller of `p.lower` and `p.upper`. As the region
#'   shrinks to a point it becomes the two-sided test of that point.
#'
#' Each test is decided by the position of a confidence interval, and the row
#' reports that interval: of level `1 - 2 * alpha` for the one-sided and
#' equivalence tests, and `1 - alpha` for the two-sided and minimal-effect
#' tests, which split `alpha` between two tails. With the default
#' `alpha = 0.05`, an equivalence test reports a 90% interval, and H0 is
#' rejected exactly when that interval lies inside the region.
#'
#' # The correlation
#'
#' The bounds of `rho` are given, and its estimate and interval reported, on
#' the correlation scale. The test runs on the Fisher-\eqn{z} scale, where the
#' sampling distribution of the estimate is closer to normal:
#' \deqn{z = \frac{\mathrm{atanh}(\hat\rho) - \mathrm{atanh}(c)}{SE /
#' (1 - \hat\rho^2)}.}{z = (atanh(rho) - atanh(c)) / (SE / (1 - rho^2)).}
#' Because `atanh()` is increasing, \eqn{\rho > c} and
#' \eqn{\mathrm{atanh}(\rho) > \mathrm{atanh}(c)}{atanh(rho) > atanh(c)} are
#' the same statement, so the p-value and the decision are those of the
#' hypothesis about `rho`, and the interval, transformed back, stays inside
#' (-1, 1).
#'
#' Every p-value comes from the interval reported with it. [usdt_tests()]
#' instead tests a null association on the covariance that `rho` and the slope
#' share, so its p-values for `rho = 0` and `slope = 0` can differ slightly
#' from these.
#'
#' # The intercept and the slope
#'
#' Both divide by the variance of the direct sensitivity. When the direct
#' sensitivity barely varies across subjects, their standard errors are large
#' and their Wald intervals inaccurate, so their equivalence tests need many
#' subjects to succeed, and the bootstrap of [usdt_boot()] gives more accurate
#' intervals.
#'
#' # Boundary fits and the bootstrap
#'
#' The correlation, slope and intercept need the full covariance matrix of the
#' fit, which is unavailable when the sensitivity covariance is on the
#' boundary (a variance estimated as zero, or a correlation of -1 or 1).
#' Their rows are then marked `"not estimable"`, while the difference, which
#' only needs the fixed effects, is still tested. After [usdt_boot()], every
#' test uses the bootstrap replicates and the interval type chosen there. The
#' p-value at a bound is the share of the bootstrap distribution beyond it,
#' counted as \eqn{(k + 1) / (B + 1)} as in [usdt_boot()], so the decision
#' agrees with the reported interval up to that correction.
#'
#' @references
#' Lakens, D., Scheel, A. M., & Isager, P. M. (2018). Equivalence testing for
#' psychological research: A tutorial. *Advances in Methods and Practices in
#' Psychological Science*, 1(2), 259-269. \doi{10.1177/2515245918770963}
#'
#' Schuirmann, D. J. (1987). A comparison of the two one-sided tests procedure
#' and the power approach for assessing the equivalence of average
#' bioavailability. *Journal of Pharmacokinetics and Biopharmaceutics*, 15(6),
#' 657-680. \doi{10.1007/BF01068419}
#'
#' @seealso [usdt_tests()], [usdt_boot()], [usdt_bf()]
#'
#' @examples
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
#' # Mean sensitivities equivalent within 0.2 d' units, and an unconscious
#' # effect at zero awareness smaller than 0.2
#' usdt_equivalence(m, c("diff in [-0.2, 0.2]", "intercept < 0.2"))
#'
#' # The bounds of the correlation are correlations
#' usdt_equivalence(m, "rho in [-0.3, 0.3]")
#'
#' @export
usdt_equivalence <- function(fit, hypothesis, direct = "d_D",
                             indirect = "d_I", alpha = 0.05) {
  if (.is_posterior(fit)) {
    .usdt_stop("usdt_equivalence() tests maximum-likelihood fits. Test a ",
               "Bayesian fit with usdt_bf().")
  }
  if (!inherits(fit, c("hsdt", "glmerMod"))) {
    .usdt_stop("`fit` must come from hsdt() or lme4::glmer(), not a plain ",
               class(fit)[1L], ".")
  }
  if (missing(hypothesis)) {
    .usdt_stop("`hypothesis` has no default: state the bounds to test, e.g. ",
               "\"diff in [-0.2, 0.2]\" or \"intercept < 0.2\".")
  }
  if (!is.character(hypothesis) || !length(hypothesis) || anyNA(hypothesis)) {
    .usdt_stop("`hypothesis` must be one or more strings such as ",
               "\"diff in [-0.2, 0.2]\".")
  }
  .check_scalar_number(alpha, "alpha", lower = 0, upper = 0.5,
                       open_lower = TRUE, open_upper = TRUE)
  tests <- lapply(hypothesis, .parse_hypothesis)

  # The replicates of a bootstrapped fit stay usable where the Wald
  # approximation fails, so they test every quantity once they exist.
  boot <- if (inherits(fit, "hsdt")) fit$boot
  sampling <- if (isTRUE(boot$summary_available)) {
    .equiv_boot(boot)
  } else {
    .equiv_wald(.pars_of(fit, direct, indirect),
                bootable = inherits(fit, "hsdt"))
  }

  out <- do.call(rbind, lapply(tests, function(h) {
    .equiv_one(h, sampling[[h$quantity]], alpha)
  }))
  rownames(out) <- NULL
  structure(out, alpha = alpha, class = c("usdt_equivalence", "data.frame"))
}

#' @param x A `usdt_equivalence` object.
#' @param digits Number of decimals printed.
#' @param width Width of the printed report, in characters.
#' @param ... Ignored.
#' @rdname usdt_equivalence
#' @export
print.usdt_equivalence <- function(x, digits = 3L, width = 80L, ...) {
  # A subset that lost columns is printed as the data frame it now is.
  needed <- c("test", "H1", "H0", "estimate", "se", "conf.level", "conf.low",
              "conf.high", "ci_method", "p.lower", "p.upper", "p.value",
              "decision", "status", "reason")
  if (!all(needed %in% names(x))) return(NextMethod())

  width <- max(60L, as.integer(width))
  chars <- .usdt_chars()
  alpha <- attr(x, "alpha")
  number <- function(v) {
    if (is.na(v)) return("NA")
    if (is.infinite(v)) return(if (v > 0) "Inf" else "-Inf")
    formatC(v, format = "f", digits = digits)
  }
  p <- function(v) trimws(.fmt_p(v))
  # Long values, such as the reason a row is not estimable, wrap under their
  # own column.
  row <- function(label, value) {
    lines <- strwrap(value, width = width - 32L)
    cat(sprintf("  %-27s : %s\n", label, lines[1L]))
    if (length(lines) > 1L) {
      cat(paste0(strrep(" ", 32L), lines[-1L], "\n"), sep = "")
    }
  }

  # Every hypothesis gets a complete card, so mixed test types stay explicit.
  for (i in seq_len(nrow(x))) {
    if (i > 1L) cat("\n")
    cat(" uSDT hypothesis ", chars$h, " ", x$test[[i]], "\n", sep = "")
    cat(strrep(chars$hh, width), "\n", sep = "")
    row("Hypothesis (H1)", x$H1[[i]])
    row("Hypothesis (H0)", x$H0[[i]])
    row("Estimate", number(x$estimate[[i]]))
    if (x$status[[i]] != "ok") {
      row("Not estimable", x$reason[[i]])
    } else {
      row("Standard error", number(x$se[[i]]))
      row(sprintf("%g%% CI, %s", 100 * x$conf.level[[i]], x$ci_method[[i]]),
          sprintf("[%s, %s]", number(x$conf.low[[i]]),
                  number(x$conf.high[[i]])))
      if (!is.na(x$p.lower[[i]])) {
        row("p-value, lower bound", p(x$p.lower[[i]]))
        row("p-value, upper bound", p(x$p.upper[[i]]))
      }
      row("p-value", p(x$p.value[[i]]))
      row(sprintf("Decision (alpha = %g)", alpha),
          if (x$decision[[i]] == "reject H0") "H0 rejected, in favour of H1"
          else "H0 not rejected")
    }
    cat(strrep(chars$h, width), "\n", sep = "")
  }
  invisible(x)
}

# Internal functions

# This function describes the Wald sampling distribution of each tested
# quantity from the fitted parameters `p`: normal around the estimate, and for
# the correlation normal on the Fisher-z scale, so that its bounds, estimate
# and interval stay on the correlation scale. `bootable` adds to the reason of
# an unavailable row that usdt_boot() can test it.
.equiv_wald <- function(p, bootable = FALSE) {
  rows <- rbind(.diff_rows(p, 0.95), .cor_rows(p, 0.95), .reg_rows(p, 0.95))
  quantities <- c("diff", "rho", "intercept", "slope")
  identity_link <- list(to = identity, from = identity)

  out <- lapply(seq_len(nrow(rows)), function(i) {
    r <- rows[i, ]
    rho <- quantities[i] == "rho"
    link <- if (rho) .fisher_link else identity_link
    se <- if (rho) .fisher_z(r$estimate, r$se)$se_z else r$se
    if (r$status != "ok" || !is.finite(se) || se <= 0) {
      reason <- if (is.na(r$reason)) "the standard error is invalid" else
        r$reason
      if (bootable) {
        reason <- paste0(reason, "; usdt_boot() can test it on bootstrap ",
                         "replicates")
      }
      return(list(ok = FALSE, estimate = r$estimate, reason = reason))
    }
    centre <- link$to(r$estimate)
    list(ok = TRUE, estimate = r$estimate, se = r$se,
         method = if (rho) "Fisher-z" else "Wald",
         cdf = function(x) stats::pnorm(link$to(x), centre, se),
         sf = function(x) stats::pnorm(link$to(x), centre, se,
                                       lower.tail = FALSE),
         ci = function(level) {
           link$from(centre + c(-1, 1) * stats::qnorm(1 - (1 - level) / 2) * se)
         })
  })
  stats::setNames(out, quantities)
}

# This function describes the bootstrap distribution of each tested quantity
# from the `boot` element of usdt_boot(): the distribution whose quantiles are
# the interval of the chosen type, so that each test agrees with the interval
# it reports. The basic and normal intervals of the correlation are built on
# the Fisher-z scale, and so are its tests.
.equiv_boot <- function(boot) {
  t <- boot$t[boot$ok, , drop = FALSE]
  type <- boot$type
  identity_link <- list(to = identity, from = identity)

  out <- lapply(c("diff", "rho", "intercept", "slope"), function(q) {
    v <- t[, q]
    est <- boot$t0[[q]]
    link <- if (q == "rho") .fisher_link else NULL
    to <- if (q == "rho" && type != "perc") link$to else identity_link$to
    if (!all(is.finite(to(c(est, v))))) {
      return(list(ok = FALSE, estimate = est, reason = paste0(
        "the correlation reaches -1 or 1, where the Fisher-z scale of the \"",
        type, "\" bootstrap is undefined; the percentile bootstrap (type = ",
        "\"perc\") can test it")))
    }

    # The percentile interval is read from the replicates themselves, the
    # basic one from their reflection about the estimate, and the normal one
    # from a normal distribution with their bias and spread.
    centre <- to(est)
    if (type == "norm") {
      m <- 2 * centre - mean(to(v))
      s <- stats::sd(to(v))
      cdf <- function(x) stats::pnorm(to(x), m, s)
      sf <- function(x) stats::pnorm(to(x), m, s, lower.tail = FALSE)
    } else {
      d <- if (type == "perc") v else 2 * centre - to(v)
      cdf <- function(x) (sum(d <= to(x)) + 1) / (length(d) + 1)
      sf <- function(x) (sum(d >= to(x)) + 1) / (length(d) + 1)
    }
    list(ok = TRUE, estimate = est, se = stats::sd(v),
         method = paste0("bootstrap (", type, ")"), cdf = cdf, sf = sf,
         ci = function(level) .boot_ci(v, est, level, type, link = link)[1L, ])
  })
  stats::setNames(out, c("diff", "rho", "intercept", "slope"))
}

# This function tests one hypothesis on the sampling distribution `s` of its
# quantity. `cdf` and `sf` give the probability below and above a value, so
# each one-sided p-value is the probability beyond its bound.
.equiv_one <- function(h, s, alpha) {
  name <- .bf_name(h$quantity)
  statements <- .hypothesis_statements(name, h)
  test <- switch(h$op, "=" = "Two-sided test", "<" = , ">" = "One-sided test",
                 "in" = "Equivalence test (TOST)", "Minimal-effect test")

  # The two-sided and minimal-effect tests split alpha between two tails; the
  # others put all of it in one tail at each bound.
  level <- if (h$op %in% c("=", "out")) 1 - alpha else 1 - 2 * alpha
  row <- data.frame(
    hypothesis = .bf_label(name, h), test = test,
    H1 = statements[["H1"]], H0 = statements[["H0"]],
    estimate = s$estimate, se = NA_real_, conf.level = level,
    conf.low = NA_real_, conf.high = NA_real_, ci_method = "not available",
    p.lower = NA_real_, p.upper = NA_real_, p.value = NA_real_,
    decision = NA_character_, status = "not estimable",
    reason = if (s$ok) NA_character_ else s$reason,
    stringsAsFactors = FALSE)
  if (!s$ok) return(row)

  a <- h$value[1L]
  b <- h$value[length(h$value)]
  p <- switch(h$op,
    "=" = c(NA, NA, min(1, 2 * min(s$cdf(a), s$sf(a)))),
    ">" = c(NA, NA, s$cdf(a)),
    "<" = c(NA, NA, s$sf(a)),
    "in" = c(s$cdf(a), s$sf(b), max(s$cdf(a), s$sf(b))),
    out = c(s$sf(a), s$cdf(b), min(1, 2 * min(s$sf(a), s$cdf(b)))))
  ci <- s$ci(level)

  row$se <- s$se
  row$conf.low <- ci[1L]
  row$conf.high <- ci[2L]
  row$ci_method <- s$method
  row$p.lower <- p[1L]
  row$p.upper <- p[2L]
  row$p.value <- p[3L]
  row$decision <- if (p[3L] < alpha) "reject H0" else "do not reject H0"
  row$status <- "ok"
  row
}
