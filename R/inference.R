# inference.R
# This script builds confidence intervals for the reported estimates.
# Author: Ricardo Rey-Sáez
# Last modified: 04-09-2026

# Shared calculations

# This function calculates approximate standard errors. `grad` holds one
# gradient, or one gradient per row when a whole curve is needed at once.
.delta_se <- function(grad, V) {
  g <- if (is.matrix(grad)) grad else matrix(grad, nrow = 1L)
  v <- rowSums((g %*% V) * g)
  v[!is.finite(v) | v < 0] <- NA_real_
  sqrt(v)
}

# This function returns a Wald test and interval.
.wald <- function(est, se, level = 0.95, null = 0) {
  if (!is.finite(est) || !is.finite(se) || se <= 0) {
    return(list(estimate = est, se = se, statistic = NA_real_,
                p.value = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
                ci_method = "not available"))
  }
  z  <- qnorm(1 - (1 - level) / 2)
  st <- (est - null) / se
  list(estimate = est, se = se, statistic = st,
       p.value  = 2 * stats::pnorm(-abs(st)),
       conf.low = est - z * se, conf.high = est + z * se,
       ci_method = "Wald")
}

# Fisher-z interval

# This function builds an interval on the Fisher-z scale.
.fisher_z <- function(rho, se_rho, level = 0.95) {

  # The transformation is undefined at the boundary.
  if (!is.finite(rho) || !is.finite(se_rho) || abs(rho) >= 1) {
    return(list(conf.low = NA_real_, conf.high = NA_real_, se_z = NA_real_,
                ci_method = "Fisher-z"))
  }

  # The interval returns to the correlation scale.
  z    <- qnorm(1 - (1 - level) / 2)
  zf   <- atanh(rho)
  se_z <- se_rho / (1 - rho^2)
  list(conf.low = tanh(zf - z * se_z), conf.high = tanh(zf + z * se_z),
       se_z = se_z, ci_method = "Fisher-z")
}

# This pair moves a correlation to an unbounded scale and back.
.fisher_link <- list(to = atanh, from = tanh)

# Bootstrap intervals

# This function builds one interval for each column of bootstrap replicates,
# so a single call covers one estimate or every point of a curve.
.boot_ci <- function(v, estimate, level, type, link = NULL) {
  v <- as.matrix(v)
  a <- (1 - level) / 2

  # A monotone link leaves the percentile interval unchanged.
  if (type == "perc") {
    ci <- t(apply(v, 2L, stats::quantile, probs = c(a, 1 - a), names = FALSE))
  } else {
    if (!is.null(link)) {
      v        <- link$to(v)
      estimate <- link$to(estimate)
    }
    ci <- switch(type,
      basic = 2 * estimate - t(apply(v, 2L, stats::quantile,
                                     probs = c(1 - a, a), names = FALSE)),
      norm  = 2 * estimate - colMeans(v) +
              outer(apply(v, 2L, stats::sd), c(-1, 1) * qnorm(1 - a)))
    if (!is.null(link)) ci <- link$from(ci)
  }

  # Boundary values have no finite centre on the transformed scale, and the
  # replicate names carry no meaning for an interval.
  ci[!is.finite(ci)] <- NA_real_
  dimnames(ci) <- NULL
  ci
}

# Result rows

# This function creates one result row.
.row <- function(term, estimate, se = NA_real_, statistic = NA_real_,
                 p.value = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
                 ci_method = NA_character_, status = "ok",
                 reason = NA_character_) {
  data.frame(term = term, estimate = estimate, se = se, statistic = statistic,
             p.value = p.value, conf.low = conf.low, conf.high = conf.high,
             ci_method = ci_method, status = status, reason = reason,
             stringsAsFactors = FALSE)
}

# This function marks an unavailable result and gives its reason.
.row_na <- function(term, estimate, reason) {
  .row(term, estimate, ci_method = "not available", status = "not estimable",
       reason = reason)
}
