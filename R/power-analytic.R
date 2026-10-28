# power-analytic.R
# Approximate power of the uSDT hypotheses from the sampling distribution of
# the estimates
# Author: Ricardo Rey-Sáez
# Last modified: 28-10-2026

# Internal functions

# Both engines treat each subject's d' as its true d' plus a measurement error
# whose variance depends on the trials, independent between tasks. The observed
# d' then have covariance S = Sigma + diag(vD, vI), and the sample mean and
# covariance of N subjects are sufficient for the five values the hypotheses
# test. The closed form evaluates the tests at the population values; the
# moment engine draws the sample mean and covariance and decides each test as
# the fitted model would.

# Number of draws of the sample mean and covariance, and the seed that gives
# every design the same underlying draws, so that power changes smoothly with
# the design.
.moment_draws <- 2e4L
.moment_seed <- 20261028L

# This function runs `code` with a fixed seed and restores the random number
# stream of the session afterwards, so the approximate power is reproducible
# without changing the user's simulations.
.with_seed <- function(seed, code) {
  env <- globalenv()
  had_seed <- exists(".Random.seed", envir = env, inherits = FALSE)
  if (had_seed) old <- get(".Random.seed", envir = env, inherits = FALSE)
  on.exit(if (had_seed) assign(".Random.seed", old, envir = env) else
    rm(".Random.seed", envir = env))
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  code
}

# This function gives the nodes and weights of Gauss-Hermite quadrature for an
# expectation over a standard normal variable (Golub-Welsch): the weights sum
# to one.
.gauss_hermite <- function(n = 20L) {
  i <- seq_len(n - 1L)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- J[cbind(i + 1L, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  list(z = sqrt(2) * e$values, w = e$vectors[1L, ]^2)
}

# This function gives the measurement error variance of one subject's d' in
# each task, averaged over the subjects of the population. Each subject's
# variance is that of usdt_reliability() (.sensitivity_variance()) at its d'
# and criterion; the average is taken on the information scale (a harmonic
# mean), which weighs subjects as the likelihood does. `pop` is one scenario,
# `cells` the signal and noise trials of each task, and `split` the tasks
# whose criterion the median split fixes at zero.
.error_variance <- function(pop, cells, split) {
  gh <- .gauss_hermite()
  out <- vapply(c("direct", "indirect"), function(k) {
    s <- if (k == "direct") "D" else "I"
    d <- pop[[paste0("mu_", s)]] + pop[[paste0("sigma_", s)]] * gh$z
    sd_c <- pop[[paste0("sd_criterion_", s)]]
    if (split[[k]]) {
      criterion <- 0
      w_c <- 1
    } else if (sd_c > 0) {
      criterion <- pop[[paste0("criterion_", s)]] + sd_c * gh$z
      w_c <- gh$w
    } else {
      criterion <- pop[[paste0("criterion_", s)]]
      w_c <- 1
    }

    # Every pair of d' and criterion nodes is one kind of subject. The model
    # writes the criterion with the opposite sign, as the intercept of the
    # task.
    dd <- rep(d, times = length(criterion))
    cc <- rep(criterion, each = length(d))
    weight <- rep(gh$w, times = length(criterion)) * rep(w_c, each = length(d))
    eta <- cbind(-cc - dd / 2, -cc + dd / 2)
    task_cells <- data.frame(sig = c(FALSE, TRUE),
                             n = c(cells$noise[[k]], cells$signal[[k]]))
    v <- .sensitivity_variance(task_cells, eta, free = !split[[k]])

    # A subject at an extreme rate carries no information about its d'.
    information <- 1 / v
    information[is.na(information)] <- 0
    1 / sum(weight * information)
  }, 0)
  if (any(!is.finite(out))) {
    .usdt_stop("the trials carry no information about the d' of the ",
               names(out)[!is.finite(out)][1L], " task: its rates are at 0 ",
               "or 1 for every subject. Check the means and SDs of d' and of ",
               "the criterion.")
  }
  out
}

# This function gives the variance of g' x, where x holds the estimates of
# the five values (gamma_D, gamma_I, s2_D, s2_I, s_DI) of N subjects whose
# observed d' have covariance S (entries S11, S22, S12). The means have
# covariance S / N and the covariance entries Cov(s_ij, s_kl) =
# (S_ik S_jl + S_il S_jk) / N, independent of the means. `g` has one row per
# set of values, and S one entry per row or a single one.
.moment_quad <- function(g, S11, S22, S12, N) {
  g <- matrix(g, ncol = 5L)
  means <- g[, 1L]^2 * S11 + 2 * g[, 1L] * g[, 2L] * S12 + g[, 2L]^2 * S22
  covariances <- 2 * g[, 3L]^2 * S11^2 + 2 * g[, 4L]^2 * S22^2 +
    g[, 5L]^2 * (S11 * S22 + S12^2) + 4 * g[, 3L] * g[, 4L] * S12^2 +
    4 * g[, 3L] * g[, 5L] * S11 * S12 + 4 * g[, 4L] * g[, 5L] * S22 * S12
  (means + covariances) / N
}

# This function gives the estimate and standard error of each tested quantity
# from the five values `P` (a matrix with columns gamma_D, gamma_I, s2_D, s2_I
# and s_DI) and the observed covariance S, with the gradients the fitted model
# uses (.diff_rows(), .cor_rows(), .reg_rows() and .cov_test()). The
# correlation also has its Fisher-z estimate and standard error.
.moment_tests <- function(P, S11, S22, S12, N) {
  gD <- P[, "gamma_D"]; gI <- P[, "gamma_I"]
  s2D <- P[, "s2_D"]; s2I <- P[, "s2_I"]; sDI <- P[, "s_DI"]
  n <- nrow(P)
  zero <- rep(0, n)

  # A draw on the boundary has no correlation, and its Fisher-z needs one
  # inside (-1, 1); those draws are not estimable anyway.
  positive <- s2D > 0 & s2I > 0
  root <- sqrt(ifelse(positive, s2D * s2I, NA_real_))
  rho <- sDI / root
  inside <- !is.na(rho) & abs(rho) < 1
  fisher <- rep(NA_real_, n)
  fisher[inside] <- atanh(rho[inside])
  slope <- sDI / s2D

  se <- function(g) sqrt(pmax(.moment_quad(g, S11, S22, S12, N), 0))
  se_rho <- se(cbind(zero, zero, -rho / (2 * s2D), -rho / (2 * s2I), 1 / root))
  list(
    diff = list(estimate = gI - gD, se = se(cbind(-1, 1, zero, zero, zero))),
    covariance = list(estimate = sDI, se = se(cbind(zero, zero, zero, zero, 1))),
    rho = list(estimate = rho, se = se_rho),
    fisher = list(estimate = fisher, se = se_rho / (1 - rho^2)),
    intercept = list(estimate = gI - slope * gD,
                     se = se(cbind(-slope, 1, gD * sDI / s2D^2, zero, -gD / s2D))),
    slope = list(estimate = slope,
                 se = se(cbind(zero, zero, -sDI / s2D^2, zero, 1 / s2D))))
}

# This function tells, for each set of five values, whether the fitted model
# would report the correlation, slope and intercept: hsdt() needs the
# sensitivity covariance inside its valid range and every variance parameter
# away from its bound (.usdt_pars() and .full_vcov()).
.moment_estimable <- function(P) {
  s2D <- P[, "s2_D"]; s2I <- P[, "s2_I"]; sDI <- P[, "s_DI"]
  inside <- s2D > 1e-8 & s2I > 1e-8
  rho <- ifelse(inside, sDI / sqrt(pmax(s2D, 1e-300) * pmax(s2I, 1e-300)), 1)
  inside & abs(rho) < 1 - 1e-6 &
    sqrt(pmax(s2D, 0)) >= 1e-3 & sqrt(pmax(s2I * (1 - rho^2), 0)) >= 1e-3
}

# This function names the test the analysis runs for one hypothesis from
# .parse_hypothesis(). usdt_tests(), which hsdt() reports, tests a null
# correlation or slope on the covariance they share; every other statement is
# tested as in usdt_equivalence(), the correlation on the Fisher-z scale.
.power_test <- function(h) {
  if (h$op == "=" && h$value == 0 && h$quantity %in% c("rho", "slope")) {
    return("covariance")
  }
  if (h$quantity == "rho") "fisher" else h$quantity
}

# This function gives the bounds of a hypothesis on the scale of its test.
.power_bounds <- function(h) {
  if (.power_test(h) == "fisher") atanh(h$value) else h$value
}

# This function decides each test from its estimates and standard errors, as
# usdt_tests() and usdt_equivalence() do at level alpha.
.power_reject <- function(op, bounds, estimate, se, alpha) {
  z1 <- stats::qnorm(1 - alpha)
  z2 <- stats::qnorm(1 - alpha / 2)
  t <- function(b) (estimate - b) / se
  out <- switch(op,
    "=" = abs(t(bounds)) > z2,
    ">" = t(bounds) > z1,
    "<" = t(bounds) < -z1,
    "in" = t(bounds[1L]) > z1 & t(bounds[2L]) < -z1,
    out = t(bounds[1L]) < -z2 | t(bounds[2L]) > z2)
  out & is.finite(se) & se > 0
}

# This function gives the power of a test whose estimate is normal around
# `theta` with standard error `se`, for each kind of statement.
.power_normal <- function(op, bounds, theta, se, alpha) {
  z1 <- stats::qnorm(1 - alpha)
  z2 <- stats::qnorm(1 - alpha / 2)
  a <- bounds[1L]
  b <- bounds[length(bounds)]
  switch(op,
    "=" = stats::pnorm((theta - a) / se - z2) + stats::pnorm(-(theta - a) / se - z2),
    ">" = stats::pnorm((theta - a) / se - z1),
    "<" = stats::pnorm((a - theta) / se - z1),
    "in" = pmax(0, stats::pnorm((b - theta) / se - z1) -
                  stats::pnorm((a - theta) / se + z1)),
    out = stats::pnorm((a - theta) / se - z2) + stats::pnorm((theta - b) / se - z2))
}

# This function tells whether the population value of a quantity lies in the
# null hypothesis of a statement, where the rejection rate is the type I error
# rate.
.power_in_null <- function(op, value, theta, tol = 1e-10) {
  a <- value[1L]
  b <- value[length(value)]
  switch(op,
    "=" = abs(theta - a) < tol,
    ">" = theta <= a + tol,
    "<" = theta >= a - tol,
    "in" = theta <= a + tol | theta >= b - tol,
    out = theta >= a - tol & theta <= b + tol)
}

# This function gives the standard normal and uniform draws that the moment
# engine turns into sample means and covariances. The uniforms become the
# chi-square variables of the Bartlett decomposition through their quantile
# function, so a change in N moves every draw smoothly.
.moment_base <- function(R = .moment_draws) {
  .with_seed(.moment_seed, list(
    z1 = stats::rnorm(R), z2 = stats::rnorm(R), n21 = stats::rnorm(R),
    u1 = stats::runif(R), u2 = stats::runif(R)))
}

# This function gives the two chi variables of the Bartlett decomposition for
# N subjects. They depend on N only, so a grid of designs computes them once
# for each number of subjects.
.moment_chi <- function(base, N) {
  list(c1 = sqrt(stats::qchisq(base$u1, N - 1)),
       c2 = sqrt(stats::qchisq(base$u2, N - 2)))
}

# This function evaluates every hypothesis in `tests` (from .parse_hypothesis())
# for one scenario `pop`, N subjects and the error variances `v`. It returns,
# for each hypothesis, the population value of its quantity, the closed-form
# standard error and power, and the share of draws of the moment engine that
# the analysis could estimate and that reject H0; and the share that rejects
# every hypothesis at once.
.power_approximate <- function(pop, N, v, tests, alpha, base,
                               chi = .moment_chi(base, N)) {
  s2D <- pop$sigma_D^2; s2I <- pop$sigma_I^2
  sDI <- pop$rho * pop$sigma_D * pop$sigma_I
  S11 <- s2D + v[["direct"]]; S22 <- s2I + v[["indirect"]]; S12 <- sDI

  # The closed form evaluates the tests at the population values.
  truth <- cbind(gamma_D = pop$mu_D, gamma_I = pop$mu_I, s2_D = s2D,
                 s2_I = s2I, s_DI = sDI)
  expected <- .moment_tests(truth, S11, S22, S12, N)

  # The sample mean and covariance of each draw (Bartlett decomposition of
  # the Wishart distribution with N - 1 degrees of freedom).
  df <- N - 1
  l11 <- sqrt(S11); l21 <- S12 / l11; l22 <- sqrt(S22 - l21^2)
  a11 <- chi$c1^2; a12 <- chi$c1 * base$n21; a22 <- base$n21^2 + chi$c2^2
  H11 <- l11^2 * a11 / df
  H12 <- l11 * (l21 * a11 + l22 * a12) / df
  H22 <- (l21^2 * a11 + 2 * l21 * l22 * a12 + l22^2 * a22) / df
  P <- cbind(gamma_D = pop$mu_D + l11 * base$z1 / sqrt(N),
             gamma_I = pop$mu_I + (l21 * base$z1 + l22 * base$z2) / sqrt(N),
             s2_D = H11 - v[["direct"]], s2_I = H22 - v[["indirect"]],
             s_DI = H12)
  drawn <- .moment_tests(P, H11, H22, H12, N)
  estimable <- .moment_estimable(P)

  values <- c(diff = pop$diff, rho = pop$rho, intercept = pop$intercept,
              slope = pop$slope)
  decisions <- matrix(FALSE, length(base$z1), length(tests))
  rows <- vector("list", length(tests))
  for (j in seq_along(tests)) {
    h <- tests[[j]]
    test <- .power_test(h)
    bounds <- .power_bounds(h)
    theta <- values[[h$quantity]]

    # The difference is a fixed effect, which the fitted model always reports.
    ok <- if (h$quantity == "diff") rep(TRUE, length(estimable)) else estimable
    reject <- ok & .power_reject(h$op, bounds, drawn[[test]]$estimate,
                                 drawn[[test]]$se, alpha)
    reject[is.na(reject)] <- FALSE
    decisions[, j] <- reject

    # The reported standard error is that of the quantity on its own scale.
    centre <- switch(test, covariance = sDI, fisher = atanh(pop$rho), theta)
    named <- switch(test, covariance = h$quantity, fisher = "rho", test)
    rows[[j]] <- data.frame(
      true = theta, se = expected[[named]]$se,
      power_closed = .power_normal(h$op, bounds, centre, expected[[test]]$se,
                                   alpha),
      power = mean(reject), estimable = mean(ok),
      in_null = .power_in_null(h$op, h$value, theta))
  }
  needs_covariance <- any(vapply(tests, function(h) h$quantity != "diff", TRUE))
  list(rows = do.call(rbind, rows),
       joint = mean(rowSums(decisions) == length(tests)),
       joint_estimable = if (needs_covariance) mean(estimable) else 1)
}
