# helper-trials.R
# This script builds trial-level fixtures for the tests. It is not part of the
# package: it only gives the tests a data set with a known structure.
# Author: Ricardo Rey-Sáez
# Last modified: 09-09-2026

# This function draws two correlated sensitivities per subject and returns one
# balanced block of trials for each task. With `rt = TRUE` the influence task
# returns lognormal times that fall as the evidence rises, so a faster response
# is the signal response and a median split has something to work on.
make_trials <- function(n_subj = 50, n_trials = 100,
                        gamma_A = 0.8, gamma_I = 0.3,
                        sd_A = 0.5, sd_I = 0.3, rho = 0.5,
                        crit_A = 0, sd_crit = 0.3,
                        crit_I = 0, sd_crit_I = 0.3,
                        rt = FALSE) {

  # The model draws two related sensitivities for each subject.
  S <- matrix(c(sd_A^2, rho * sd_A * sd_I,
                rho * sd_A * sd_I, sd_I^2), nrow = 2L)
  u  <- matrix(stats::rnorm(n_subj * 2L), n_subj) %*% chol(S)

  # A median split leaves the influence criterion at zero without variation.
  if (rt) { crit_I <- 0; sd_crit_I <- 0 }
  uc <- list(A = stats::rnorm(n_subj, 0, sd_crit),
             I = stats::rnorm(n_subj, 0, sd_crit_I))

  # Each subject receives a balanced block for both tasks.
  half <- n_trials %/% 2L
  cond <- rep(c(0L, 1L), each = half)
  dev  <- cond - 0.5

  out <- vector("list", 2L * n_subj)
  k <- 1L
  for (task in c("A", "I")) {
    g  <- if (task == "A") gamma_A else gamma_I
    uu <- if (task == "A") u[, 1L] else u[, 2L]
    cc <- if (task == "A") crit_A + uc$A else crit_I + uc$I

    for (i in seq_len(n_subj)) {

      # The model calculates the evidence for every trial.
      eta <- cc[i] + (g + uu[i]) * dev
      row <- data.frame(subj = i, task = task, cond = cond,
                        response = NA_integer_, stringsAsFactors = FALSE)

      # The task returns either binary responses or response times. The two
      # constants of the time only place the median at exp(6.2), about 493 ms.
      # A median split is invariant under any monotone transform, so neither
      # of them reaches the recovered d'.
      if (rt && task == "I") {
        lat <- stats::rnorm(n_trials, mean = eta, sd = 1)
        row$rt <- exp(6.2 - 0.25 * lat)
      } else {
        row$response <- stats::rbinom(n_trials, 1L, stats::pnorm(eta))
        if (rt) row$rt <- NA_real_
      }
      out[[k]] <- row
      k <- k + 1L
    }
  }

  # The function joins all trial blocks.
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}
