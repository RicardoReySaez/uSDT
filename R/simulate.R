# simulate.R
# This script simulates trial data from the hierarchical uSDT model.
# Author: Ricardo Rey-Sáez
# Last modified: 04-09-2026

# Public functions

#' Simulate data from a hierarchical SDT model
#'
#' Generates trial-level data for a direct and an indirect task whose
#' subject-level sensitivities are correlated, which is the structure
#' [hsdt()] estimates.
#'
#' @param n_subj Number of subjects.
#' @param n_trials Number of trials per subject and task, split evenly
#'   between the signal and the noise condition.
#' @param gamma_D,gamma_I Group-level sensitivity (d') of the direct and the
#'   indirect task.
#' @param sd_D,sd_I Between-subject standard deviation of each sensitivity.
#' @param rho Latent correlation between the two sensitivities.
#' @param crit_D,sd_crit Group-level *model intercept* of the direct task and
#'   its between-subject standard deviation. The classical criterion is its
#'   negative, `c = -crit_D`. See Details.
#' @param crit_I,sd_crit_I The same for the indirect task. Both are forced to
#'   zero when `rt = TRUE`, because a Meyen median split leaves no criterion to
#'   estimate. Leaving them at zero for a *binary* indirect task produces data
#'   that no model with an indirect criterion can fit without becoming
#'   singular, so the defaults mirror the direct task.
#' @param rt If `TRUE`, the indirect task is returned as response times rather
#'   than as a binary response, so that the median split of [meyen_split()] can
#'   be exercised. The times are a monotone decreasing transform of the latent
#'   evidence, so faster responses correspond to signal.
#'
#' @return A data frame with one row per trial and columns `subj`, `task`
#'   (`"D"` or `"I"`), `cond` (`1` signal, `0` noise) and `response`. When
#'   `rt = TRUE` a numeric `rt` column is added, and `response` is `NA` for the
#'   indirect task.
#'
#' @details
#' The linear predictor is `crit_D + d' * S`, so `crit_D` is the model
#' intercept, the quantity [hsdt()] estimates as `c_D`. The classical
#' criterion runs the other way: `c = -crit_D` under deviation coding, which is
#' what [sdt_moments()] returns and what the printed summary reports. Simulating
#' `crit_D = 0.5` therefore recovers a criterion of `-0.5`.
#'
#' @seealso [hsdt()], [sdt_moments()]
#'
#' @examples
#' set.seed(1)
#' df <- usdt_simulate(n_subj = 20, n_trials = 60)
#' head(df)
#' table(df$task, df$cond)
#'
#' @export
usdt_simulate <- function(n_subj = 50, n_trials = 100,
                          gamma_D = 0.8, gamma_I = 0.3,
                          sd_D = 0.5, sd_I = 0.3, rho = 0.5,
                          crit_D = 0, sd_crit = 0.3,
                          crit_I = 0, sd_crit_I = 0.3,
                          rt = FALSE) {

  # The function checks the parameters that have a restricted range.
  .check_scalar_number(n_subj, "n_subj", lower = 1, whole = TRUE)
  .check_scalar_number(n_trials, "n_trials", lower = 2, whole = TRUE)
  for (arg in c("gamma_D", "gamma_I", "crit_D", "crit_I")) {
    .check_scalar_number(get(arg), arg)
  }
  .check_scalar_number(rho, "rho", lower = -1, upper = 1,
                       open_lower = TRUE, open_upper = TRUE)
  .check_scalar_number(sd_D, "sd_D", lower = 0, open_lower = TRUE)
  .check_scalar_number(sd_I, "sd_I", lower = 0, open_lower = TRUE)
  .check_scalar_number(sd_crit, "sd_crit", lower = 0)
  .check_scalar_number(sd_crit_I, "sd_crit_I", lower = 0)
  if (!is.logical(rt) || length(rt) != 1L || is.na(rt)) {
    .usdt_stop("`rt` must be `TRUE` or `FALSE`.")
  }
  if (n_trials %% 2L != 0L) {
    .usdt_stop("`n_trials` must be even so that the two conditions are balanced.")
  }

  # The model draws two related sensitivities for each subject.
  S <- matrix(c(sd_D^2, rho * sd_D * sd_I,
                rho * sd_D * sd_I, sd_I^2), nrow = 2L)
  u  <- matrix(stats::rnorm(n_subj * 2L), n_subj) %*% chol(S)

  # A median split leaves the indirect criterion at zero without variation.
  if (rt) { crit_I <- 0; sd_crit_I <- 0 }
  uc <- list(D = stats::rnorm(n_subj, 0, sd_crit),
             I = stats::rnorm(n_subj, 0, sd_crit_I))

  # Each subject receives a balanced block for both tasks.
  half <- n_trials %/% 2L
  cond <- rep(c(0L, 1L), each = half)
  dev  <- cond - 0.5

  out <- vector("list", 2L * n_subj)
  k <- 1L
  for (task in c("D", "I")) {
    g  <- if (task == "D") gamma_D else gamma_I
    uu <- if (task == "D") u[, 1L] else u[, 2L]
    cc <- if (task == "D") crit_D + uc$D else crit_I + uc$I

    for (i in seq_len(n_subj)) {

      # The model calculates the evidence for every trial.
      eta <- cc[i] + (g + uu[i]) * dev
      row <- data.frame(subj = i, task = task, cond = cond,
                        response = NA_integer_, stringsAsFactors = FALSE)

      # The task returns either binary responses or response times.
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
