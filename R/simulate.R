# simulate.R
# This script simulates trial data from the hierarchical uSDT model.
# Author: Ricardo Rey-Sáez
# Last modified: 08-09-2026

# Public functions

#' Simulate data from a hierarchical SDT model
#'
#' Generates trial-level data for a direct and an indirect task. Every subject
#' has one sensitivity in each task, and the two sensitivities correlate across
#' subjects. This is the structure that [hsdt()] estimates, so the function is
#' useful to check an analysis before collecting data, or to study the power of
#' a planned design.
#'
#' @param n_subj Number of subjects.
#' @param n_trials Number of trials per subject and task. Half of them belong
#'   to the signal condition and half to the noise condition.
#' @param gamma_D,gamma_I Average sensitivity (d') of the direct and the
#'   indirect task.
#' @param sd_D,sd_I How much each sensitivity varies between subjects.
#' @param rho Correlation between the two sensitivities across subjects.
#' @param crit_D,sd_crit Model intercept of the direct task and how much it
#'   varies between subjects. The classical criterion has the opposite sign,
#'   `c = -crit_D`. See Details.
#' @param crit_I,sd_crit_I The same two values for the indirect task. The
#'   function sets both to zero when `rt = TRUE`, because a median split leaves
#'   no criterion to estimate. For a binary indirect task, keep `sd_crit_I`
#'   above zero. A value of zero gives every subject the same criterion, and
#'   the model then becomes singular when it tries to estimate that variation.
#' @param rt If `TRUE`, the indirect task returns response times instead of
#'   binary responses, so that the median split of [meyen_split()] has
#'   something to work on. The times fall as the evidence for signal rises, so
#'   faster responses correspond to signal. See Details.
#'
#' @return A data frame with one row per trial and the columns `subj`, `task`
#'   (`"D"` or `"I"`), `cond` (`1` for signal, `0` for noise) and `response`.
#'   With `rt = TRUE` it also has a numeric `rt` column, and `response` is `NA`
#'   in the indirect task.
#'
#' @details
#' The linear predictor is `crit_D + d' * S`, so `crit_D` is the model
#' intercept, which [hsdt()] reports as `c_D`. The classical criterion has the
#' opposite sign under deviation coding, `c = -crit_D`. That is the value
#' [sdt_moments()] returns and the value the printed summary shows. Simulating
#' `crit_D = 0.5` therefore gives a criterion of `-0.5`.
#'
#' # Simulated response times
#'
#' With `rt = TRUE` the indirect task returns `exp(6.2 - 0.25 * e)`, where `e`
#' is the latent evidence of that trial. The logarithm of the time is normal,
#' so the time itself follows a lognormal distribution, and it decreases as the
#' evidence grows. A faster response is therefore the signal response.
#'
#' The two constants only set the scale in milliseconds. They place the median
#' at `exp(6.2)`, around 493 ms, with the middle 95% of times between 302 and
#' 806 ms. A median split gives the same result under any transformation that
#' preserves the order of the values, so the recovered sensitivity does not
#' depend on them.
#'
#' The simulated times carry no variation beyond the evidence itself. This
#' makes `gamma_I` exactly the sensitivity of the dichotomized measure. Real
#' response times also vary for reasons unrelated to the discrimination, and
#' that extra variation would lower the sensitivity recovered from the split.
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

        # The exponent is the trial's own latent evidence, negated, so the
        # time is lognormal and falls as the evidence rises: a faster response
        # is the signal response. A median split is invariant under any
        # monotone transform, so neither constant reaches the recovered d'.
        # They only place the median at exp(6.2), about 493 ms, with the
        # middle 95% between 302 and 806 ms.
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
