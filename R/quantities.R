# quantities.R
# This script defines the quantities tested by the three hypotheses.
# Author: Ricardo Rey-Sáez
# Last modified: 03-10-2026

# Internal functions

# This function derives the tested quantities from the five values that
# describe the two sensitivities: their means (gamma_D, gamma_I), variances
# (s2_D, s2_I) and covariance (s_DI). `P` is a named vector or a matrix with one
# row per set of values, so the same definitions serve the fitted estimates and
# every resampled or simulated set.
.usdt_quantities <- function(P) {
  if (is.null(dim(P))) P <- t(P)
  slope <- P[, "s_DI"] / P[, "s2_D"]
  cbind(diff      = P[, "gamma_I"] - P[, "gamma_D"],
        rho       = pmax(-1, pmin(1, P[, "s_DI"] / sqrt(P[, "s2_D"] * P[, "s2_I"]))),
        intercept = P[, "gamma_I"] - slope * P[, "gamma_D"],
        slope     = slope)
}
