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

# This function writes a bivariate normal given by its means, SDs and
# correlation (columns mu_D, mu_I, sigma_D, sigma_I and rho) as the five values
# read by .usdt_quantities().
.bivariate_primitives <- function(B) {
  cbind(gamma_D = B[, "mu_D"], gamma_I = B[, "mu_I"],
        s2_D = B[, "sigma_D"]^2, s2_I = B[, "sigma_I"]^2,
        s_DI = B[, "rho"] * B[, "sigma_D"] * B[, "sigma_I"])
}
