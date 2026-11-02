# quantities.R
# This script defines the quantities tested by the three hypotheses.
# Author: Ricardo Rey-Sáez
# Last modified: 03-10-2026

# Internal functions

# This function derives the tested quantities from the five values that
# describe the two sensitivities: their means (gamma_A, gamma_I), variances
# (s2_A, s2_I) and covariance (s_AI). `P` is a named vector or a matrix with one
# row per set of values, so the same definitions serve the fitted estimates and
# every resampled or simulated set.
.usdt_quantities <- function(P) {
  if (is.null(dim(P))) P <- t(P)
  slope <- P[, "s_AI"] / P[, "s2_A"]
  cbind(diff      = P[, "gamma_I"] - P[, "gamma_A"],
        rho       = pmax(-1, pmin(1, P[, "s_AI"] /
                                       sqrt(P[, "s2_A"] * P[, "s2_I"]))),
        intercept = P[, "gamma_I"] - slope * P[, "gamma_A"],
        slope     = slope)
}

# This function writes a bivariate normal given by its means, SDs and
# correlation (columns mu_A, mu_I, sigma_A, sigma_I and rho) as the five values
# read by .usdt_quantities().
.bivariate_primitives <- function(B) {
  cbind(gamma_A = B[, "mu_A"], gamma_I = B[, "mu_I"],
        s2_A = B[, "sigma_A"]^2, s2_I = B[, "sigma_I"]^2,
        s_AI = B[, "rho"] * B[, "sigma_A"] * B[, "sigma_I"])
}
