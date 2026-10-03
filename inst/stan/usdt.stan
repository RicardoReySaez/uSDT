// usdt.stan
// Hierarchical signal detection theory model 
// Author: Ricardo Rey-Sáez
// Last modified: 03-10-2026

functions {
  // Scaled beta log-probability density function
  real scaled_beta_lpdf(real rho, real a, real b) {
    return (a - 1) * log1p(rho) + (b - 1) * log1m(rho)
           - (a + b - 1) * log2() - lbeta(a, b);
  }
}

data {
  int<lower=1> I;                          // Number of subjects
  array[I, 2] int<lower=0> hit;            // Hits per subject and task
  array[I, 2] int<lower=0> fa;             // False alarms per subject and task
  array[I, 2] int<lower=0> N_hit;          // Signal trials (0 if the subject lacks the task)
  array[I, 2] int<lower=0> N_fa;           // Noise trials (0 if the subject lacks the task)
  int<lower=0, upper=1> UV;                // Unequal variances indicator
  array[2] int<lower=0, upper=1> free_c;   // 0: criterion fixed by the Meyen split
  int<lower=0, upper=1> prior_only;        // 1: sample from the priors alone

  // Prior hyperparameters
  vector[2] mu_d_loc;                      // Normal priors on the mean d'
  vector<lower=0>[2] mu_d_scale;
  vector<lower=0>[2] sd_d_df;              // Half-Student-t priors on the d' SDs
  vector<lower=0>[2] sd_d_scale;
  real<lower=0> rho_d_a;                   // Scaled-beta prior on the d' correlation
  real<lower=0> rho_d_b;
  vector[sum(free_c)] mu_c_loc;            // Priors for the estimated criteria only
  vector<lower=0>[sum(free_c)] mu_c_scale;
  vector<lower=0>[sum(free_c)] sd_c_scale;
  real<lower=0> rho_c_a;                   // Used only when both criteria are free
  real<lower=0> rho_c_b;
  vector<lower=0>[UV ? 2 : 0] sd_s_scale;  // Lognormal scales of the signal SDs
}

parameters {
  // Population-level distribution of d'
  vector[2] mu_d;
  vector<lower=0>[2] sigma_d;
  real<lower=-1, upper=1> rho_d;

  // Population-level distribution of the estimated criteria
  vector[sum(free_c)] mu_c;
  vector<lower=0>[sum(free_c)] sigma_c;
  array[sum(free_c) == 2] real<lower=-1, upper=1> rho_c;

  // Signal SDs, present only under unequal variances
  vector<lower=0>[UV ? 2 : 0] sigma_s_free;

  // Standardized subject-level random effects
  matrix[I, 2] z_d;
  matrix[I, sum(free_c)] z_c;
}

transformed parameters {
  vector[2] sigma_s = UV ? sigma_s_free : rep_vector(1, 2);

  // Subject-level d': the 2 x 2 Cholesky factor written out
  matrix[I, 2] d;
  d[, 1] = mu_d[1] + sigma_d[1] * z_d[, 1];
  d[, 2] = mu_d[2] + sigma_d[2] * 
          (rho_d * z_d[, 1] + sqrt(1 - square(rho_d)) * z_d[, 2]);

  // Subject-level criteria. A task split at its median has HR + FAR = 1, which
  // fixes c = d (1 - sigma_s) / (2 * (1 + sigma_s)), that is c = 0 with EV
  matrix[I, 2] c;
  {
    int k = 0;
    for (j in 1:2) {
      if (free_c[j]) {
        k += 1;
        if (k == 1) {
          c[, j] = mu_c[k] + sigma_c[k] * z_c[, k];
        } else {
          c[, j] = mu_c[k] + sigma_c[k] * 
                  (rho_c[1] * z_c[, 1] + sqrt(1 - square(rho_c[1])) * z_c[, 2]);
        }
      } else {
        c[, j] = d[, j] * (1 - sigma_s[j]) / (2 * (1 + sigma_s[j]));
      }
    }
  }
}

model {
  // Model priors
  mu_d    ~ normal(mu_d_loc, mu_d_scale);
  sigma_d ~ student_t(sd_d_df, 0, sd_d_scale);
  rho_d   ~ scaled_beta(rho_d_a, rho_d_b);
  mu_c    ~ normal(mu_c_loc, mu_c_scale);
  sigma_c ~ normal(0, sd_c_scale);
  for (k in 1:size(rho_c)) {
    rho_c[k] ~ scaled_beta(rho_c_a, rho_c_b);
  }
  sigma_s_free   ~ lognormal(0, sd_s_scale);
  to_vector(z_d) ~ std_normal();
  to_vector(z_c) ~ std_normal();

  // Model Log-Likelihood on aggregated counts
  if (!prior_only) {
    for (j in 1:2) {
      hit[, j] ~ binomial(N_hit[, j], Phi((0.5 * d[, j] - c[, j]) / sigma_s[j]));
      fa[, j]  ~ binomial(N_fa[, j],  Phi(-0.5 * d[, j] - c[, j]));
    }
  }
}

generated quantities {
  // Group-level sensitivity difference and latent regression
  real delta = mu_d[2] - mu_d[1];
  real slope = rho_d * sigma_d[2] / sigma_d[1];
  real intercept = mu_d[2] - slope * mu_d[1];
}
