# brms-model.R
# Write the Bayesian hierarchical SDT model as a brms model
# Author: Ricardo Rey-Sáez
# Last modified: 15-10-2026

# Internal functions

# This function builds the brms model of a design (from .bayes_design()): its
# non-linear formula, its data and its priors.
.brms_model <- function(data, design, priors) {
  list(formula = .brms_formula(design), data = .brms_data(data$agg),
       prior = .brms_prior(priors, design))
}

# This function writes the aggregated counts with the columns the formula
# reads. Each row is one task and condition of one subject: d_D and d_I are
# +0.5 (signal) or -0.5 (noise) in the rows of their task and 0 elsewhere;
# c_D and c_I are -1 in the rows of their task, so that their coefficients are
# criteria on the classical scale; sig, D and I mark signal rows and tasks.
# The criteria vary over subj_c, a copy of subj, so that brms gives them a
# correlation prior of their own.
.brms_data <- function(agg) {
  D <- as.numeric(agg$task == "D")
  I <- 1 - D
  condition <- ifelse(agg$sig, 0.5, -0.5)
  subj <- factor(agg$subj)
  data.frame(y = agg$y, n = agg$n, subj = subj, subj_c = subj,
             sig = as.numeric(agg$sig), D = D, I = I, c_D = -D, c_I = -I,
             d_D = condition * D, d_I = condition * I)
}

# This function writes the model as a brms non-linear formula on the counts.
# The linear predictor of a row is d/2 - c for signal and -d/2 - c for noise,
# and under unequal variances a signal row is divided by the signal SD of its
# task. A criterion fixed by the Meyen split is 0 under equal variances and
# d (1 - s) / (2 (1 + s)) under unequal ones, the value that makes
# HR + FAR = 1. The correlation parameterization gives the two sensitivities
# correlated subject effects; the latent regression writes the indirect d' of
# each subject as b0 + b1 * dD + e, with e an independent subject effect.
.brms_formula <- function(design) {
  tasks <- c("D", "I")
  regression <- design$parameterization == "regression"
  uv <- design$unequal_variances
  d <- c(D = "dD", I = if (regression) "(b0 + b1 * dD + e)" else "dI")
  criterion <- vapply(seq_along(tasks), function(j) {
    k <- tasks[j]
    if (design$free_c[j] == 1L) return(paste0("c", k))
    if (uv) return(paste0(d[[k]], " * (1 - s", k, ") / (2 * (1 + s", k, "))"))
    NA_character_
  }, "")
  eta <- paste(c(paste0("d_", tasks, " * ", d),
                 paste0("c_", tasks, " * ", criterion)[!is.na(criterion)]),
               collapse = " + ")
  if (uv) {
    eta <- paste0("(", eta, ") / (1 + sig * (D * (sD - 1) + I * (sI - 1)))")
  }

  free <- tasks[design$free_c == 1L]
  nlpars <- c(
    if (regression) {
      c("dD ~ 1 + (1 | subj)", "b0 ~ 1", "b1 ~ 1", "e ~ 0 + (1 | subj)")
    } else {
      c("dD ~ 1 + (1 | p | subj)", "dI ~ 1 + (1 | p | subj)")
    },
    sprintf("c%s ~ 1 + (1 | %ssubj_c)", free,
            if (length(free) == 2L) "q | " else ""),
    if (uv) c("sD ~ 1", "sI ~ 1"))
  brms::bf(stats::as.formula(paste("y | trials(n) ~", eta), env = baseenv()),
           flist = lapply(nlpars, stats::as.formula, env = baseenv()),
           nl = TRUE)
}

# This function writes each prior a design uses as a brms prior on its own
# parameter, so that the prior table of the fit lists one per parameter and
# task. A signal SD is a population-level coefficient, bounded at zero.
.brms_prior <- function(priors, design) {
  used <- .used_priors(priors, design)
  target <- function(parameter, task) {
    k <- c(direct = "D", indirect = "I", both = "")[[task]]
    switch(parameter,
      dprime = list(class = "b", nlpar = paste0("d", k)),
      sd_dprime = list(class = "sd", group = "subj", nlpar = paste0("d", k)),
      cor_dprime = list(class = "cor", group = "subj"),
      criterion = list(class = "b", nlpar = paste0("c", k)),
      sd_criterion = list(class = "sd", group = "subj_c",
                          nlpar = paste0("c", k)),
      cor_criterion = list(class = "cor", group = "subj_c"),
      sd_signal = list(class = "b", nlpar = paste0("s", k), lb = 0),
      intercept = list(class = "b", nlpar = "b0"),
      slope = list(class = "b", nlpar = "b1"),
      sd_residual = list(class = "sd", group = "subj", nlpar = "e"))
  }
  do.call(c, lapply(seq_len(nrow(used)), function(i) {
    do.call(brms::set_prior, c(list(prior = used$prior[i]),
                               target(used$parameter[i], used$task[i])))
  }))
}

# This function names, in a brms fit of the design, the population parameters
# the model samples, by their uSDT names: the means and SDs of d' and their
# correlation (under the latent regression, the direct mean and SD with the
# intercept, slope and residual SD), the means, SDs and correlation of the
# estimated criteria, and the signal SDs.
.brms_names <- function(design) {
  tasks <- c("D", "I")
  free <- tasks[design$free_c == 1L]
  b <- function(nlpar) paste0("b_", nlpar, "_Intercept")
  sd <- function(group, nlpar) paste0("sd_", group, "__", nlpar, "_Intercept")
  sensitivities <- if (design$parameterization == "regression") {
    c("mu_d[1]" = b("dD"), "sigma_d[1]" = sd("subj", "dD"),
      intercept = b("b0"), slope = b("b1"), sigma_e = sd("subj", "e"))
  } else {
    c("mu_d[1]" = b("dD"), "mu_d[2]" = b("dI"),
      "sigma_d[1]" = sd("subj", "dD"), "sigma_d[2]" = sd("subj", "dI"),
      rho_d = "cor_subj__dD_Intercept__dI_Intercept")
  }
  k <- seq_along(free)
  criteria <- c(
    stats::setNames(b(paste0("c", free)), sprintf("mu_c[%d]", k)),
    stats::setNames(sd("subj_c", paste0("c", free)), sprintf("sigma_c[%d]", k)),
    if (length(free) == 2L) c("rho_c[1]" = "cor_subj_c__cD_Intercept__cI_Intercept"))
  signal <- if (design$unequal_variances) {
    c("sigma_s[1]" = b("sD"), "sigma_s[2]" = b("sI"))
  }
  c(sensitivities, criteria, signal)
}
