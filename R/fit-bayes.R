# fit-bayes.R
# Fit the Bayesian hierarchical SDT model with Stan
# Author: Ricardo Rey-Sáez
# Last modified: 03-10-2026

# Internal functions

# The variables read from every fit: the population distributions and the
# subject-level sensitivities and criteria.
.bayes_variables <- c("mu_d", "sigma_d", "rho_d", "mu_c", "sigma_c", "rho_c",
                      "sigma_s", "d", "c")

# This function fits the Bayesian model and builds its hsdt object.
.hsdt_bayes <- function(data, fix_criteria, level, priors, unequal_variances,
                        backend, dots, call) {
  settings <- c(.bayes_settings(dots), backend = backend)
  rlang::check_installed(c(backend, "posterior"),
                         reason = "to fit the Bayesian model.")

  free_c <- .bayes_free_criteria(data, fix_criteria)
  stan <- .stan_data(data, free_c, unequal_variances, priors)
  sampled <- .stan_sample(.stan_model(settings$backend), stan$data, settings,
                          .bayes_variables)
  tests <- .bayes_tests(sampled$draws, level)
  diagnostics <- .bayes_diagnostics(sampled, tests, settings)

  if (length(diagnostics$issues)) {
    .usdt_warn("the posterior may be unreliable: ",
               paste(diagnostics$issues, collapse = "; "), ".\n",
               "  Increase `iter`, or raise `control = list(adapt_delta = ...)` ",
               "towards 1 when there are divergent transitions.")
  }

  structure(list(fit = sampled$fit, draws = sampled$draws,
                 subjects = stan$subjects, tests = tests, priors = priors,
                 design = list(criteria = c("c_D", "c_I")[free_c == 1L],
                               free_c = free_c,
                               unequal_variances = unequal_variances),
                 data = data, diagnostics = diagnostics, settings = settings,
                 call = call, level = level, estimation = "bayesian"),
            class = "hsdt")
}

# This function completes the sampling options given through `...`.
.bayes_settings <- function(dots) {
  defaults <- list(chains = 4L, iter = 3500L, warmup = 1000L,
                   cores = getOption("mc.cores", 1L), seed = NULL,
                   control = list(adapt_delta = 0.95, max_treedepth = 10L),
                   refresh = 0L)
  unknown <- setdiff(names(dots), names(defaults))
  if (length(unknown)) {
    .usdt_stop("`", unknown[1L], "` is not a sampling option. A Bayesian fit ",
               "accepts ", paste0("`", names(defaults), "`", collapse = ", "),
               ".")
  }
  unknown <- setdiff(names(dots$control), names(defaults$control))
  if (length(unknown)) {
    .usdt_stop("`control` accepts `adapt_delta` and `max_treedepth`, not `",
               unknown[1L], "`.")
  }
  s <- utils::modifyList(defaults, dots)

  for (arg in c("chains", "iter", "warmup", "cores")) {
    .check_scalar_number(s[[arg]], arg, lower = 1, whole = TRUE)
  }
  if (s$iter <= s$warmup) {
    .usdt_stop("`iter` counts the warmup too, so it must exceed `warmup`.")
  }
  if (!is.null(s$seed)) .check_scalar_number(s$seed, "seed", lower = 0, whole = TRUE)
  .check_scalar_number(s$refresh, "refresh", lower = 0, whole = TRUE)
  .check_scalar_number(s$control$adapt_delta, "adapt_delta", lower = 0,
                       upper = 1, open_lower = TRUE, open_upper = TRUE)
  .check_scalar_number(s$control$max_treedepth, "max_treedepth", lower = 1,
                       whole = TRUE)
  s
}

# This function decides which criteria the model estimates. The Stan model
# uses the classical criterion whatever the coding of the data, so a criterion
# fixed by the Meyen split is judged on the deviation scale.
.bayes_free_criteria <- function(data, fix_criteria) {
  if (fix_criteria == "none") return(c(1L, 1L))
  check <- .criterion_check(data$agg, "deviation")
  fixed <- vapply(c("direct", "indirect"), function(k) {
    data$meta$tasks[[k]]$dichotomized && check[[k]]$negligible
  }, TRUE)
  as.integer(!fixed)
}

# This function arranges the counts as subject x task matrices for Stan. A
# subject missing from one task keeps zero trials there, which adds nothing to
# the likelihood.
.stan_data <- function(data, free_c, unequal_variances, priors) {
  agg <- data$agg
  subjects <- unique(as.character(agg$subj))
  where <- cbind(match(as.character(agg$subj), subjects),
                 match(as.character(agg$task), c("D", "I")))
  counts <- function(signal, column) {
    m <- matrix(0L, length(subjects), 2L)
    rows <- agg$sig == signal
    m[where[rows, , drop = FALSE]] <- as.integer(agg[[column]][rows])
    m
  }
  list(subjects = subjects,
       data = c(list(I = length(subjects),
                     hit = counts(TRUE, "y"), fa = counts(FALSE, "y"),
                     N_hit = counts(TRUE, "n"), N_fa = counts(FALSE, "n"),
                     UV = as.integer(unequal_variances), free_c = free_c),
                .stan_priors(priors, free_c, unequal_variances)))
}

# This function draws the bivariate normal of the sensitivities from a fit and
# derives the tested quantities draw by draw.
.bayes_quantities <- function(draws) {
  m <- posterior::as_draws_matrix(
    posterior::subset_draws(draws, variable = c("mu_d", "sigma_d", "rho_d")))
  column <- function(v) as.numeric(m[, v])
  B <- cbind(mu_D = column("mu_d[1]"), mu_I = column("mu_d[2]"),
             sigma_D = column("sigma_d[1]"), sigma_I = column("sigma_d[2]"),
             rho = column("rho_d"))
  .usdt_quantities(.bivariate_primitives(B))
}

# This function summarises the posterior of the three hypotheses: the posterior
# mean and SD, the central credible interval, and the two-sided posterior
# p-value 2 min{P(q > 0), P(q < 0)}. That p-value falls below 1 - level exactly
# when zero lies outside the interval.
.bayes_tests <- function(draws, level) {
  Q <- .bayes_quantities(draws)
  shape <- c(posterior::niterations(draws), posterior::nchains(draws))
  a <- (1 - level) / 2
  terms <- c(diff = "d'(indirect) - d'(direct)", rho = "correlation",
             intercept = "intercept", slope = "slope")
  hypotheses <- c(diff = "H1", rho = "H2", intercept = "H3", slope = "H3")
  rows <- lapply(names(terms), function(q) {
    x <- Q[, q]
    chains <- matrix(x, shape[1L], shape[2L])
    limits <- stats::quantile(x, c(a, 1 - a), names = FALSE)
    data.frame(hypothesis = hypotheses[[q]], term = terms[[q]],
               estimate = mean(x), est.error = stats::sd(x),
               conf.low = limits[1L], conf.high = limits[2L],
               ci_method = "quantile",
               p.value = 2 * min(mean(x > 0), mean(x < 0)),
               rhat = posterior::rhat(chains),
               ess_bulk = posterior::ess_bulk(chains),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# This function collects the sampler diagnostics and names any problem.
.bayes_diagnostics <- function(sampled, tests, settings) {
  population <- posterior::subset_draws(
    sampled$draws, variable = intersect(
      c("mu_d", "sigma_d", "rho_d", "mu_c", "sigma_c", "rho_c", "sigma_s"),
      unique(sub("\\[.*$", "", posterior::variables(sampled$draws)))))
  summary <- posterior::summarise_draws(population, "rhat", "ess_bulk",
                                        "ess_tail")
  max_rhat <- max(summary$rhat, tests$rhat, na.rm = TRUE)
  min_bulk <- min(summary$ess_bulk, tests$ess_bulk, na.rm = TRUE)
  min_tail <- min(summary$ess_tail, na.rm = TRUE)

  issues <- c(
    if (sampled$divergent > 0)
      sprintf("%d divergent transitions", sampled$divergent),
    if (sampled$treedepth > 0)
      sprintf("%d transitions hit the maximum tree depth", sampled$treedepth),
    if (max_rhat > 1.01) sprintf("R-hat reaches %.3f", max_rhat),
    if (min(min_bulk, min_tail) < 400)
      sprintf("the effective sample size falls to %.0f", min(min_bulk, min_tail))
  )
  list(backend = settings$backend, chains = settings$chains,
       draws = posterior::ndraws(sampled$draws), max_rhat = max_rhat,
       min_ess_bulk = min_bulk, min_ess_tail = min_tail,
       divergent = sampled$divergent, treedepth = sampled$treedepth,
       issues = issues)
}
