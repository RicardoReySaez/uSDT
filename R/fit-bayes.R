# fit-bayes.R
# Fit the Bayesian hierarchical SDT model with brms
# Author: Ricardo Rey-Sáez
# Last modified: 15-10-2026

# Internal functions

# This function fits the Bayesian model with brms and builds its hsdt object.
.hsdt_bayes <- function(data, fix_criteria, level, priors, unequal_variances,
                        parameterization, backend, dots, call) {
  settings <- c(.bayes_settings(dots), backend = backend)
  .check_backend(backend)

  design <- .bayes_design(data, fix_criteria, unequal_variances,
                          parameterization)
  model <- .brms_model(data, design, priors)
  if (settings$refresh > 0) .bayes_header(settings)
  fit <- .brms_sample(.brms_compiled(model, backend), model, settings)
  object <- structure(
    list(fit = fit, subjects = unique(as.character(data$agg$subj)),
         tests = NULL, priors = priors, design = design, data = data,
         diagnostics = NULL, settings = settings, call = call, level = level,
         estimation = "bayesian"),
    class = "hsdt")
  population <- .population_draws(object)
  tests <- .bayes_tests(.usdt_bivariate(population), level)
  diagnostics <- .bayes_diagnostics(fit, population, tests, settings)
  object$tests <- tests
  object$diagnostics <- diagnostics
  if (settings$refresh > 0) .bayes_report(diagnostics, settings)

  # The warning repeats any problem, so code that runs quietly still sees it,
  # and gives the remedy of each one.
  if (length(diagnostics$issues)) {
    remedies <- c(
      if (diagnostics$divergent > 0)
        "raise `control = list(adapt_delta = ...)` towards 1",
      if (diagnostics$treedepth > 0)
        "raise `control = list(max_treedepth = ...)`",
      if (diagnostics$max_rhat > 1.01 ||
          min(diagnostics$min_ess_bulk, diagnostics$min_ess_tail) < 400)
        "increase `iter`")
    .usdt_warn("the posterior may be unreliable: ",
               paste(diagnostics$issues, collapse = "; "), ".\n",
               "  To fix it, ", paste(remedies, collapse = ", or "), ".")
  }

  object
}

# This function checks that the chosen backend can compile the model. brms
# brings rstan, which compiles at run time against the headers of BH and
# RcppEigen that a binary install of rstan does not bring; cmdstanr, which is
# not on CRAN, needs CmdStan itself.
.check_backend <- function(backend) {
  needs <- switch(backend, rstan = c("BH", "RcppEigen"), cmdstanr = "cmdstanr")
  rlang::check_installed(needs, reason = paste0("to compile the model with ",
                                                backend, "."))
  if (backend == "cmdstanr" &&
      is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
    .usdt_stop("cmdstanr cannot find CmdStan. Install it once with ",
               "cmdstanr::install_cmdstan(), or use backend = \"rstan\".")
  }
  invisible(TRUE)
}

# This function announces the fit and how its chains run.
.bayes_header <- function(settings) {
  cli::cli_rule(left = "Fitting the Bayesian hierarchical SDT model",
                right = settings$backend)
  draws <- settings$iter - settings$warmup
  at_once <- min(settings$cores, settings$chains)
  cli::cli_alert_info(paste0(
    "{settings$chains} chain{?s} of {draws} draws after {settings$warmup} ",
    "warmup iterations, {at_once} at a time."))
}

# This function reports, after the progress bar, the sampler and convergence
# checks that .bayes_diagnostics() also turns into a warning.
.bayes_report <- function(diagnostics, settings) {
  depth <- settings$control$max_treedepth
  n <- diagnostics$divergent
  if (n == 0) {
    cli::cli_alert_success("No divergent transitions.")
  } else {
    cli::cli_alert_warning("{n} divergent transition{?s}.")
  }
  n <- diagnostics$treedepth
  if (n == 0) {
    cli::cli_alert_success(
      "No transition reached the maximum tree depth ({depth}).")
  } else {
    cli::cli_alert_warning(
      "{n} transition{?s} reached the maximum tree depth ({depth}).")
  }
  rhat <- sprintf("%.3f", diagnostics$max_rhat)
  if (diagnostics$max_rhat <= 1.01) {
    cli::cli_alert_success("Largest R-hat {rhat}, below 1.01.")
  } else {
    cli::cli_alert_warning("Largest R-hat {rhat}, above 1.01.")
  }
  smallest <- min(diagnostics$min_ess_bulk, diagnostics$min_ess_tail)
  ess <- round(smallest)
  if (smallest >= 400) {
    cli::cli_alert_success("Smallest effective sample size {ess}, above 400.")
  } else {
    cli::cli_alert_warning("Smallest effective sample size {ess}, below 400.")
  }
}

# This function returns how many chains run at once by default: up to four,
# leaving two cores free, and no more than two while R CMD check limits the
# cores a package may use.
.default_cores <- function() {
  n <- tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
  cores <- if (is.na(n)) 1L else max(1L, min(n - 2L, 4L))
  limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(limit) && !identical(tolower(limit), "false")) {
    cores <- min(cores, 2L)
  }
  as.integer(cores)
}

# This function completes the sampling options given through `...`. Each of
# the four chains keeps 5000 draws after 1000 warmup iterations, and reports
# its progress every 1% of its iterations.
.bayes_settings <- function(dots) {
  defaults <- list(chains = 4L, iter = 6000L, warmup = 1000L,
                   cores = .default_cores(), seed = NULL,
                   control = list(adapt_delta = 0.95, max_treedepth = 10L),
                   refresh = NULL)
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
  if (is.null(s$refresh)) s$refresh <- max(1L, s$iter %/% 100L)
  .check_scalar_number(s$refresh, "refresh", lower = 0, whole = TRUE)
  .check_scalar_number(s$control$adapt_delta, "adapt_delta", lower = 0,
                       upper = 1, open_lower = TRUE, open_upper = TRUE)
  .check_scalar_number(s$control$max_treedepth, "max_treedepth", lower = 1,
                       whole = TRUE)
  s
}

# This function describes the Bayesian model: which criteria it estimates
# (`free_c`, 1 for estimated, and their names), whether the signal SDs are
# free, and how it parameterizes the two sensitivities ("correlation" or
# "regression"). The model uses the classical criterion whatever the coding of
# the data, so a criterion fixed by the Meyen split is judged on the deviation
# scale.
.bayes_design <- function(data, fix_criteria, unequal_variances,
                          parameterization) {
  free_c <- if (fix_criteria == "none") c(1L, 1L) else {
    check <- .criterion_check(data$agg, "deviation")
    as.integer(!vapply(c("direct", "indirect"), function(k) {
      data$meta$tasks[[k]]$dichotomized && check[[k]]$negligible
    }, TRUE))
  }
  list(criteria = c("c_D", "c_I")[free_c == 1L], free_c = free_c,
       unequal_variances = unequal_variances,
       parameterization = parameterization)
}

# Compiled models stay in the session, so later fits skip the cache.
.brms_models <- new.env(parent = emptyenv())

# This function returns a brms fit without draws that holds the compiled model:
# from the session, from the cache directory, or compiled once. The key hashes
# the Stan code brms writes, which holds the priors and the brms version, and
# the backend. Each compiled model has its own folder, which for cmdstanr also
# holds the executable; using a model refreshes its folder, and only the most
# recently used ones are kept. The option `uSDT.cache_dir` moves the cache,
# for instance to a temporary folder.
.brms_compiled <- function(model, backend) {
  family <- stats::binomial("probit")
  code <- brms::stancode(model$formula, data = model$data, family = family,
                         prior = model$prior, backend = backend)
  key <- paste0("brms-", backend, "-", rlang::hash(code))
  if (!is.null(.brms_models[[key]])) return(.brms_models[[key]])

  cache <- getOption("uSDT.cache_dir", tools::R_user_dir("uSDT", "cache"))
  dir <- file.path(cache, key)
  file <- file.path(dir, "compiled.rds")
  compiled <- if (file.exists(file)) {
    tryCatch(readRDS(file), error = function(e) NULL)
  }
  if (is.null(compiled)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    cli::cli_alert_info(paste("Compiling the Stan model; this happens once",
                              "per model and set of priors."))
    old <- options(cmdstanr_write_stan_file_dir = dir)
    on.exit(options(old), add = TRUE)
    # brms checks the compiled model with one iteration, whose sampler warnings
    # say nothing about the fit.
    utils::capture.output(compiled <- suppressMessages(suppressWarnings(
      brms::brm(model$formula, data = model$data, family = family,
                prior = model$prior, chains = 0, backend = backend,
                silent = 2, refresh = 0,
                stan_model_args = if (backend == "cmdstanr") list(dir = dir)
                else list()))))
    saveRDS(compiled, file)
  }
  Sys.setFileTime(dir, Sys.time())
  .prune_cache(cache)
  .brms_models[[key]] <- compiled
  compiled
}

# This function keeps the most recently used compiled models in the cache and
# removes the others, with the files of the Stan model of earlier versions of
# uSDT. Nothing else in the folder is touched.
.prune_cache <- function(cache, keep = 20L) {
  entries <- list.files(cache, full.names = TRUE)
  models <- entries[startsWith(basename(entries), "brms-") & dir.exists(entries)]
  recent <- utils::head(models[order(file.mtime(models), decreasing = TRUE)],
                        keep)
  unlink(c(setdiff(models, recent),
           entries[startsWith(basename(entries), "usdt-")]), recursive = TRUE)
}

# This function samples each chain of the model in its own background session
# with brms, behind the progress bar of .run_chains(), and combines the chains
# into one brmsfit. Chains share one seed; the chain number gives each its own
# stream.
.brms_sample <- function(compiled, model, settings) {
  seed <- settings$seed %||% sample.int(.Machine$integer.max, 1L)
  args <- list(newdata = model$data, recompile = FALSE, chains = 1L,
               cores = 1L, iter = settings$iter, warmup = settings$warmup,
               seed = seed, control = settings$control,
               refresh = settings$refresh, silent = 1L)
  id <- if (settings$backend == "rstan") "chain_id" else "chain_ids"
  chains <- .run_chains(settings, function(chain) {
    callr::r_bg(function(fit, args) {
      loadNamespace("brms")
      suppressWarnings(do.call(stats::update, c(list(fit), args)))
    }, args = list(fit = compiled,
                   args = c(args, stats::setNames(list(chain), id))))
  })
  # Combining the chains loads rstan, which brms announces.
  suppressMessages(brms::combine_models(mlist = chains, check_data = FALSE))
}

# This function returns the draws of the population parameters of a Bayesian
# fit under their uSDT names: mu_d, sigma_d and rho_d for the two
# sensitivities, mu_c, sigma_c and rho_c for the estimated criteria, and
# sigma_s under unequal variances. The latent regression samples the direct
# mean and SD with the intercept, slope and residual SD, from which the
# indirect mean and SD and the correlation follow draw by draw.
.population_draws <- function(object) {
  names <- .brms_names(object$design)
  draws <- posterior::as_draws_array(object$fit, variable = unname(names))
  posterior::variables(draws) <- names(names)
  if (object$design$parameterization == "regression") {
    get <- function(name) {
      as.vector(posterior::extract_variable_matrix(draws, name))
    }
    slope <- get("slope")
    sigma_D <- get("sigma_d[1]")
    sigma_I <- sqrt(slope^2 * sigma_D^2 + get("sigma_e")^2)
    derived <- posterior::draws_array(
      `mu_d[2]` = get("intercept") + slope * get("mu_d[1]"),
      `sigma_d[2]` = sigma_I, rho_d = slope * sigma_D / sigma_I,
      .nchains = posterior::nchains(draws))
    keep <- c("mu_d[1]", "mu_d[2]", "sigma_d[1]", "sigma_d[2]", "rho_d",
              setdiff(names(names), c("mu_d[1]", "sigma_d[1]", "intercept",
                                      "slope", "sigma_e")))
    draws <- posterior::subset_draws(posterior::bind_draws(draws, derived),
                                     variable = keep)
  }
  draws
}

# The names given to the bivariate normal of the sensitivities, whatever model
# the draws come from.
.bivariate_names <- c("mu_D", "mu_I", "sigma_D", "sigma_I", "rho")

# This function selects the bivariate normal of the sensitivities from draws,
# given the names of its five variables in that order, and gives them the
# common names.
.select_bivariate <- function(draws, variables) {
  B <- posterior::subset_draws(draws, variable = variables)
  posterior::variables(B) <- .bivariate_names
  B
}

# These functions give, draw by draw, the five values that describe the two
# sensitivities and the tested quantities derived from them. The draws are
# stacked chain after chain.
.draw_primitives <- function(bivariate) {
  .bivariate_primitives(unclass(posterior::as_draws_matrix(bivariate)))
}
.draw_quantities <- function(bivariate) {
  .usdt_quantities(.draw_primitives(bivariate))
}

# This function returns the draws of one subject-level variable of a Bayesian
# fit, `d` or `c`, in one task (1 direct, 2 indirect): one column per subject,
# in the order of `subjects`. A subject's value is the population coefficient
# plus the subject's effect; under the latent regression the indirect d' is
# b0 + b1 * d'(direct) + e. A criterion fixed by the Meyen split makes
# HR + FAR = 1, which is 0 under equal variances.
.subject_draws <- function(object, variable, task, subjects = object$subjects) {
  design <- object$design
  k <- c("D", "I")[task]
  read <- function(variables) {
    unclass(posterior::as_draws_matrix(object$fit, variable = variables))
  }
  effects <- function(nlpar, group) {
    sprintf("r_%s__%s[%s,Intercept]", group, nlpar, subjects)
  }
  coefficient <- function(nlpar, group = "subj") {
    m <- read(c(paste0("b_", nlpar, "_Intercept"), effects(nlpar, group)))
    m[, -1L, drop = FALSE] + m[, 1L]
  }

  out <- if (variable == "d" && task == 2L &&
             design$parameterization == "regression") {
    m <- read(c("b_b0_Intercept", "b_b1_Intercept", effects("e", "subj")))
    m[, 1L] + m[, 2L] * coefficient("dD") + m[, -(1:2), drop = FALSE]
  } else if (variable == "d") {
    coefficient(paste0("d", k))
  } else if (design$free_c[task] == 1L) {
    coefficient(paste0("c", k), "subj_c")
  } else {
    s <- if (design$unequal_variances) {
      as.vector(read(paste0("b_s", k, "_Intercept")))
    } else {
      1
    }
    .subject_draws(object, "d", task, subjects) * (1 - s) / (2 * (1 + s))
  }
  dimnames(out) <- list(NULL, subjects)
  out
}

# This function selects the bivariate normal of the sensitivities from the
# population draws of a Bayesian fit.
.usdt_bivariate <- function(population) {
  .select_bivariate(population, c("mu_d[1]", "mu_d[2]", "sigma_d[1]",
                                  "sigma_d[2]", "rho_d"))
}

# This function summarises the posterior of the three hypotheses from the draws
# of the bivariate normal: the posterior mean and SD, the central credible
# interval, and the two-sided posterior p-value 2 min{P(q > 0), P(q < 0)}. That
# p-value falls below 1 - level exactly when zero lies outside the interval.
.bayes_tests <- function(bivariate, level) {
  Q <- .draw_quantities(bivariate)
  shape <- c(posterior::niterations(bivariate), posterior::nchains(bivariate))
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

# This function collects the sampler diagnostics of the fit and the
# convergence of its population parameters and tested quantities, and names
# any problem.
.bayes_diagnostics <- function(fit, population, tests, settings) {
  sampler <- brms::nuts_params(fit)
  divergent <- sum(sampler$Value[sampler$Parameter == "divergent__"])
  treedepth <- sum(sampler$Value[sampler$Parameter == "treedepth__"] >=
                     settings$control$max_treedepth)
  summary <- posterior::summarise_draws(population, "rhat", "ess_bulk",
                                        "ess_tail")
  max_rhat <- max(summary$rhat, tests$rhat, na.rm = TRUE)
  min_bulk <- min(summary$ess_bulk, tests$ess_bulk, na.rm = TRUE)
  min_tail <- min(summary$ess_tail, na.rm = TRUE)

  issues <- c(
    if (divergent > 0) sprintf("%d divergent transitions", divergent),
    if (treedepth > 0)
      sprintf("%d transitions hit the maximum tree depth", treedepth),
    if (max_rhat > 1.01) sprintf("R-hat reaches %.3f", max_rhat),
    if (min(min_bulk, min_tail) < 400)
      sprintf("the effective sample size falls to %.0f", min(min_bulk, min_tail))
  )
  list(backend = settings$backend, chains = settings$chains,
       draws = posterior::ndraws(population), max_rhat = max_rhat,
       min_ess_bulk = min_bulk, min_ess_tail = min_tail,
       divergent = divergent, treedepth = treedepth, issues = issues)
}

# This function runs every chain in its own background R session, at most
# `cores` at a time, and returns what each session returned. Stan's messages
# and warnings stay in those sessions; the iterations they report move a single
# progress bar, which `refresh = 0` hides. Interrupting the fit stops them.
.run_chains <- function(settings, start) {
  chains <- settings$chains
  done <- integer(chains)
  results <- vector("list", chains)
  queue <- seq_len(chains)
  running <- list()
  on.exit(for (p in running) p$kill(), add = TRUE)

  # The bar redraws itself in place, so it shows only in a console that can
  # do that; in knitr documents and logs a single line reports the time taken.
  # It starts with the first progress report after the first iteration, so
  # the seconds the sessions take to start do not distort its estimate of the
  # time left. It shows at once, in light blue where the console has colours,
  # and names the phase the chains are in.
  bar <- NULL
  label <- if (chains == 1L) "1 chain" else paste(chains, "chains")
  dynamic <- settings$refresh > 0 && cli::is_dynamic_tty()
  started <- Sys.time()
  if (dynamic) {
    old <- options(cli.progress_show_after = 0,
                   cli.progress_bar_style = .bar_style())
    on.exit(options(old), add = TRUE)
  }
  while (length(queue) || length(running)) {
    while (length(queue) && length(running) < settings$cores) {
      running[[as.character(queue[1L])]] <- start(queue[1L])
      queue <- queue[-1L]
    }
    for (key in names(running)) {
      p <- running[[key]]
      chain <- as.integer(key)
      alive <- p$is_alive()
      lines <- p$read_output_lines()
      p$read_error_lines()
      at <- regmatches(lines, regexpr("Iteration:\\s*[0-9]+", lines))
      if (length(at)) done[chain] <- max(as.integer(sub("\\D+", "", at)))
      if (!alive) {
        results[[chain]] <- tryCatch(p$get_result(), error = function(e) {
          .usdt_stop("chain ", chain, " failed: ",
                     conditionMessage(e$parent %||% e))
        })
        done[chain] <- settings$iter
        running[[key]] <- NULL
      }
    }
    phase <- if (all(done > settings$warmup)) "Sampling" else "Warmup  "
    if (dynamic && is.null(bar) && max(done) > 1L) {
      bar <- cli::cli_progress_bar(
        total = chains * settings$iter, clear = FALSE, auto_terminate = FALSE,
        status = phase,
        format = paste("{cli::pb_spin} {cli::pb_status} {cli::pb_bar}",
                       "{cli::pb_percent} | ETA: {cli::pb_eta}"),
        format_done = paste("{cli::col_green(cli::symbol$tick)} Sampled",
                            "{label} in {cli::pb_elapsed}."))
    }
    if (!is.null(bar)) {
      cli::cli_progress_update(id = bar, set = sum(done), status = phase)
    }
    if (length(running)) Sys.sleep(0.1)
  }
  if (!is.null(bar)) {
    cli::cli_progress_done(id = bar)
  } else if (settings$refresh > 0) {
    secs <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    elapsed <- if (secs < 60) sprintf("%.1fs", secs) else
      sprintf("%dm %ds", as.integer(secs %/% 60), as.integer(round(secs %% 60)))
    cli::cli_alert_success("Sampled {label} in {elapsed}.")
  }
  results
}

# This function returns the bar: light blue squares for the iterations run and
# a grey line for those to come. Their shapes differ, so the bar reads on light
# and dark consoles and without colours; a console that cannot draw them keeps
# cli's plain bar (NULL).
.bar_style <- function() {
  if (!cli::is_utf8_output()) return(NULL)
  done <- cli::make_ansi_style("#6CB4EE")("\u25a0")
  list(complete = done, current = done,
       incomplete = cli::make_ansi_style("#9AA0A6")("\u2500"))
}
