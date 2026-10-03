# stan.R
# Compile, cache and sample the Stan model of uSDT
# Author: Ricardo Rey-Sáez
# Last modified: 03-10-2026

# Internal functions

# Compiled models stay in the session, so later fits skip the cache.
.stan_models <- new.env(parent = emptyenv())

# This function returns the compiled model for a backend. The model compiles
# once per machine into the uSDT cache directory and is reused until its code
# changes; files left by older versions of the model are then removed. The
# option `uSDT.cache_dir` moves the cache, for instance to a temporary folder.
.stan_model <- function(backend) {
  source <- system.file("stan", "usdt.stan", package = "uSDT", mustWork = TRUE)
  stem <- paste0("usdt-", substr(unname(tools::md5sum(source)), 1L, 12L))
  key <- paste(backend, stem)
  if (!is.null(.stan_models[[key]])) return(.stan_models[[key]])

  dir <- getOption("uSDT.cache_dir", tools::R_user_dir("uSDT", "cache"))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  cached <- list.files(dir, pattern = "^usdt-", full.names = TRUE)
  unlink(cached[!startsWith(basename(cached), stem)], recursive = TRUE)
  file <- file.path(dir, paste0(stem, ".stan"))
  if (!file.exists(file)) file.copy(source, file)

  compiling <- "compiling the Stan model; this happens once per installation."
  model <- switch(backend,
    rstan = {
      rds <- file.path(dir, paste0(stem, ".rds"))
      model <- if (file.exists(rds)) {
        tryCatch(readRDS(rds), error = function(e) NULL)
      }
      if (is.null(model)) {
        .usdt_msg(compiling)
        model <- rstan::stan_model(file, model_name = "usdt")
        saveRDS(model, rds)
      }
      model
    },
    cmdstanr = {
      exe <- file.path(dir, paste0(stem, if (.Platform$OS.type == "windows")
        ".exe" else ""))
      if (!file.exists(exe)) .usdt_msg(compiling)
      cmdstanr::cmdstan_model(file, dir = dir, quiet = TRUE)
    })
  .stan_models[[key]] <- model
  model
}

# This function samples the model and returns the backend fit, the draws of the
# requested variables with their chains, and the sampler diagnostics.
.stan_sample <- function(model, data, settings, variables) {
  if (settings$backend == "rstan") {
    args <- list(object = model, data = data, chains = settings$chains,
                 iter = settings$iter, warmup = settings$warmup,
                 cores = settings$cores, control = settings$control,
                 refresh = settings$refresh, pars = c("z_d", "z_c"),
                 include = FALSE, show_messages = FALSE)
    if (!is.null(settings$seed)) args$seed <- settings$seed

    # uSDT reports its own diagnostics, so rstan's warnings would repeat them.
    fit <- suppressWarnings(do.call(rstan::sampling, args))
    params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    divergent <- sum(vapply(params, function(p) sum(p[, "divergent__"]), 0))
    treedepth <- sum(vapply(params, function(p) {
      sum(p[, "treedepth__"] >= settings$control$max_treedepth)
    }, 0))
    draws <- posterior::as_draws_array(as.array(fit))
  } else {
    fit <- model$sample(
      data = data, chains = settings$chains,
      parallel_chains = settings$cores, iter_warmup = settings$warmup,
      iter_sampling = settings$iter - settings$warmup, seed = settings$seed,
      adapt_delta = settings$control$adapt_delta,
      max_treedepth = settings$control$max_treedepth,
      refresh = settings$refresh, show_messages = FALSE,
      show_exceptions = FALSE)
    summary <- fit$diagnostic_summary(quiet = TRUE)
    divergent <- sum(summary$num_divergent)
    treedepth <- sum(summary$num_max_treedepth)

    # cmdstanr reads its output files lazily, so only uSDT's variables are read.
    written <- unique(sub("\\[.*$", "", fit$metadata()$variables))
    draws <- fit$draws(variables = intersect(variables, written))
  }

  # Only the variables uSDT reads are kept, whatever their size in this fit.
  present <- posterior::variables(draws)
  keep <- present[sub("\\[.*$", "", present) %in% variables]
  list(fit = fit, draws = posterior::subset_draws(draws, variable = keep),
       divergent = divergent, treedepth = treedepth)
}
