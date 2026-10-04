# stan.R
# Compile, cache and sample the Stan model of uSDT
# Author: Ricardo Rey-Sáez
# Last modified: 04-10-2026

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

  compiling <- function() {
    cli::cli_alert_info(
      "Compiling the Stan model; this happens once per installation.")
  }
  model <- switch(backend,
    rstan = {
      rds <- file.path(dir, paste0(stem, ".rds"))
      model <- if (file.exists(rds)) {
        tryCatch(readRDS(rds), error = function(e) NULL)
      }
      if (is.null(model)) {
        compiling()
        model <- rstan::stan_model(file, model_name = "usdt")
        saveRDS(model, rds)
      }
      model
    },
    cmdstanr = {
      exe <- file.path(dir, paste0(stem, if (.Platform$OS.type == "windows")
        ".exe" else ""))
      if (!file.exists(exe)) compiling()
      cmdstanr::cmdstan_model(file, dir = dir, quiet = TRUE)
    })
  .stan_models[[key]] <- model
  model
}

# This function samples the model and returns the backend fit, the draws of the
# requested variables with their chains, and the sampler diagnostics.
.stan_sample <- function(model, data, settings, variables) {
  # Chains share one seed; the chain number gives each its own stream.
  seed <- settings$seed %||% sample.int(.Machine$integer.max, 1L)

  if (settings$backend == "rstan") {
    args <- list(data = data, chains = 1L, iter = settings$iter,
                 warmup = settings$warmup, seed = seed,
                 control = settings$control, refresh = settings$refresh,
                 pars = c("z_d", "z_c"), include = FALSE,
                 show_messages = FALSE)
    chains <- .run_chains(settings, function(chain) {
      callr::r_bg(function(model, args) {
        loadNamespace("rstan")
        args <- c(list(object = model), args)
        suppressWarnings(do.call(rstan::sampling, args))
      }, args = list(model = model, args = c(args, chain_id = chain)))
    })
    fit <- rstan::sflist2stanfit(chains)
    params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    divergent <- sum(vapply(params, function(p) sum(p[, "divergent__"]), 0))
    treedepth <- sum(vapply(params, function(p) {
      sum(p[, "treedepth__"] >= settings$control$max_treedepth)
    }, 0))
    draws <- posterior::as_draws_array(as.array(fit))
  } else {
    # The chains write their draws to this session's temporary folder, which
    # outlives the processes that sampled them.
    args <- list(data = data, chains = 1L, seed = seed,
                 iter_warmup = settings$warmup,
                 iter_sampling = settings$iter - settings$warmup,
                 adapt_delta = settings$control$adapt_delta,
                 max_treedepth = settings$control$max_treedepth,
                 refresh = settings$refresh, show_messages = TRUE,
                 show_exceptions = FALSE, output_dir = tempdir())
    path <- cmdstanr::cmdstan_path()
    files <- .run_chains(settings, function(chain) {
      callr::r_bg(function(model, args, path) {
        suppressMessages(cmdstanr::set_cmdstan_path(path))
        do.call(model$sample, args)$output_files()
      }, args = list(model = model, args = c(args, chain_ids = chain),
                     path = path))
    })
    fit <- cmdstanr::as_cmdstan_fit(unlist(files), check_diagnostics = FALSE)
    summary <- fit$diagnostic_summary(quiet = TRUE)
    divergent <- sum(summary$num_divergent)
    treedepth <- sum(summary$num_max_treedepth)
    written <- unique(sub("\\[.*$", "", fit$metadata()$variables))
    draws <- fit$draws(variables = intersect(variables, written))
  }

  # Only the variables uSDT reads are kept, whatever their size in this fit.
  present <- posterior::variables(draws)
  keep <- present[sub("\\[.*$", "", present) %in% variables]
  list(fit = fit, draws = posterior::subset_draws(draws, variable = keep),
       divergent = divergent, treedepth = treedepth)
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
  done <- cli::make_ansi_style("#6CB4EE")("■")
  list(complete = done, current = done,
       incomplete = cli::make_ansi_style("#9AA0A6")("─"))
}
