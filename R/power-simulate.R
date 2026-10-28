# power-simulate.R
# Power of the uSDT hypotheses by simulating and analysing whole studies
# Author: Ricardo Rey-Sáez
# Last modified: 28-10-2026

# Internal functions

# This function simulates one study: N subjects from one scenario of a
# population, each with the given signal and noise trials per task. Each
# subject draws its two d' from the bivariate normal of the population and its
# criteria from theirs. A task split at the median draws one latent score per
# trial, N(d'/2, sd_signal^2) for a signal trial and N(-d'/2, 1) for a noise
# trial, and meyen_split() turns the higher half of each subject's scores into
# signal responses, as with response times. Any other task draws its counts
# from the binomial distribution its criterion gives. The result is a
# usdt_data object built as usdt_data_tasks() builds one, so hsdt() fits the
# same model it fits to real data, criterion check included.
.simulate_usdt_data <- function(pop, N, cells, split, labels) {
  w <- matrix(stats::rnorm(2L * N), N)
  d <- cbind(direct = pop$mu_D + pop$sigma_D * w[, 1L],
             indirect = pop$mu_I + pop$sigma_I *
               (pop$rho * w[, 1L] + sqrt(1 - pop$rho^2) * w[, 2L]))
  u <- matrix(stats::rnorm(2L * N), N)
  r <- pop$cor_criterion
  criterion <- cbind(
    direct = pop$criterion_D + pop$sd_criterion_D * u[, 1L],
    indirect = pop$criterion_I + pop$sd_criterion_I *
      (r * u[, 1L] + sqrt(1 - r^2) * u[, 2L]))
  sigma <- c(direct = pop$sd_signal_D, indirect = pop$sd_signal_I)
  subjects <- formatC(seq_len(N), width = nchar(N), flag = "0")

  counts <- lapply(c("direct", "indirect"), function(k) {
    ns <- cells$signal[[k]]
    nn <- cells$noise[[k]]
    if (split[[k]]) {
      signal <- d[, k] / 2 + sigma[[k]] * matrix(stats::rnorm(N * ns), N)
      noise <- -d[, k] / 2 + matrix(stats::rnorm(N * nn), N)
      scores <- cbind(signal, noise)
      responses <- matrix(meyen_split(as.vector(scores),
                                      by = rep(subjects, times = ns + nn),
                                      signal = "slower"), N)
      hits <- rowSums(responses[, seq_len(ns), drop = FALSE])
      false_alarms <- rowSums(responses[, ns + seq_len(nn), drop = FALSE])
    } else {
      hits <- stats::rbinom(N, ns, stats::pnorm((d[, k] / 2 - criterion[, k]) /
                                                  sigma[[k]]))
      false_alarms <- stats::rbinom(N, nn, stats::pnorm(-d[, k] / 2 -
                                                          criterion[, k]))
    }
    list(hits = hits, false_alarms = false_alarms, signal = ns, noise = nn)
  })
  names(counts) <- c("direct", "indirect")

  # The aggregated rows follow the layout of usdt_data_tasks(), in deviation
  # coding.
  agg <- do.call(rbind, lapply(c("direct", "indirect"), function(k) {
    x <- counts[[k]]
    data.frame(subj = rep(subjects, each = 2L),
               task = if (k == "direct") "D" else "I",
               sig = rep(c(FALSE, TRUE), N),
               y = as.integer(as.vector(rbind(x$false_alarms, x$hits))),
               n = as.integer(rep(c(x$noise, x$signal), N)),
               stringsAsFactors = FALSE)
  }))
  agg$task <- factor(agg$task, levels = c("D", "I"))
  agg$cond <- ifelse(agg$sig, 0.5, -0.5)
  agg$c_D <- as.integer(agg$task == "D")
  agg$c_I <- as.integer(agg$task == "I")
  agg$d_D <- agg$cond * (agg$task == "D")
  agg$d_I <- agg$cond * (agg$task == "I")
  agg <- agg[order(agg$subj, agg$task, agg$sig), ]
  rownames(agg) <- NULL

  # A split task has its criterion fixed only if the split leaves it at zero,
  # as usdt_data_tasks() decides.
  check <- .criterion_check(agg, "deviation")
  tasks <- lapply(c("direct", "indirect"), function(k) {
    total <- cells$signal[[k]] + cells$noise[[k]]
    list(dichotomized = split[[k]], subject = "subject",
         condition_col = "condition", response_col = "response",
         condition_levels = c(signal = "signal", noise = "noise"),
         response_levels = c(signal = 1, noise = 0),
         dic = list(odd = if (total %% 2L == 1L) N else 0L,
                    p = rep(1 / 2, N)))
  })
  names(tasks) <- c("direct", "indirect")
  structure(list(
    agg = agg,
    meta = list(
      input = "simulated study", entry = "usdt_power()", coding = "deviation",
      ties = c(direct = "noise", indirect = "noise"), labels = labels,
      tasks = tasks,
      criterion_zero = c(direct = split[["direct"]] && check$direct$negligible,
                         indirect = split[["indirect"]] &&
                           check$indirect$negligible),
      criterion = check, n_subj = N, n_both = N, n_only = 0L,
      n_trials = sum(agg$n), n_rows = nrow(agg),
      descriptives = .descriptives(agg, "deviation"))),
    class = "usdt_data")
}

# This function fits one simulated study. It only calls hsdt(), quietly, so
# the tests can replace it with a fake that returns a fit without fitting.
.power_fit <- function(data, level) {
  suppressWarnings(suppressMessages(hsdt(data, level = level)))
}

# Columns kept for every replicate: the status of the fit, the estimate,
# standard error and interval of each quantity, and one decision per
# hypothesis.
.replicate_quantities <- c("diff", "rho", "intercept", "slope")

# This function simulates, fits and tests one study, with its own random
# number stream. It returns a numeric vector: whether the fit worked, whether
# the correlation and latent regression are estimable, the estimate, standard
# error and (1 - alpha) interval of each quantity from usdt_tests(), and the
# decision on each hypothesis, as usdt_tests() decides a point null at zero
# and usdt_equivalence() any other statement.
.power_replicate <- function(stream, pop, N, cells, split, labels, tests,
                             alpha) {
  env <- globalenv()
  had_seed <- exists(".Random.seed", envir = env, inherits = FALSE)
  if (had_seed) old <- get(".Random.seed", envir = env, inherits = FALSE)
  on.exit(if (had_seed) assign(".Random.seed", old, envir = env) else
    rm(".Random.seed", envir = env))
  assign(".Random.seed", stream, envir = env)

  q <- .replicate_quantities
  out <- c(fitted = 0, estimable = 0,
           stats::setNames(rep(NA_real_, 4L * length(q)),
                           paste0(rep(c("est_", "se_", "lo_", "hi_"),
                                      each = length(q)), q)),
           stats::setNames(rep(0, length(tests)),
                           paste0("reject_", seq_along(tests))))
  data <- .simulate_usdt_data(pop, N, cells, split, labels)
  m <- tryCatch(.power_fit(data, 1 - alpha), error = function(e) NULL)
  if (is.null(m)) return(out)

  # The rows of usdt_tests() give every quantity its interval.
  rows <- m$tests
  at <- match(c("d'(indirect) - d'(direct)", "correlation", "intercept",
                "slope"), rows$term)
  ok <- rows$status[at] == "ok"
  out[["fitted"]] <- 1
  out[["estimable"]] <- as.numeric(isTRUE(m$pars$joint_ok))
  out[paste0("est_", q)] <- rows$estimate[at]
  out[paste0("se_", q)] <- ifelse(ok, rows$se[at], NA_real_)
  out[paste0("lo_", q)] <- ifelse(ok, rows$conf.low[at], NA_real_)
  out[paste0("hi_", q)] <- ifelse(ok, rows$conf.high[at], NA_real_)

  sampling <- NULL
  for (j in seq_along(tests)) {
    h <- tests[[j]]
    if (.power_test(h) == "covariance" ||
        (h$op == "=" && h$value == 0 && h$quantity %in% c("diff", "intercept"))) {
      p <- rows$p.value[at[match(h$quantity, q)]]
      reject <- isTRUE(p < alpha)
    } else {
      if (is.null(sampling)) sampling <- .equiv_wald(m$pars)
      row <- .equiv_one(h, sampling[[h$quantity]], alpha)
      reject <- identical(row$decision, "reject H0")
    }
    out[[paste0("reject_", j)]] <- as.numeric(reject)
  }
  out
}

# This function gives one random number stream per replicate, so a run gives
# the same studies whatever the number of cores. The streams follow from
# `seed`, or from one draw of the session's own stream when it is NULL.
.power_streams <- function(n, seed) {
  env <- globalenv()
  had_seed <- exists(".Random.seed", envir = env, inherits = FALSE)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  if (had_seed) old <- get(".Random.seed", envir = env, inherits = FALSE)
  on.exit(if (had_seed) assign(".Random.seed", old, envir = env) else
    rm(".Random.seed", envir = env))
  set.seed(seed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion",
           sample.kind = "Rejection")
  streams <- vector("list", n)
  streams[[1L]] <- get(".Random.seed", envir = env)
  for (i in seq_len(n - 1L)) {
    streams[[i + 1L]] <- parallel::nextRNGStream(streams[[i]])
  }
  streams
}

# This function starts the worker processes and loads uSDT in each: the
# installed package, or the source tree when the session runs a development
# copy loaded with pkgload.
.power_cluster <- function(ncores) {
  cl <- parallel::makePSOCKcluster(as.integer(ncores))
  path <- getNamespaceInfo(asNamespace("uSDT"), "path")
  source_tree <- dir.exists(file.path(path, "R")) &&
    !file.exists(file.path(path, "R", "uSDT.rdb"))
  ready <- tryCatch(unlist(parallel::clusterCall(cl, function(path, source_tree) {
    if (source_tree) {
      load_all <- getExportedValue("pkgload", "load_all")
      suppressMessages(load_all(path, quiet = TRUE, export_all = TRUE,
                                helpers = FALSE, attach_testthat = FALSE))
    } else {
      loadNamespace("uSDT", lib.loc = dirname(path))
    }
    TRUE
  }, path, source_tree)), error = function(e) FALSE)
  if (!all(ready)) {
    parallel::stopCluster(cl)
    .usdt_stop("the parallel workers could not load uSDT.")
  }
  cl
}

# This function runs `nsim` replicates of one scenario and design, on the
# workers of `cl` when there are any, and returns one row per replicate. The
# progress bar `pb` counts replicates across every design of a run.
.power_simulate <- function(streams, pop, N, cells, split, labels, tests,
                            alpha, cl = NULL, pb = NULL, done = 0L) {
  run <- function(s) {
    .power_replicate(s, pop, N, cells, split, labels, tests, alpha)
  }
  batch <- if (is.null(cl)) 10L else 10L * length(cl)
  starts <- seq.int(1L, length(streams), by = batch)
  rows <- vector("list", length(starts))
  for (b in seq_along(starts)) {
    idx <- starts[b]:min(starts[b] + batch - 1L, length(streams))
    rows[[b]] <- if (is.null(cl)) {
      lapply(streams[idx], run)
    } else {
      parallel::parLapply(cl, streams[idx], .power_replicate, pop = pop,
                          N = N, cells = cells, split = split,
                          labels = labels, tests = tests, alpha = alpha)
    }
    if (!is.null(pb)) utils::setTxtProgressBar(pb, done + max(idx))
  }
  out <- do.call(rbind, unlist(rows, recursive = FALSE))
  rownames(out) <- NULL
  out
}

# This function summarises the replicates of one scenario and design: the
# share that rejects each hypothesis, with its Monte Carlo error and Wilson
# interval, and the share in which the analysis could estimate its quantity.
.power_from_replicates <- function(reps, tests, level = 0.95) {
  n <- nrow(reps)
  z <- stats::qnorm(1 - (1 - level) / 2)
  wilson <- function(p) {
    centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
    half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
    c(centre - half, centre + half)
  }
  decisions <- reps[, paste0("reject_", seq_along(tests)), drop = FALSE] == 1
  rows <- lapply(seq_along(tests), function(j) {
    p <- mean(decisions[, j])
    estimable <- if (tests[[j]]$quantity == "diff") reps[, "fitted"] else
      reps[, "estimable"]
    ci <- wilson(p)
    data.frame(power = p, mc_se = sqrt(p * (1 - p) / n), power.low = ci[1L],
               power.high = ci[2L], estimable = mean(estimable == 1))
  })
  joint <- mean(rowSums(decisions) == length(tests))
  only_diff <- all(vapply(tests, function(h) h$quantity == "diff", TRUE))
  list(rows = do.call(rbind, rows), joint = joint,
       joint_mc_se = sqrt(joint * (1 - joint) / n), joint_ci = wilson(joint),
       joint_estimable = mean(reps[, if (only_diff) "fitted" else
         "estimable"] == 1),
       failures = sum(reps[, "fitted"] == 0))
}

# This function gives the rows of a population the approximation cannot
# describe (unequal variances): the population value of each quantity and
# whether it lies in H0, without standard errors.
.power_truth <- function(pop, tests) {
  values <- c(diff = pop$diff, rho = pop$rho, intercept = pop$intercept,
              slope = pop$slope)
  rows <- do.call(rbind, lapply(tests, function(h) {
    theta <- values[[h$quantity]]
    data.frame(true = theta, se = NA_real_, power = NA_real_,
               estimable = NA_real_,
               in_null = .power_in_null(h$op, h$value, theta))
  }))
  list(rows = rows, joint = NA_real_, joint_estimable = NA_real_)
}

# This function puts the simulated power of a design in place of the
# approximate one, keeping the population values and expected standard errors.
.power_merge <- function(res, sim) {
  res$rows$power <- sim$rows$power
  res$rows$estimable <- sim$rows$estimable
  res$rows$mc_se <- sim$rows$mc_se
  res$rows$power.low <- sim$rows$power.low
  res$rows$power.high <- sim$rows$power.high
  res$joint <- sim$joint
  res$joint_mc_se <- sim$joint_mc_se
  res$joint_ci <- sim$joint_ci
  res$joint_estimable <- sim$joint_estimable
  res$failures <- sim$failures
  res
}
