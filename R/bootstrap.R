# bootstrap.R
# This script runs a parametric bootstrap for fitted uSDT models.
# Author: Ricardo Rey-Sáez
# Last modified: 08-09-2026

# Public functions

# Fewer usable replicates than this cannot support a two-sided interval.
.boot_min <- 500L

#' Bootstrap intervals for a hierarchical SDT model
#'
#' Simulates many datasets from the fitted model, refits the model to each one,
#' and builds the intervals of the three hypotheses from the results. This is
#' useful when the ordinary intervals are unavailable or hard to trust, which
#' happens when the model reaches a boundary. The work is done by
#' `lme4::bootMer()`.
#'
#' @param object An `hsdt` object from [hsdt()].
#' @param nsim Number of usable replicates to reach. It must be at least 500.
#' @param ncores Number of cores to use. Values above one run the replicates in
#'   parallel through the `parallel` package, which comes with R. The temporary
#'   cluster behaves the same way on Windows, macOS and Linux, and it closes
#'   when the bootstrap ends.
#' @param max_attempts Largest number of replicates to fit. The default allows
#'   two attempts for every usable replicate requested.
#' @param seed Seed for the simulated datasets, so the result can be
#'   reproduced.
#' @param progress Show a progress bar. It appears by default in interactive
#'   sessions.
#' @param level Confidence level.
#' @param type Type of interval. `"perc"` takes the percentiles of the
#'   replicates, `"norm"` builds a normal interval around the bias-corrected
#'   estimate, and `"basic"` reflects the percentiles around the estimate. The
#'   last two work on the Fisher-z scale for the correlation, which keeps their
#'   limits inside its range. They need a finite centre on that scale, so a
#'   correlation that sits on the boundary reports them as missing. `"perc"`
#'   stays available in that case.
#'
#' @return The `hsdt` object, with the interval columns of its `tests` table
#'   replaced by the bootstrap results. The new `boot` element holds the
#'   replicates of the three hypotheses in `t`, the sensitivity variances in
#'   `variance`, the average task parameters in `population`, and the estimates
#'   of every subject in `subjects`. It also holds the counts and diagnostics
#'   of the run.
#'
#' @details
#' The function drops a replicate only when the model fails to fit or fails to
#' converge. It keeps singular and boundary replicates, because they are the
#' answer the model gives for a difficult dataset, and removing them would make
#' the intervals narrower than they should be. `boot$retained` reports how many
#' there were. The run continues until it reaches `nsim` usable replicates or
#' `max_attempts` fitted samples.
#'
#' An incomplete run still returns the object, with every attempt and its
#' diagnostics in `boot`, and it gives a warning. Bootstrap summaries replace
#' the original intervals only from 500 usable replicates onwards.
#'
#' The point estimates do not change. A bootstrap describes how much an
#' estimate would vary from sample to sample, and the estimate itself remains
#' the one the model produced. The bootstrap p-value compares the fitted
#' estimate in absolute value with the centred distribution of the replicates.
#' The count adds one to the numerator and the denominator, so a finite
#' simulation never returns a p-value of zero.
#'
#' The `population` and `subjects` components keep four parameters from every
#' attempted replicate, the two criterion intercepts `c_D` and `c_I` and the
#' two sensitivities `d_D` and `d_I`. A criterion that the model fixed is
#' stored as zero, and the values of a subject combine the refitted average
#' with that subject's own departure from it. Their first dimension follows
#' `boot$ok`, so the same usable replicates can be selected again.
#'
#' @seealso [hsdt()], [usdt_tests()]
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' df <- usdt_simulate(n_subj = 40, n_trials = 100)
#' d  <- usdt_data_long(df, task_col = "task",
#'                      task_levels      = c(direct = "D", indirect = "I"),
#'                      subject_col      = "subj",
#'                      condition_col    = "cond",
#'                      condition_levels = c(signal = 1, noise = 0),
#'                      response_col     = "response",
#'                      response_levels  = c(signal = 1, noise = 0))
#' m <- hsdt(d)
#' b <- usdt_boot(m, nsim = 500)
#' summary(b)
#' }
#'
#' @export
usdt_boot <- function(object, nsim = 1000, ncores = 1L,
                      max_attempts = 2 * nsim, seed = NULL,
                      progress = interactive(), level = 0.95,
                      type = c("perc", "norm", "basic")) {

  # The function requires a fitted uSDT model.
  if (!inherits(object, "hsdt")) {
    .usdt_stop("`object` must come from hsdt(), not a plain ",
               class(object)[1L], ".")
  }
  .check_scalar_number(nsim, "nsim", lower = .boot_min, whole = TRUE)
  .check_scalar_number(ncores, "ncores", lower = 1, whole = TRUE)
  .check_scalar_number(max_attempts, "max_attempts", lower = nsim, whole = TRUE)
  .check_confidence_level(level)
  if (!is.null(seed)) {
    .check_scalar_number(seed, "seed", lower = 0, whole = TRUE)
  }
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    .usdt_stop("`progress` must be `TRUE` or `FALSE`.")
  }
  type <- match.arg(type)
  nsim <- as.integer(nsim)
  max_attempts <- as.integer(max_attempts)

  # Each refitted model returns the reported and task-specific estimates.
  stat <- .boot_statistic(object$fit, "d_D", "d_I")
  layout <- attr(stat, "layout")
  t0 <- stat(object$fit)

  # One cluster serves every batch on every operating system.
  par_args <- list()
  if (ncores > 1L) {
    cl <- parallel::makePSOCKcluster(as.integer(ncores))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    ready <- unlist(parallel::clusterCall(
      cl, requireNamespace, package = "lme4", quietly = TRUE
    ))
    if (!all(ready)) {
      .usdt_stop("the parallel workers could not load lme4.")
    }
    par_args <- list(parallel = "snow", ncpus = as.integer(ncores), cl = cl)
  }

  # The seed starts one continuous sequence across all batches.
  if (!is.null(seed)) {
    set.seed(seed)
  } else if (!exists(".Random.seed", envir = .GlobalEnv)) {
    stats::runif(1)
  }

  # The progress bar counts usable samples.
  pb <- NULL
  if (progress) {
    pb <- utils::txtProgressBar(min = 0, max = nsim, style = 3)
    on.exit(if (!is.null(pb)) close(pb), add = TRUE)
  }

  # Small batches keep the progress visible and reuse the same workers.
  t <- matrix(numeric(0), nrow = 0L, ncol = length(t0),
              dimnames = list(NULL, names(t0)))
  ok <- logical(0)
  messages <- list()
  batch_error <- NULL
  attempted <- 0L
  usable <- 0L
  batch_size <- min(100L, nsim)

  while (usable < nsim && attempted < max_attempts) {
    batch_n <- min(batch_size, nsim - usable, max_attempts - attempted)
    bt <- tryCatch(
      suppressWarnings(do.call(
        lme4::bootMer,
        c(list(x = object$fit, FUN = stat, nsim = batch_n,
               type = "parametric", use.u = FALSE, seed = NULL), par_args)
      )),
      error = function(e) e
    )
    if (inherits(bt, "error")) {
      batch_error <- conditionMessage(bt)
      break
    }

    filter <- .boot_filter(bt$t)
    t <- rbind(t, bt$t)
    ok <- c(ok, filter$ok)
    messages[[length(messages) + 1L]] <- attr(bt, "boot.all.msgs")
    attempted <- nrow(t)
    usable <- sum(ok)
    if (!is.null(pb)) utils::setTxtProgressBar(pb, usable)
  }
  if (!is.null(pb)) {
    close(pb)
    pb <- NULL
  }

  # Every attempted sample remains available in the returned object.
  filter <- .boot_filter(t)
  usable <- sum(filter$ok)
  complete <- usable >= nsim
  summary_available <- usable >= .boot_min
  if (summary_available) {
    object$tests <- .boot_tests(object$tests,
                                t[filter$ok, , drop = FALSE], level, type)
    object$level <- level
  }
  effects <- .boot_effects(t0, t, layout)
  object$boot <- list(
    t0 = effects$summary$t0, t = effects$summary$t,
    variance = effects$variance,
    population = effects$population, subjects = effects$subjects,
    ok = filter$ok, requested = nsim,
    attempted = attempted, usable = usable, max_attempts = max_attempts,
    complete = complete, summary_available = summary_available,
    level = level, type = type,
    failures = c(non_finite = sum(!filter$finite),
                 non_converged = sum(!filter$converged)),
    retained = c(singular = sum(filter$ok & filter$singular),
                 boundary = sum(filter$ok & filter$boundary)),
    messages = messages, error = batch_error
  )

  # A warning preserves the result while making weak output explicit.
  if (!complete) {
    .usdt_warn("bootstrap returned ", usable, " of ", nsim,
               " requested usable replicates after ", attempted, " attempts",
               if (!is.null(batch_error)) paste0(". The last batch failed: ",
                                                  batch_error) else ".",
               if (summary_available)
                 " Hypothesis summaries use the available replicates."
               else
                 paste0(" Fewer than ", .boot_min, " were usable, so the ",
                        "original hypothesis results remain."),
               " Inspect `object$boot` for all attempts and diagnostics.")
  } else if (usable / attempted < 0.8) {
    .usdt_warn("bootstrap reached ", nsim, " usable replicates after ",
               attempted, " attempts, but more than 20% were discarded. ",
               "Inspect `object$boot$failures` before interpreting the results.")
  }
  object
}

# Internal functions

# This function creates the calculation used by each bootstrap sample.
.boot_statistic <- function(fit, direct, indirect) {
  blk <- .locate_block(fit, direct, indirect)
  fixed <- names(lme4::fixef(fit))
  subjects <- rownames(lme4::ranef(fit, condVar = FALSE)[["subj"]])
  parameters <- c("c_D", "c_I", direct, indirect)
  summary_names <- c("diff", "rho", "intercept", "slope",
                     "converged", "singular", "boundary")
  variance_names <- c("s2_D", "s2_I")
  variance_output_names <- paste0("variance.", variance_names)
  population_names <- paste0("population.", parameters)
  subject_names <- unlist(lapply(seq_along(parameters), function(j) {
    paste0("subject.", parameters[j], ".", seq_along(subjects))
  }), use.names = FALSE)
  output_names <- c(summary_names, variance_output_names,
                    population_names, subject_names)
  stat <- function(f, block_index, block_size, direct_random,
                   indirect_random, direct_fixed, indirect_fixed,
                   subjects, parameters, output_names) {
    theta <- lme4::getME(f, "theta")
    beta <- lme4::fixef(f)
    L <- matrix(0, block_size, block_size)
    L[lower.tri(L, diag = TRUE)] <- theta[block_index]
    S <- tcrossprod(L)

    gamma_D <- unname(beta[direct_fixed])
    gamma_I <- unname(beta[indirect_fixed])
    s2_D <- S[direct_random, direct_random]
    s2_I <- S[indirect_random, indirect_random]
    s_DI <- S[direct_random, indirect_random]
    slope <- s_DI / s2_D
    # The worker repeats the correlation limit in its own environment.
    rho <- max(-1, min(1, s_DI / sqrt(s2_D * s2_I)))

    opt <- f@optinfo$conv$opt
    opt_ok <- !length(opt) || (is.numeric(opt) && all(opt == 0))
    messages <- f@optinfo$conv$lme4$messages
    if (length(messages)) {
      messages <- messages[!grepl("boundary.*singular", messages,
                                  ignore.case = TRUE)]
    }

    population <- stats::setNames(rep.int(0, length(parameters)), parameters)
    in_model <- intersect(parameters, names(beta))
    population[in_model] <- beta[in_model]

    random <- lme4::ranef(f, condVar = FALSE)[["subj"]]
    subject_row <- match(subjects, rownames(random))
    individual <- matrix(0, nrow = length(subjects), ncol = length(parameters),
                         dimnames = list(subjects, parameters))
    random_parameters <- intersect(parameters, colnames(random))
    individual[, random_parameters] <-
      as.matrix(random[subject_row, random_parameters, drop = FALSE])
    individual <- sweep(individual, 2L, population, "+")

    finite <- all(is.finite(c(gamma_D, gamma_I, S, slope, rho,
                              population, individual)))
    boundary <- TRUE
    if (finite) {
      values <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
      boundary <- min(s2_D, s2_I) < 1e-8 || min(values) < 1e-8 ||
        abs(rho) > 1 - 1e-6
    }

    stats::setNames(
      c(gamma_D - gamma_I, rho, gamma_I - slope * gamma_D, slope,
        as.numeric(opt_ok && !length(messages)),
        as.numeric(lme4::isSingular(f, tol = 1e-4)),
        as.numeric(boundary), s2_D, s2_I,
        population, as.vector(individual)),
      output_names
    )
  }
  args <- formals(stat)
  args$block_index <- blk$idx
  args$block_size <- blk$q
  args$direct_random <- blk$iD
  args$indirect_random <- blk$iI
  args$direct_fixed <- match(direct, fixed)
  args$indirect_fixed <- match(indirect, fixed)
  args$subjects <- subjects
  args$parameters <- parameters
  args$output_names <- output_names
  formals(stat) <- args
  environment(stat) <- baseenv()
  attr(stat, "layout") <- list(
    summary = seq_along(summary_names),
    variance = length(summary_names) + seq_along(variance_names),
    population = length(summary_names) + length(variance_names) +
      seq_along(parameters),
    subjects = length(summary_names) + length(variance_names) +
      length(parameters) +
      seq_len(length(subjects) * length(parameters)),
    subject_ids = subjects,
    parameters = parameters,
    variance_names = variance_names
  )
  stat
}

# This function separates compact bootstrap arrays from the worker output.
.boot_effects <- function(t0, t, layout) {

  # This block keeps the hypothesis estimates and diagnostics.
  summary <- list(
    t0 = t0[layout$summary],
    t = t[, layout$summary, drop = FALSE]
  )

  # This block keeps both sensitivity variances.
  variance <- list(
    t0 = stats::setNames(unname(t0[layout$variance]), layout$variance_names),
    t = t[, layout$variance, drop = FALSE]
  )
  colnames(variance$t) <- layout$variance_names

  # This block keeps one fixed value for each task parameter.
  population <- list(
    t0 = stats::setNames(unname(t0[layout$population]), layout$parameters),
    t = t[, layout$population, drop = FALSE]
  )
  colnames(population$t) <- layout$parameters

  # This block keeps the subject estimates from each replicate.
  subjects <- list(
    t0 = matrix(unname(t0[layout$subjects]),
                nrow = length(layout$subject_ids),
                dimnames = list(layout$subject_ids, layout$parameters)),
    t = array(unname(t[, layout$subjects, drop = FALSE]),
              dim = c(nrow(t), length(layout$subject_ids),
                      length(layout$parameters)),
              dimnames = list(NULL, layout$subject_ids, layout$parameters))
  )
  list(summary = summary, variance = variance,
       population = population, subjects = subjects)
}

# This function separates failed samples from fitted samples.
.boot_filter <- function(t) {
  estimates <- c("diff", "rho", "intercept", "slope")
  if (!nrow(t)) {
    empty <- logical(0)
    return(list(ok = empty, finite = empty, converged = empty,
                singular = empty, boundary = empty))
  }
  finite <- rowSums(!is.finite(t[, estimates, drop = FALSE])) == 0
  converged <- !is.na(t[, "converged"]) & t[, "converged"] == 1
  list(ok = finite & converged,
       finite = finite, converged = converged,
       singular = !is.na(t[, "singular"]) & t[, "singular"] == 1,
       boundary = !is.na(t[, "boundary"]) & t[, "boundary"] == 1)
}

# This function adds bootstrap results to the hypothesis table.
.boot_tests <- function(tests, t, level, type) {
  map <- c("d'(direct) - d'(indirect)" = "diff", correlation = "rho",
           intercept = "intercept", slope = "slope")

  for (i in seq_len(nrow(tests))) {
    col <- map[[tests$term[i]]]
    v   <- t[, col]
    est <- tests$estimate[i]

    # The centred replicates stand in for the distribution under the null.
    centred <- v - est

    # A correlation needs the Fisher-z scale for the location intervals.
    ci <- .boot_ci(v, est, level, type,
                   link = if (col == "rho") .fisher_link else NULL)

    tests$se[i]        <- stats::sd(v)
    tests$statistic[i] <- NA_real_
    tests$p.value[i]   <- (sum(abs(centred) >= abs(est)) + 1) / (length(v) + 1)
    tests$conf.low[i]  <- ci[1L]
    tests$conf.high[i] <- ci[2L]
    tests$ci_method[i] <- paste0("bootstrap (", type, ")")
    if (tests$status[i] == "not estimable") {
      tests$status[i] <- "ok"
      tests$reason[i] <- NA_character_
    }
  }
  tests
}
