# reliability.R
# Estimate task-level and subject-level reliability for uSDT models
# Author: Ricardo Rey-Sáez
# Last modified: 03-10-2026

#' Reliability of direct and indirect task measures
#'
#' Estimates the reliability of sensitivity (\eqn{d'}) by separating true
#' variance across participants from sampling noise caused by finite trial counts.
#' Values close to 1 indicate that the measure reliably separates participants,
#' whereas values near 0 indicate that observed differences are mostly
#' measurement error.
#'
#' @param object An `hsdt` object from [hsdt()], frequentist or Bayesian. A
#'   frequentist fit may also contain bootstrap results from [usdt_boot()].
#'
#' @return An object of class `usdt_reliability` containing:
#' * `$tasks`: Overall reliability and variance components for each task.
#' * `$subjects`: Participant-level \eqn{d'}, error variances, and individual
#'   reliabilities.
#' * Bootstrap intervals for both components when available in `object`.
#'
#' For a Bayesian fit, both tables hold posterior means and add the columns
#' `sd`, `conf.low` and `conf.high`, and `$posterior$tasks` keeps the draws of
#' each task's reliability.
#'
#' @details
#' For participant \eqn{i} in task \eqn{j}, reliability is defined as:
#' \deqn{\frac{\tau_j^2}{\tau_j^2 + v_{ij}}}
#' where \eqn{\tau_j^2} is the true variance in sensitivity across participants
#' from the model's random effects, and \eqn{v_{ij}} is the squared standard
#' error of \eqn{d'} for that participant. This error variance reflects how
#' precisely their trials determine sensitivity, accounting for trial count,
#' performance level on the probit curve, and uncertainty in the criterion.
#'
#' This variance is calculated using the large-sample Fisher information formula
#' from Gourevitch and Galanter (1967). Unlike [sdt_moments()], which evaluates
#' that formula at raw empirical proportions (`var_gg`), `usdt_reliability()`
#' evaluates it at the response probabilities predicted by the fitted
#' hierarchical model.
#'
#' Overall task reliability averages \eqn{v_{ij}} across participants,
#' representing the expected proportion of true variance for a participant
#' drawn at random from the sample.
#'
#' When `object` includes bootstrap replicates from [usdt_boot()], confidence
#' intervals for reliability are computed automatically across all retained
#' samples.
#'
#' For a Bayesian fit, reliability is computed draw by draw: each draw gives
#' \eqn{\tau_j^2} and, at that draw's subject sensitivities and criteria, each
#' \eqn{v_{ij}}. The tables report posterior means, posterior standard
#' deviations and central credible intervals. Under unequal variances the
#' information of the signal trials is scaled by the signal standard deviation
#' \eqn{\sigma_s}, so that \eqn{v_{ij} = \sigma_s^2 / w_{\mathrm{signal}} +
#' 1 / w_{\mathrm{noise}}} with an estimated criterion, which reduces to the
#' Gourevitch and Galanter formula when \eqn{\sigma_s = 1}.
#'
#' @references
#' Gourevitch, V., & Galanter, E. (1967). A significance test for one parameter
#' isosensitivity functions. \emph{Psychometrika}, 32(1), 25--33.
#' \doi{10.1007/BF02289402}
#'
#' @seealso [hsdt()], [usdt_boot()], [sdt_moments()]
#'
#' @examples
#' # Contextual cuing data from Vadillo et al. (2025)
#' d <- usdt_data_tasks(
#'   direct   = vadillo_awareness,
#'   indirect = vadillo_cuing,
#'   subject_col      = "subj",
#'   condition_col    = "condition",
#'   condition_levels = c(signal = "old", noise = "new"),
#'   response_col     = list(direct = "judged.old", indirect = "rt"),
#'   response_levels  = list(direct   = c(signal = 1, noise = 0),
#'                           indirect = c(signal = "faster", noise = "slower")),
#'   dichotomize      = list(direct = FALSE, indirect = TRUE)
#' )
#'
#' r <- usdt_reliability(hsdt(d))
#' r
#'
#' # Participant-level estimates (one row per subject and task)
#' head(r$subjects)
#'
#' @export
usdt_reliability <- function(object) {

  # Reliability needs the fitted subject variation.
  if (!inherits(object, "hsdt")) {
    .usdt_stop("`object` must come from hsdt(), not a plain ",
               class(object)[1L], ".\n  Reliability needs the fitted ",
               "between-subject variance and the information each subject ",
               "contributes.")
  }
  if (.is_bayes(object)) return(.reliability_bayes(object, match.call()))

  agg     <- object$data$agg
  labels  <- object$data$meta$labels
  effects <- .task_effects(object$fit, unique(as.character(agg$subj)))$subjects

  results <- lapply(c("D", "I"), function(task) {
    slope   <- paste0("d_", task)
    columns <- c(intersect(paste0("c_", task), object$design$criteria), slope)
    tau2    <- object$pars$est[[if (task == "D") "s2_D" else "s2_I"]]
    rows    <- agg[agg$task == task, ]
    who     <- unique(as.character(rows$subj))

    # Each subject supplies the cells that inform their own sensitivity.
    values <- vapply(who, function(s) {
      cells <- rows[as.character(rows$subj) == s, , drop = FALSE]
      coef  <- effects[s, ][columns]
      eta   <- drop(as.matrix(cells[, columns, drop = FALSE]) %*% coef)
      c(dprime   = unname(coef[[slope]]),
        variance = .sensitivity_variance(cells, eta, length(columns) == 2L))
    }, c(dprime = 0, variance = 0))

    label <- labels[[if (task == "D") "direct" else "indirect"]]
    list(
      task = data.frame(task = label, subjects = length(who), tau2 = tau2,
                        mean_variance = mean(values["variance", ]),
                        reliability = tau2 / (tau2 + mean(values["variance", ])),
                        stringsAsFactors = FALSE),
      subjects = data.frame(task = label, subj = who,
                            dprime = unname(values["dprime", ]),
                            variance = unname(values["variance", ]),
                            reliability = tau2 / (tau2 + values["variance", ]),
                            stringsAsFactors = FALSE)
    )
  })

  tasks    <- do.call(rbind, lapply(results, `[[`, "task"))
  subjects <- do.call(rbind, lapply(results, `[[`, "subjects"))
  rownames(tasks) <- rownames(subjects) <- NULL

  boot <- NULL
  if (!is.null(object$boot) && isTRUE(object$boot$summary_available)) {
    boot <- .reliability_boot(object, tasks, subjects)
    task_columns <- setdiff(names(boot$tasks), "task")
    subject_columns <- setdiff(names(boot$subjects), c("task", "subj"))
    tasks[task_columns] <- boot$tasks[
      match(tasks$task, boot$tasks$task), task_columns, drop = FALSE]
    subject_key <- paste(subjects$task, subjects$subj, sep = "\r")
    boot_key <- paste(boot$subjects$task, boot$subjects$subj, sep = "\r")
    subjects[subject_columns] <- boot$subjects[
      match(subject_key, boot_key), subject_columns, drop = FALSE]
  }

  # A sensitivity variance on the boundary leaves nothing to be reliable about.
  flat <- tasks$task[tasks$tau2 < 1e-8]
  if (length(flat)) {
    .usdt_warn("the between-subject variance of ", paste(flat, collapse = ", "),
               " is on the boundary, so its reliability is zero by ",
               "construction rather than by measurement.")
  }

  structure(list(tasks = tasks, subjects = subjects, boot = boot,
                 call = match.call()),
            class = "usdt_reliability")
}

# Internal functions

# This function gives the sampling variance of one subject's d' from the
# expected probit information of its signal and noise cells. `eta` holds the
# linear predictor of each cell, as a vector or as a matrix with one row per
# replicate or draw. With an estimated criterion the variance is
# sigma^2 / w_signal + 1 / w_noise; with the criterion fixed by the Meyen split
# it is (1 + sigma)^2 / (w_signal + w_noise). Equal variances (sigma = 1) give
# the formula of Gourevitch and Galanter (1967).
.sensitivity_variance <- function(cells, eta, free, sigma = 1) {
  eta <- matrix(eta, ncol = nrow(cells))
  p <- stats::pnorm(eta)
  w <- sweep(stats::dnorm(eta)^2 / (p * (1 - p)), 2L, cells$n, "*")
  w_signal <- w[, cells$sig]
  w_noise <- w[, !cells$sig]
  variance <- if (free) {
    sigma^2 / w_signal + 1 / w_noise
  } else {
    (1 + sigma)^2 / (w_signal + w_noise)
  }
  variance[!is.finite(variance) | variance <= 0] <- NA_real_
  variance
}

# This function computes reliability draw by draw for a Bayesian fit and
# summarises it with posterior means, SDs and central credible intervals.
.reliability_bayes <- function(object, call) {
  agg <- object$data$agg
  labels <- object$data$meta$labels
  free <- object$design$free_c
  a <- (1 - object$level) / 2
  draw <- function(name) {
    as.numeric(posterior::extract_variable(object$draws, name))
  }
  interval <- function(x) stats::quantile(x, c(a, 1 - a), names = FALSE)

  results <- lapply(1:2, function(j) {
    label <- labels[[c("direct", "indirect")[j]]]
    rows <- agg[agg$task == c("D", "I")[j], ]
    who <- unique(as.character(rows$subj))
    d <- .subject_draws(object, "d", j, who)
    criterion <- .subject_draws(object, "c", j, who)
    sigma <- if (object$design$unequal_variances) draw(sprintf("sigma_s[%d]", j)) else 1
    tau2 <- draw(sprintf("sigma_d[%d]", j))^2

    # Each draw places the subject's cells on its own probit curve.
    variance <- vapply(seq_along(who), function(k) {
      cells <- rows[as.character(rows$subj) == who[k], , drop = FALSE]
      eta <- vapply(cells$sig, function(signal) {
        if (signal) (d[, k] / 2 - criterion[, k]) / sigma else
          -d[, k] / 2 - criterion[, k]
      }, numeric(nrow(d)))
      .sensitivity_variance(cells, eta, free[j] == 1L, sigma)
    }, numeric(nrow(d)))

    reliability <- tau2 / (tau2 + variance)
    task_draws <- tau2 / (tau2 + rowMeans(variance))
    limits <- interval(task_draws)
    subject_limits <- apply(reliability, 2L, interval)
    list(
      task = data.frame(task = label, subjects = length(who), tau2 = mean(tau2),
                        mean_variance = mean(rowMeans(variance)),
                        reliability = mean(task_draws), sd = stats::sd(task_draws),
                        conf.low = limits[1L], conf.high = limits[2L],
                        stringsAsFactors = FALSE),
      subjects = data.frame(task = label, subj = who, dprime = colMeans(d),
                            variance = colMeans(variance),
                            reliability = colMeans(reliability),
                            sd = apply(reliability, 2L, stats::sd),
                            conf.low = subject_limits[1L, ],
                            conf.high = subject_limits[2L, ],
                            stringsAsFactors = FALSE),
      draws = task_draws)
  })

  tasks <- do.call(rbind, lapply(results, `[[`, "task"))
  subjects <- do.call(rbind, lapply(results, `[[`, "subjects"))
  rownames(tasks) <- rownames(subjects) <- NULL
  draws <- vapply(results, `[[`, numeric(length(results[[1L]]$draws)), "draws")
  colnames(draws) <- tasks$task
  structure(list(tasks = tasks, subjects = subjects, boot = NULL,
                 posterior = list(tasks = draws, level = object$level),
                 call = call),
            class = "usdt_reliability")
}

# This function recalculates reliability across bootstrap samples.
.reliability_boot <- function(object, tasks, subjects) {
  boot <- object$boot
  if (is.null(boot$variance) || is.null(boot$subjects)) {
    .usdt_stop("the bootstrap object lacks reliability parameters. Rerun ",
               "usdt_boot() with the current uSDT version.")
  }

  keep <- which(boot$ok)
  theta <- boot$variance$t[keep, , drop = FALSE]
  estimates <- boot$subjects$t[keep, , , drop = FALSE]
  n_boot <- length(keep)
  task_samples <- matrix(NA_real_, n_boot, 2L,
                         dimnames = list(NULL, tasks$task))
  subject_samples <- matrix(
    NA_real_, n_boot, nrow(subjects),
    dimnames = list(NULL, paste(subjects$task, subjects$subj, sep = "\r")))

  for (task in c("D", "I")) {
    label <- object$data$meta$labels[[if (task == "D") "direct" else "indirect"]]
    slope <- paste0("d_", task)
    columns <- c(intersect(paste0("c_", task), object$design$criteria), slope)
    rows <- object$data$agg[object$data$agg$task == task, ]
    who <- unique(as.character(rows$subj))
    estimate_rows <- match(who, dimnames(estimates)[[2L]])
    result_rows <- match(paste(label, who, sep = "\r"),
                         colnames(subject_samples))
    if (anyNA(estimate_rows) || anyNA(result_rows)) {
      .usdt_stop("the bootstrap subjects do not match the fitted data.")
    }

    local_variance <- matrix(NA_real_, n_boot, length(who))
    for (j in seq_along(who)) {
      cells <- rows[as.character(rows$subj) == who[j], , drop = FALSE]
      coefficient <- matrix(
        estimates[, estimate_rows[j], columns, drop = FALSE],
        nrow = n_boot, ncol = length(columns))
      eta <- coefficient %*% t(as.matrix(cells[, columns, drop = FALSE]))
      local_variance[, j] <- .sensitivity_variance(
        cells, eta, length(columns) == 2L)
    }

    tau2 <- theta[, paste0("s2_", task)]
    reliability <- tau2 / sweep(local_variance, 1L, tau2, "+")
    task_samples[, label] <- tau2 / (tau2 + rowMeans(local_variance))
    subject_samples[, result_rows] <- reliability
  }

  task_summary <- do.call(rbind, lapply(seq_len(ncol(task_samples)), function(i) {
    cbind(task = colnames(task_samples)[i],
          .reliability_boot_stats(task_samples[, i], tasks$reliability[i],
                                  object$level, boot$type))
  }))
  subject_summary <- do.call(rbind, lapply(seq_len(ncol(subject_samples)), function(i) {
    key <- strsplit(colnames(subject_samples)[i], "\r", fixed = TRUE)[[1L]]
    estimate <- subjects$reliability[
      subjects$task == key[1L] & subjects$subj == key[2L]]
    cbind(task = key[1L], subj = key[2L],
          .reliability_boot_stats(subject_samples[, i], estimate,
                                  object$level, boot$type))
  }))
  rownames(task_summary) <- rownames(subject_summary) <- NULL

  list(tasks = task_summary, subjects = subject_summary,
       samples = list(tasks = task_samples, subjects = subject_samples),
       usable = n_boot, level = object$level, type = boot$type)
}

# This function summarises one bootstrap distribution.
.reliability_boot_stats <- function(values, estimate, level, type) {
  values <- values[is.finite(values)]
  if (!length(values)) {
    return(data.frame(boot_median = NA_real_, boot_se = NA_real_,
                      conf.low = NA_real_, conf.high = NA_real_, boot_n = 0L))
  }
  link <- list(to = stats::qlogis, from = stats::plogis)
  ci <- .boot_ci(values, estimate, level, type, link = link)
  data.frame(boot_median = stats::median(values),
             boot_se = stats::sd(values),
             conf.low = ci[1L], conf.high = ci[2L],
             boot_n = length(values))
}

# Printed output

#' @param x A `usdt_reliability` object.
#' @param ... Ignored.
#' @rdname usdt_reliability
#' @export
summary.usdt_reliability <- function(object, ...) {

  # The table shows the group estimate and the subject distribution.
  cat(.rule("Reliability summary"), "\n\n")
  if (!is.null(object$posterior)) {
    cat(sprintf("  %-12s %8s %10s %8s  %-18s   %s\n",
                "Task", "Subjects", "Group mean", "SD",
                sprintf("%.0f%% CrI", 100 * object$posterior$level),
                "By-subject mean median [min, max]"))
    for (i in seq_len(nrow(object$tasks))) {
      task <- object$tasks$task[i]
      values <- object$subjects$reliability[object$subjects$task == task]
      cat(sprintf("  %-12s %8s %10s %8s  %-18s   %s\n",
                  task, .fmt_int(object$tasks$subjects[i]),
                  .fmt_reliability(object$tasks$reliability[i], 10L),
                  .fmt_reliability(object$tasks$sd[i], 8L),
                  .fmt_ci(object$tasks$conf.low[i], object$tasks$conf.high[i]),
                  .fmt_reliability_range(values)))
    }
    cat("\n  Means, SDs and credible intervals summarise the posterior draws.\n")
  } else if (is.null(object$boot)) {
    cat(sprintf("  %-12s %8s %20s   %s\n",
                "Task", "Subjects", "Group-level estimate",
                "By-subject estimate median [min, max]"))
    for (i in seq_len(nrow(object$tasks))) {
      task <- object$tasks$task[i]
      values <- object$subjects$reliability[object$subjects$task == task]
      cat(sprintf("  %-12s %8s %20s   %s\n",
                  task, .fmt_int(object$tasks$subjects[i]),
                  .fmt_reliability(object$tasks$reliability[i], 20L),
                  .fmt_reliability_range(values)))
    }
  } else {
    type <- switch(object$boot$type,
                   perc = "percentile", norm = "normal", basic = "basic")
    cat(sprintf("  %-12s %8s %17s %8s  %-25s   %s\n",
                "Task", "Subjects", "Group boot median", "Boot SE",
                sprintf("Boot %.0f%% CI (%s)", 100 * object$boot$level, type),
                "By-subject boot median [min, max]"))
    for (i in seq_len(nrow(object$tasks))) {
      task <- object$tasks$task[i]
      values <- object$subjects$boot_median[object$subjects$task == task]
      cat(sprintf("  %-12s %8s %17s %8s  %-25s   %s\n",
                  task, .fmt_int(object$tasks$subjects[i]),
                  .fmt_reliability(object$tasks$boot_median[i], 17L),
                  .fmt_reliability(object$tasks$boot_se[i], 8L),
                  .fmt_ci(object$tasks$conf.low[i], object$tasks$conf.high[i]),
                  .fmt_reliability_range(values)))
    }
  }

  invisible(object)
}

#' @rdname usdt_reliability
#' @export
print.usdt_reliability <- function(x, ...) summary.usdt_reliability(x, ...)

# This function prints one reliability value.
.fmt_reliability <- function(x, width = 0L) {
  value <- if (is.na(x)) "NA" else
    sub("^(-?)0[.]", "\\1.", formatC(x, format = "f", digits = 3L))
  formatC(value, width = width)
}

# This function prints the middle and range of subject values.
.fmt_reliability_range <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return("NA")
  values <- c(stats::median(x), min(x), max(x))
  sprintf("%s [%s, %s]",
          .fmt_reliability(values[1L]),
          .fmt_reliability(values[2L]),
          .fmt_reliability(values[3L]))
}
