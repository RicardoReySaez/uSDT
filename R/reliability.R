# reliability.R
# This script estimates reliability for each task.
# Author: Ricardo Rey-Sáez
# Last modified: 08-09-2026

#' Reliability of the direct and indirect measures
#'
#' Estimates how much of the spread in `d'` comes from real differences between
#' subjects, and how much comes from the noise of a limited number of trials. A
#' value close to one means that the task separates subjects well. A value
#' close to zero means that most of the observed spread is measurement error.
#'
#' @param object An `hsdt` object from [hsdt()], which may also carry the
#'   results of [usdt_boot()]. The summary method takes the
#'   `usdt_reliability` object that this function returns.
#'
#' @return An object of class `usdt_reliability`. Its `tasks` data frame gives
#'   one reliability per task, with the two variances it comes from. Its
#'   `subjects` data frame gives the `d'`, the measurement variance and the
#'   reliability of every subject. Bootstrap intervals appear in both when they
#'   are available.
#'
#' @details
#' For subject `j` in task `t` the reported value is
#' `tau2_t / (tau2_t + v_tj)`. The first term, `tau2_t`, is the variance of the
#' sensitivity between subjects, which the model estimates. The second term,
#' `v_tj`, is the variance of the sensitivity that the trials of that subject
#' can support on their own.
#'
#' `v_tj` depends on the number of trials and on the position of the subject on
#' the response curve. A cell contributes `n * dnorm(eta)^2 / (p * (1 - p))`,
#' so two subjects with the same number of trials can differ in precision. The
#' calculation also discounts the information spent on estimating the
#' criterion. When the model has fixed a criterion to zero there is nothing to
#' discount, and the same formula applies.
#'
#' This measure looks at each task alone and never uses the other task, so it
#' describes the information the data of one subject actually carry. The
#' estimates drawn by [plot.hsdt()] are more precise than this, because the
#' model there does borrow information across subjects and tasks.
#'
#' The value reported for a whole task replaces `v_tj` by its average. For a
#' subject drawn at random, the variance of a single measurement is
#' `tau2_t + E(v_tj)`, so the ratio says which share of that spread is real.
#' Averaging the reliabilities of the individual subjects would answer a
#' different question. The average over subjects treats them all as equally
#' important, which suits a sample that represents the population of interest.
#'
#' `v_tj` is the classical sampling variance of `d'` of Gourevitch and Galanter
#' (1967), which [sdt_moments()] reports as `var_gg`. It also equals the
#' standard error that a probit regression would give for the sensitivity of
#' that subject alone. The three are the same formula, evaluated at different
#' points. `var_gg` uses the observed rates of the subject, while `v_tj` uses
#' the rates that the hierarchical model predicts.
#'
#' When `object` carries a usable bootstrap, the function recomputes the
#' reliability in every retained replicate, each one with its own variance and
#' its own subject estimates. Percentile intervals remain available at a
#' boundary. Normal and basic intervals work on the logit scale, so they are
#' missing when a reliability reaches zero or one.
#'
#' @references
#' Gourevitch, V., & Galanter, E. (1967). A significance test for one parameter
#' isosensitivity functions. *Psychometrika*.
#'
#' @seealso [hsdt()], [usdt_boot()], [sdt_moments()]
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' df <- usdt_simulate(n_subj = 30, n_trials = 80)
#' data <- usdt_data_long(
#'   df,
#'   task_col = "task",
#'   task_levels = c(direct = "D", indirect = "I"),
#'   subject_col = "subj",
#'   condition_col = "cond",
#'   condition_levels = c(signal = 1, noise = 0),
#'   response_col = "response",
#'   response_levels = c(signal = 1, noise = 0)
#' )
#' usdt_reliability(hsdt(data))
#' }
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

  agg    <- object$data$agg
  labels <- object$data$meta$labels
  fixed  <- lme4::fixef(object$fit)
  random <- lme4::ranef(object$fit, condVar = FALSE)[["subj"]]

  results <- lapply(c("D", "I"), function(task) {
    slope   <- paste0("d_", task)
    columns <- c(intersect(paste0("c_", task), object$design$criteria), slope)
    tau2    <- object$pars$est[[if (task == "D") "s2_D" else "s2_I"]]
    rows    <- agg[agg$task == task, ]
    who     <- unique(as.character(rows$subj))

    # Each subject supplies the cells that inform their own sensitivity.
    values <- vapply(who, function(s) {
      cells <- rows[as.character(rows$subj) == s, , drop = FALSE]
      coef  <- fixed[columns] + unlist(random[s, columns, drop = FALSE],
                                       use.names = FALSE)
      eta   <- drop(as.matrix(cells[, columns, drop = FALSE]) %*% coef)
      c(dprime   = unname(coef[[slope]]),
        variance = .sensitivity_variance(cells, columns, slope, eta))
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

# This function estimates a subject's trial-level sensitivity variance.
.sensitivity_variance <- function(cells, columns, slope, eta) {

  # The expected information of a binomial probit weights each cell.
  p <- stats::pnorm(eta)
  w <- cells$n * stats::dnorm(eta)^2 / (p * (1 - p))

  D <- as.matrix(cells[, columns, drop = FALSE])
  tryCatch(solve(crossprod(D * sqrt(w)))[slope, slope],
           error = function(e) NA_real_)
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
      local_variance[, j] <- .sensitivity_variance_many(
        cells, columns, slope, eta)
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

# This function estimates several local variances at once.
.sensitivity_variance_many <- function(cells, columns, slope, eta) {
  p <- stats::pnorm(eta)
  weight <- sweep(stats::dnorm(eta)^2 / (p * (1 - p)),
                  2L, cells$n, "*")
  design <- as.matrix(cells[, columns, drop = FALSE])
  information <- function(i, j) {
    rowSums(sweep(weight, 2L, design[, i] * design[, j], "*"))
  }

  slope_column <- match(slope, columns)
  if (length(columns) == 1L) {
    variance <- 1 / information(slope_column, slope_column)
  } else {
    other <- setdiff(seq_along(columns), slope_column)
    iss <- information(slope_column, slope_column)
    ioo <- information(other, other)
    iso <- information(slope_column, other)
    variance <- ioo / (iss * ioo - iso^2)
  }
  variance[!is.finite(variance) | variance <= 0] <- NA_real_
  variance
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
  if (is.null(object$boot)) {
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
