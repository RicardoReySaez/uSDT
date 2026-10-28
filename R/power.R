# power.R
# Power of the tests of the three uSDT hypotheses for a planned study
# Author: Ricardo Rey-Sáez
# Last modified: 28-10-2026

# Public functions

#' Power of a planned uSDT study
#'
#' Computes the power of the tests of the three uSDT hypotheses, the
#' difference between mean sensitivities (H1), their correlation (H2) and the
#' latent regression (H3), for a study that samples a given number of subjects,
#' each with a given number of trials per task, from a population described by
#' [usdt_population()].
#'
#' @param population A `usdt_population` object from [usdt_population()]. With
#'   several scenarios, the power of each.
#' @param subjects Number of subjects. Several values give one design each.
#' @param trials Number of trials per subject in each task, signal and noise
#'   trials together: one value for both tasks, or one for each as
#'   `c(direct = , indirect = )` or `list(direct = , indirect = )`. Several
#'   values give one design each; a value shared by both tasks keeps them
#'   equal, and values given for each task are crossed.
#' @param hypothesis Character vector of hypotheses, written as in
#'   [usdt_equivalence()]: `"<quantity> <op> <value>"` with the quantity
#'   `diff`, `rho`, `slope` or `intercept` and the operator `=`, `<` or `>`, or
#'   `"<quantity> in [a, b]"` and `"<quantity> out [a, b]"`. The default tests
#'   the three point nulls that [hsdt()] reports.
#' @param alpha Significance level of every test, between 0 and 0.5 (default
#'   0.05).
#' @param signal Proportion of signal trials in each task (default 0.5), one
#'   value for both tasks or one for each. Each task must have a whole number
#'   of signal and of noise trials.
#' @param dichotomize Which tasks are split at the median, as in
#'   [usdt_data_tasks()]: `"none"`, `"indirect"`, `"direct"`, `"both"`, or one
#'   logical value for each task. The split fixes the criterion of a task at
#'   zero, so its criterion values in `population` are not used.
#' @param combine `"each"` reports the power of each hypothesis; `"all"` also
#'   reports the probability that every one of them rejects H0 in the same
#'   study.
#'
#' @return An object of class `usdt_power`: a data frame with one row per
#'   scenario, design and hypothesis, and the columns `scenario`, `subjects`,
#'   `trials_direct`, `trials_indirect`, `hypothesis`, `test`, `H1`, `H0`,
#'   `true` (the population value of the tested quantity), `se` (its expected
#'   standard error), `power`, `estimable` (the share of studies in which the
#'   analysis can estimate it), `in_null` (whether the population lies in H0,
#'   so that `power` is the type I error rate), and `reliability_direct` and
#'   `reliability_indirect`. With `combine = "all"`, a row per scenario and
#'   design whose `hypothesis` is `"all"` gives the joint power. It prints one
#'   report for a single design and a table for several.
#'
#' @details
#' # The approximation
#'
#' Each subject's estimated \eqn{d'} is treated as its true \eqn{d'} plus a
#' measurement error. The variance of that error, \eqn{v}, is the one
#' [usdt_reliability()] reports, averaged over the subjects of the population:
#' it falls with the number of trials and grows with extreme rates. The
#' observed \eqn{d'} of the subjects then have covariance
#' \eqn{S = \Sigma + \mathrm{diag}(v_D, v_I)}, where \eqn{\Sigma} is the
#' covariance of the population, and the reliability of each task is
#' \eqn{\sigma^2 / (\sigma^2 + v)}. Trials act on the tests only through
#' \eqn{v}, and subjects through their number \eqn{N}.
#'
#' With \eqn{N} subjects, the sample mean of their \eqn{d'} is normal with
#' covariance \eqn{S / N} and their sample covariance has a Wishart
#' distribution. `usdt_power()` draws both 20,000 times, subtracts the error
#' variance to estimate \eqn{\Sigma}, and in each draw runs every test as the
#' fitted model would, with the standard errors it would compute at those
#' estimates. Power is the share of draws that reject H0, with a Monte Carlo
#' error of at most 0.004. The draws are the same for every design, so power
#' changes smoothly with the design, and they do not change the random number
#' stream of the session.
#'
#' A draw whose estimated covariance of the sensitivities is not positive
#' definite (a variance of zero or a correlation of -1 or 1) is a fit on the
#' boundary. There, the fitted model reports the difference but marks the
#' correlation, slope and intercept as not estimable, so the draw counts as a
#' failure to reject them, and `estimable` reports how often that happens.
#' With a direct task of low reliability \eqn{R_D}, this probability is close
#' to \eqn{\Phi(-R_D \sqrt{N / 2})}, and more direct trials lower it faster
#' than more subjects.
#'
#' # The tests
#'
#' Each hypothesis is decided as the analysis would decide it: `"rho = 0"` and
#' `"slope = 0"` by the test on the covariance of [usdt_tests()], and every
#' other statement by the tests of [usdt_equivalence()], with the correlation
#' on the Fisher-\eqn{z} scale. `se` is the standard error at the population
#' values. When the population value lies in H0, `power` is the type I error
#' rate, and the report says so.
#'
#' # Limits
#'
#' The approximation assumes the equal-variance model that [hsdt()] fits, so
#' every `sd_signal` of `population` must be 1, and it assumes the estimated
#' standard deviations of the criteria stay away from zero. It follows the
#' fitted model closely when the reliability of both tasks is moderate or
#' high. When the direct task is unreliable, boundary fits are more frequent
#' in the fitted model than here, and the power of H2 and H3 is somewhat
#' optimistic.
#'
#' @seealso [usdt_population()], [usdt_tests()], [usdt_equivalence()],
#'   [usdt_reliability()]
#'
#' @examples
#' pop <- usdt_population(dprime    = c(direct = 0.15, indirect = 0.30),
#'                        sd_dprime = c(direct = 0.35, indirect = 0.30),
#'                        cor_dprime = 0.3, criterion = 0.1)
#'
#' # 80 subjects, 64 direct trials and 200 indirect trials split at the median
#' usdt_power(pop, subjects = 80, trials = c(direct = 64, indirect = 200),
#'            dichotomize = "indirect")
#'
#' # Several numbers of subjects, and one-sided and equivalence tests
#' usdt_power(pop, subjects = c(40, 80, 160),
#'            trials = c(direct = 64, indirect = 200),
#'            hypothesis = c("diff > 0", "intercept > 0", "rho in [-0.1, 0.1]"),
#'            dichotomize = "indirect")
#'
#' # The same design under weaker correlations
#' usdt_power(update(pop, cor_dprime = c(0.1, 0.2)), subjects = 80,
#'            trials = c(direct = 64, indirect = 200), dichotomize = "indirect")
#'
#' @export
usdt_power <- function(population, subjects, trials,
                       hypothesis  = c("diff = 0", "rho = 0", "intercept = 0"),
                       alpha       = 0.05,
                       signal      = 0.5,
                       dichotomize = c("none", "indirect", "direct", "both"),
                       combine     = c("each", "all")) {

  # The function checks the population, the design and the tests.
  if (!inherits(population, "usdt_population")) {
    .usdt_stop("`population` must come from usdt_population(), not a plain ",
               class(population)[1L], ".")
  }
  if (any(c(population$sd_signal_D, population$sd_signal_I) != 1)) {
    .usdt_stop("the power approximation assumes the equal-variance model ",
               "that hsdt() fits, so `sd_signal` must be 1.")
  }
  subjects <- .design_subjects(subjects)
  design <- .design_trials(trials, signal)
  if (!is.character(hypothesis) || !length(hypothesis) || anyNA(hypothesis)) {
    .usdt_stop("`hypothesis` must be one or more strings such as ",
               "\"intercept > 0\" or \"diff in [-0.2, 0.2]\".")
  }
  tests <- lapply(hypothesis, .parse_hypothesis)
  .check_scalar_number(alpha, "alpha", lower = 0, upper = 0.5,
                       open_lower = TRUE, open_upper = TRUE)
  split <- .check_dichotomize(dichotomize)
  combine <- match.arg(combine)

  # The error variances depend on the trials and the scenario, not on N, and
  # the chi variables of the draws on N only.
  base <- .moment_base()
  chi <- lapply(subjects, function(N) .moment_chi(base, N))
  rows <- list()
  for (s in seq_len(nrow(population))) {
    pop <- as.list(as.data.frame(population)[s, ])
    for (t in seq_len(nrow(design$trials))) {
      cells <- list(signal = unlist(design$signal[t, ]),
                    noise = unlist(design$noise[t, ]))
      v <- .error_variance(pop, cells, split)
      reliability <- unname(c(pop$sigma_D^2, pop$sigma_I^2) /
                              (c(pop$sigma_D^2, pop$sigma_I^2) + v))
      for (i in seq_along(subjects)) {
        res <- .power_approximate(pop, subjects[i], v, tests, alpha, base,
                                  chi[[i]])
        rows[[length(rows) + 1L]] <- .power_rows(
          res, tests, combine, scenario = pop$scenario, subjects = subjects[i],
          trials = unlist(design$trials[t, ]), reliability = reliability)
      }
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  structure(out, class = c("usdt_power", "data.frame"), alpha = alpha,
            method = "approximate", draws = .moment_draws, split = split,
            signal = design$proportion, combine = combine,
            labels = attr(population, "labels"),
            parameterization = attr(population, "parameterization"))
}

#' @param x A `usdt_power` object.
#' @param digits Number of decimals printed.
#' @param ... Ignored.
#' @rdname usdt_power
#' @export
print.usdt_power <- function(x, digits = 3L, ...) {
  labels <- attr(x, "labels")
  chars <- .usdt_chars()
  each <- x[x$hypothesis != "all", , drop = FALSE]
  designs <- unique(each[c("scenario", "subjects", "trials_direct",
                           "trials_indirect")])
  n_scenarios <- length(unique(each$scenario))
  right <- paste0(attr(x, "method"),
                  if (n_scenarios > 1L) paste0(" ", chars$dot, " ",
                                               n_scenarios, " scenarios"),
                  if (nrow(designs) > n_scenarios)
                    paste0(" ", chars$dot, " ", nrow(designs), " designs"))
  cat(.rule("Power", right), "\n\n")

  if (nrow(designs) == 1L) {
    .print_power_card(x, labels, digits)
  } else {
    .print_power_table(x, designs, labels, digits)
  }

  # The notes explain the marks and the approximation.
  cat("\n")
  if (any(each$in_null)) {
    cat("  ^ The population lies in H0: the rejection rate is the type I error",
        "rate.\n")
  }
  if (nrow(designs) > 1L && any(each$estimable < 0.99)) {
    cat("  * Fewer than 99% of the studies estimate it (a fit on the",
        "boundary);\n    see the `estimable` column.\n")
  }
  cat("  Approximate power: the share of ", .fmt_int(attr(x, "draws")),
      " draws of the subjects' sample\n  mean and covariance in which the ",
      "test rejects H0, at alpha = ", format(attr(x, "alpha")), ".\n",
      sep = "")
  invisible(x)
}

# A subset of the rows is a plain data frame, because the report needs every
# row of a design.
#' @export
`[.usdt_power` <- function(x, ...) {
  out <- NextMethod()
  if (is.data.frame(out)) {
    class(out) <- setdiff(class(out), "usdt_power")
    for (a in c("alpha", "method", "draws", "split", "signal", "combine",
                "labels", "parameterization")) {
      attr(out, a) <- NULL
    }
  }
  out
}

# Internal functions

# This function checks the numbers of subjects.
.design_subjects <- function(subjects) {
  if (!is.numeric(subjects) || !length(subjects) || anyNA(subjects) ||
      any(!is.finite(subjects)) || any(abs(subjects - round(subjects)) > 1e-8) ||
      any(subjects < 3)) {
    .usdt_stop("`subjects` must be one or more whole numbers of at least 3.")
  }
  unique(as.integer(round(subjects)))
}

# This function checks the trials of each task and the proportion of signal
# trials, and returns one row per design with the trials of each task and their
# signal and noise counts. A shared value keeps both tasks equal; values given
# for each task are crossed.
.design_trials <- function(trials, signal) {
  whole <- function(x, what) {
    if (is.list(x)) x <- unlist(x, use.names = FALSE)
    if (!is.numeric(x) || !length(x) || anyNA(x) || any(!is.finite(x)) ||
        any(x <= 0) || any(abs(x - round(x)) > 1e-8)) {
      .usdt_stop("`", what, "` must be one or more positive whole numbers.")
    }
    unique(as.integer(round(x)))
  }
  grid <- if (is.null(names(trials))) {
    v <- whole(trials, "trials")
    data.frame(direct = v, indirect = v)
  } else {
    t <- .per_task(trials, "trials")
    expand.grid(direct = whole(t$direct, "trials$direct"),
                indirect = whole(t$indirect, "trials$indirect"),
                KEEP.OUT.ATTRS = FALSE)
  }

  proportion <- .per_task(signal, "signal")
  for (k in c("direct", "indirect")) {
    p <- proportion[[k]]
    if (!is.numeric(p) || length(p) != 1L || is.na(p) || p <= 0 || p >= 1) {
      .usdt_stop("`signal` must be one proportion inside (0, 1) for each ",
                 "task.")
    }
  }
  proportion <- unlist(proportion)
  counts <- sweep(as.matrix(grid), 2L, proportion, "*")
  if (any(abs(counts - round(counts)) > 1e-8)) {
    at <- which(abs(counts - round(counts)) > 1e-8, arr.ind = TRUE)[1L, ]
    task <- colnames(grid)[at[["col"]]]
    .usdt_stop(grid[at[["row"]], task], " ", task, " trials with `signal = ",
               proportion[[task]], "` give ", counts[at[["row"]], task],
               " signal trials. Choose trials and `signal` that give whole ",
               "numbers of signal and noise trials.")
  }
  signal_trials <- as.data.frame(round(counts))
  noise_trials <- grid - signal_trials
  if (any(signal_trials < 1) || any(noise_trials < 1)) {
    .usdt_stop("every task needs at least one signal and one noise trial.")
  }
  list(trials = grid, signal = signal_trials, noise = noise_trials,
       proportion = proportion)
}

# This function writes the rows of one scenario and design.
.power_rows <- function(res, tests, combine, scenario, subjects, trials,
                        reliability) {
  described <- lapply(tests, function(h) {
    name <- .bf_name(h$quantity)
    statements <- .hypothesis_statements(name, h)
    data.frame(hypothesis = .bf_label(name, h),
               quantity = h$quantity,
               test = switch(h$op, "=" = "Two-sided test",
                             "<" = , ">" = "One-sided test",
                             "in" = "Equivalence test (TOST)",
                             "Minimal-effect test"),
               H1 = statements[["H1"]], H0 = statements[["H0"]],
               stringsAsFactors = FALSE)
  })
  out <- cbind(do.call(rbind, described),
               res$rows[c("true", "se", "power", "estimable", "in_null")])
  if (combine == "all") {
    out <- rbind(out, data.frame(
      hypothesis = "all", quantity = NA_character_, test = "Joint test",
      H1 = "every H1", H0 = "any H0", true = NA_real_, se = NA_real_,
      power = res$joint, estimable = res$joint_estimable, in_null = FALSE,
      stringsAsFactors = FALSE))
  }
  cbind(data.frame(scenario = scenario, subjects = subjects,
                   trials_direct = trials[["direct"]],
                   trials_indirect = trials[["indirect"]]),
        out,
        reliability_direct = reliability[1L],
        reliability_indirect = reliability[2L])
}

# This function prints the report of a single design.
.print_power_card <- function(x, labels, digits) {
  first <- x[1L, ]
  number <- function(v) {
    if (is.na(v)) return("")
    formatC(v, format = "f", digits = digits)
  }
  split <- attr(x, "split")
  task <- function(k) {
    paste0(first[[paste0("trials_", k)]], " ", labels[[k]],
           if (split[[k]]) " (split at the median)")
  }
  cat("  Design       ", first$subjects, " subjects; trials per subject: ",
      task("direct"), ", ", task("indirect"), "\n", sep = "")
  cat("  Reliability  ", labels[["direct"]], " ", number(first$reliability_direct),
      ", ", labels[["indirect"]], " ", number(first$reliability_indirect), "\n\n",
      sep = "")

  # A mark after the power keeps its own column, so the numbers stay aligned.
  marked <- any(x$in_null)
  mark <- function(flag) if (!marked) "" else if (flag) "^" else " "
  number_h <- c(diff = "H1", rho = "H2", intercept = "H3", slope = "H3")
  rows <- lapply(seq_len(nrow(x)), function(i) {
    r <- x[i, ]
    if (r$hypothesis == "all") {
      return(c("All at once", "Joint test", "", "", number(r$estimable),
               paste0(number(r$power), mark(FALSE))))
    }
    c(paste(number_h[[r$quantity]], r$hypothesis), r$test, number(r$true),
      number(r$se), number(r$estimable),
      paste0(number(r$power), mark(r$in_null)))
  })
  rows <- rbind(c("Hypothesis", "Test", "True", "SE", "Estimable",
                  paste0("Power", mark(FALSE))),
                do.call(rbind, rows))
  widths <- c(-max(nchar(rows[, 1L])), -max(nchar(rows[, 2L])), 6L, 6L, 9L, 6L)
  for (i in seq_len(nrow(rows))) {
    cat("  ", paste(.pad(rows[i, ], widths), collapse = "  "), "\n", sep = "")
  }
  invisible(NULL)
}

# This function pads each value to its own width, by characters; a negative
# width pads on the right.
.pad <- function(x, widths) {
  mapply(function(value, width) formatC(value, width = width), x, widths,
         USE.NAMES = FALSE)
}

# This function prints one row per design and one power column per
# hypothesis.
.print_power_table <- function(x, designs, labels, digits) {
  hypotheses <- unique(x$hypothesis)
  number_h <- c(diff = "H1", rho = "H2", intercept = "H3", slope = "H3")
  legend <- vapply(hypotheses, function(h) {
    if (h == "all") return("all: every hypothesis at once")
    i <- match(h, x$hypothesis)
    paste0(number_h[[x$quantity[i]]], ": ", h, ", ",
           sub("^(.)", "\\L\\1", x$test[i], perl = TRUE))
  }, "")
  cat(paste0("  ", legend, "\n"), sep = "")
  cat("\n")

  # The marks after each power share a column of their own width.
  marks <- ifelse(x$in_null, "^", "")
  marks <- paste0(marks, ifelse(x$hypothesis != "all" & x$estimable < 0.99,
                                "*", ""))
  width <- max(nchar(marks))

  scenarios <- length(unique(designs$scenario)) > 1L
  head <- c(if (scenarios) "Scenario", "Subjects", labels[["direct"]],
            labels[["indirect"]],
            paste0(hypotheses, strrep(" ", width)))
  cells <- lapply(seq_len(nrow(designs)), function(i) {
    d <- designs[i, ]
    at <- which(x$scenario == d$scenario & x$subjects == d$subjects &
                  x$trials_direct == d$trials_direct &
                  x$trials_indirect == d$trials_indirect)
    power <- vapply(hypotheses, function(h) {
      j <- at[x$hypothesis[at] == h]
      paste0(formatC(x$power[j], format = "f", digits = digits),
             formatC(marks[j], width = -width))
    }, "")
    c(if (scenarios) as.character(d$scenario), as.character(d$subjects),
      as.character(d$trials_direct), as.character(d$trials_indirect), power)
  })
  table <- rbind(head, do.call(rbind, cells))
  widths <- apply(table, 2L, function(column) max(nchar(column)))
  for (i in seq_len(nrow(table))) {
    cat("  ", paste(.pad(table[i, ], widths), collapse = "  "), "\n", sep = "")
  }
  invisible(NULL)
}
