# methods.R
# This script prints summaries for uSDT data and fitted models.
# Author: Ricardo Rey-Sáez
# Last modified: 08-09-2026

# Data summaries

#' @param x A `usdt_data` object.
#' @param ... Ignored.
#' @rdname usdt_data
#' @export
print.usdt_data <- function(x, ...) {

  m   <- x$meta
  lab <- m$labels
  fm  <- .usdt_formula(x)

  # The summary starts with the size of the data.
  cat(.rule("Data summary"), "\n\n")
  cat(sprintf("  Input:          %s (%s)\n", m$input, sub("\\(\\)$", "", m$entry)))
  cat(sprintf("  Subjects:       %d (%d in both tasks, %d in one only)\n",
              m$n_subj, m$n_both, m$n_only))
  cat(sprintf("  Trials:         %s -> %d aggregated rows (%.0f per subject)\n",
              .fmt_int(m$n_trials), m$n_rows, m$n_rows / m$n_subj))
  cat(sprintf("  Coding:         %s (condition coded %s; intercept estimates %s)\n",
              m$coding,
              if (m$coding == "deviation") "-0.5 / +0.5" else "0 / 1",
              if (m$coding == "deviation") "-c" else "z(FAR)"))
  cat(sprintf("  Parameters:     %d (%d fixed effects, %d (co)variance components)\n",
              fm$n_fixed + fm$n_var, fm$n_fixed, fm$n_var))

  # The summary shows how each column was read.
  cat("\n", .rule("Variable mapping"), "\n\n", sep = "")
  cat(sprintf("  %-14s %-10s %-15s %-13s %s\n",
              "Variable", "Task", "Column", "Signal", "Noise"))
  .print_var_rows("subject",   m, "subject",       NULL)
  .print_var_rows("condition", m, "condition_col", "condition_levels")
  .print_var_rows("response",  m, "response_col",  "response_levels")

  # The summary notes whether response times were split at the median.
  .print_split_note(m)

  # Subject summaries help reveal coding errors.
  cat("\n", .rule("Descriptives: median [min, max] across subjects"),
      "\n\n", sep = "")
  cat(sprintf("  %-10s %11s   %17s   %17s   %s\n",
              "Task", "Trials/cell", "HR", "FAR", "d' (method-of-moments)"))
  for (k in c("direct", "indirect")) {
    s <- m$descriptives[[if (k == "direct") "D" else "I"]]
    cat(sprintf("  %-10s %11.0f   %17s   %17s   %s\n", lab[[k]], s$trials,
                .fmt_rate_rng(s$hr), .fmt_rate_rng(s$far), .fmt_num_rng(s$dprime)))
  }
  edge <- sum(vapply(m$descriptives, function(s) s$edge, 0L))
  cat(if (edge == 0L) "\n  No cells at floor or ceiling.\n" else
        sprintf(paste0("\n  ! %d subject-cells at floor or ceiling. The d' above applies the\n",
                       "    Hautus correction to them; the model uses the raw counts.\n"),
                edge))

  invisible(x)
}

# Helpers for data summaries

# This function prints the role of a variable in each task.
.print_var_rows <- function(what, m, col_field, lev_field) {

  # Each task keeps its own row, so different mappings are always visible.
  for (k in c("direct", "indirect")) {
    ti  <- m$tasks[[k]]
    lev <- if (is.null(lev_field)) NULL else ti[[lev_field]]
    tag <- if (ti$dichotomized && what == "response") "  [Meyen split]" else ""
    cat(sprintf("  %-14s %-10s %-15s %-13s %s%s\n",
                if (k == "direct") what else "", m$labels[[k]],
                ti[[col_field]] %||% "-",
                if (is.null(lev)) "-" else as.character(lev[["signal"]]),
                if (is.null(lev)) "-" else as.character(lev[["noise"]]),
                tag))
  }
  invisible(NULL)
}

# This function notes an uneven median split.
.print_split_note <- function(m) {
  for (k in c("direct", "indirect")) {
    ti <- m$tasks[[k]]
    if (!ti$dichotomized || ti$dic$odd == 0L) next
    cat(sprintf(paste0("\n  ! %d of %d `%s` have an odd number of trials, so the median split\n",
                       "    is not exactly 50/50. Use ties = \"random\" to remove the bias.\n"),
                ti$dic$odd, length(ti$dic$p), ti$subject))
  }
  invisible(NULL)
}

# This function lists the estimated and fixed criteria.
.criteria_line <- function(m, fm, lab) {
  parts <- character(0)
  for (k in c("direct", "indirect")) {
    nm <- if (k == "direct") "c_D" else "c_I"
    parts <- c(parts, if (nm %in% fm$criteria)
      paste0(lab[[k]], " estimated") else
      sprintf("%s fixed to 0 (Meyen split, mean |c| = %.4f)",
              lab[[k]], m$criterion[[k]]$mean_abs))
  }
  paste(parts, collapse = ", ")
}

# This function prints the median and range of a rate.
.fmt_rate_rng <- function(v) {
  f <- function(z) sub("^0\\.", ".", formatC(z, format = "f", digits = 2L))
  sprintf("%4s [%4s, %4s]", f(stats::median(v)), f(min(v)), f(max(v)))
}

# This function prints the median and range of a signed value.
.fmt_num_rng <- function(v) {
  f <- function(z) formatC(z, format = "f", digits = 2L)
  sprintf("%5s [%5s, %5s]", f(stats::median(v)), f(min(v)), f(max(v)))
}

# This function replaces a missing value with a default.
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a

# Model summaries

#' @param object An `hsdt` object.
#' @param x An `hsdt` object.
#' @rdname hsdt
#' @export
summary.hsdt <- function(object, ...) {

  m   <- object$data$meta
  lab <- m$labels
  dg  <- object$diagnostics
  z   <- stats::qnorm(1 - (1 - object$level) / 2)

  # The summary starts with the fitted model and its diagnostics.
  cat(.rule("Model summary"), "\n\n")
  cat(sprintf("  Subjects:       %d\n", m$n_subj))
  cat(sprintf("  Observations:   %d aggregated rows (%s trials)\n",
              m$n_rows, .fmt_int(m$n_trials)))
  cat(sprintf("  Family:         binomial (probit)\n"))
  cat(sprintf("  Coding:         %s\n", m$coding))
  cat(sprintf("  Criteria:       %s\n", .criteria_line(m, object$design, lab)))
  cat(sprintf("  Estimation:     lme4::glmer (%s)%s\n", dg$optimizer,
              if (isTRUE(dg$retried))
                paste0(", tried ", paste(dg$tried, collapse = " -> ")) else ""))
  cat(sprintf("  Convergence:    %s\n", if (dg$converged) "TRUE" else "FALSE"))
  if (!is.null(object$boot)) {
    b <- object$boot
    # Boundary replicates stay in the interval, so their number is reported.
    bound <- if (b$retained[["boundary"]] > 0L)
      sprintf("; %s at the boundary", .fmt_int(b$retained[["boundary"]])) else ""
    cat(if (b$complete)
      sprintf("  Bootstrap:      %s usable replicates (%s attempts%s)\n",
              .fmt_int(b$usable), .fmt_int(b$attempted), bound) else
      sprintf("  Bootstrap:      %s of %s usable (%s attempts%s; incomplete)\n",
              .fmt_int(b$usable), .fmt_int(b$requested),
              .fmt_int(b$attempted), bound))
  }

  # The fixed-effects table includes estimates and confidence intervals.
  cat("\n", .rule("Fixed effects"), "\n\n", sep = "")
  cat(.par_header(level = object$level), sep = "")
  fx <- lme4::fixef(object$fit)
  se <- sqrt(diag(as.matrix(stats::vcov(object$fit))))
  for (nm in names(fx)) {
    p <- .term_parts(nm, lab)
    estimate <- fx[[nm]] * if (nm %in% c("c_D", "c_I")) -1 else 1
    st <- estimate / se[[nm]]
    cat(.par_row(p[1L], p[2L], estimate, se[[nm]],
                 estimate - z * se[[nm]], estimate + z * se[[nm]],
                 st, 2 * stats::pnorm(-abs(st))))
  }

  # The random-effects table names each type of estimate.
  cat("\n", .rule("Random effects"), "\n\n", sep = "")
  cat(.par_header(test = FALSE, level = object$level), sep = "")
  re <- object$pars$re
  blocks <- list(criteria = c("c_D", "c_I"), sensitivity = c("d_D", "d_I"))
  for (block in names(blocks)) {
    selected <- which(re$term %in% blocks[[block]])
    for (i in selected) {
      p  <- .random_term_parts(re$term[i], lab)
      ci <- .log_ci(re$sd[i], re$se[i], object$level)
      cat(.par_row(p[1L], p[2L], re$sd[i], re$se[i], ci[1L], ci[2L]))
    }
    if (block == "criteria" && !is.null(object$pars$criterion_cor)) {
      cr <- object$pars$criterion_cor
      ci <- .fisher_z(cr$estimate, cr$se, object$level)
      cat(.par_row("cor(c)", "both", cr$estimate, cr$se,
                   ci$conf.low, ci$conf.high))
    }
    # The table keeps the interval selected by the fit.
    if (block == "sensitivity") {
      h2 <- object$tests[object$tests$term == "correlation", ]
      cat(.par_row("cor(d')", "both", h2$estimate, h2$se,
                   h2$conf.low, h2$conf.high))
    }
  }

  # Each hypothesis has its own table.
  cat("\n", .rule("Hypotheses"), "\n\n", sep = "")
  .print_hypotheses(object, lab)

  # The notes show available follow-up analyses and missing results.
  unavailable <- unique(stats::na.omit(object$tests$reason[
    object$tests$status == "not estimable"]))
  incomplete <- !is.null(object$boot) && !isTRUE(object$boot$complete)
  if (is.null(object$boot) || incomplete || length(unavailable)) {
    cat("\n", .rule("Notes"), "\n\n", sep = "")
    if (is.null(object$boot)) {
      cat("  Use usdt_boot(fit, nsim = 1000, ncores = 4) for bootstrap CIs.\n")
    }
    if (incomplete) {
      b <- object$boot
      cat(if (isTRUE(b$summary_available))
        sprintf("  Bootstrap target not reached. Hypotheses use %s usable replicates.\n",
                .fmt_int(b$usable)) else
        "  Bootstrap target not reached. Original hypothesis results are shown.\n")
    }
    for (reason in unavailable) cat("  Not estimable: ", reason, ".\n", sep = "")
  }
  invisible(object)
}

#' @rdname hsdt
#' @export
print.hsdt <- function(x, ...) summary.hsdt(x, ...)

# Helpers for model summaries

# This function separates a parameter name from its task.
.term_parts <- function(nm, lab) {
  switch(nm,
         c_D = c("criterion", lab[["direct"]]),
         c_I = c("criterion", lab[["indirect"]]),
         d_D = c("d'",        lab[["direct"]]),
         d_I = c("d'",        lab[["indirect"]]),
         c(nm, "-"))
}

# This function names estimates from the random-effects distribution.
.random_term_parts <- function(nm, lab) {
  switch(nm,
         c_D = c("sd(criterion)", lab[["direct"]]),
         c_I = c("sd(criterion)", lab[["indirect"]]),
         d_D = c("sd(d')",        lab[["direct"]]),
         d_I = c("sd(d')",        lab[["indirect"]]),
         c(nm, "-"))
}

# This function names the columns of a parameter table.
.par_header <- function(estimate = "Estimate", test = TRUE, level = 0.95) {
  ci_label <- sprintf("%.0f%% CI", 100 * level)
  if (test) sprintf("  %-14s %-9s %9s %8s  %-18s %7s %9s\n",
                    "Parameter", "Task", estimate, "SE",
                    ci_label, "z", "p-value")
  else      sprintf("  %-14s %-9s %9s %8s  %-18s\n",
                    "Parameter", "Task", estimate, "SE", ci_label)
}

# This function formats one row of a parameter table.
.par_row <- function(parameter, task, est, se, lo, hi,
                     stat = NULL, p = NULL) {
  base <- sprintf("  %-14s %-9s %9s %8s  %-18s", parameter, task,
                  .fmt_n(est, 4L, 9L), .fmt_n(se, 4L, 8L), .fmt_ci(lo, hi))
  if (is.null(stat)) paste0(sub("\\s+$", "", base), "\n")
  else paste0(base, sprintf(" %7.2f %9s\n", stat, .fmt_p(p)))
}

# This function keeps valid intervals for random-effect deviations.
.log_ci <- function(sd, se, level) {
  if (!is.finite(sd) || !is.finite(se) || sd <= 0) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(1 - (1 - level) / 2)
  exp(log(sd) + c(-1, 1) * z * se / sd)
}

# This function prints a separate table for each hypothesis.
.print_hypotheses <- function(object, lab) {

  t <- object$tests
  delta <- .usdt_chars()$delta
  titles <- c(
    H1 = sprintf("Group-level sensitivity difference (%sd' = %s d' - %s d')",
                 delta, lab[["direct"]], lab[["indirect"]]),
    H2 = "Correlation between sensitivities across tasks",
    H3 = sprintf("Latent regression of %s d' on %s d'",
                 lab[["indirect"]], lab[["direct"]])
  )
  names <- c(correlation = "rho", intercept = "Intercept", slope = "Slope")
  names[["d'(direct) - d'(indirect)"]] <- paste0(delta, "d' (D - I)")

  for (hypothesis in c("H1", "H2", "H3")) {
    rows <- t[t$hypothesis == hypothesis, , drop = FALSE]
    cat(sprintf("%s: %s\n", hypothesis, titles[[hypothesis]]))
    cat(.hypothesis_header(object), sep = "")
    for (i in seq_len(nrow(rows))) {
      cat(.hypothesis_row(names[[rows$term[i]]], rows[i, ],
                          bootstrap = .has_boot_summary(object)))
    }
    if (hypothesis != "H3") cat("\n")
  }

  invisible(NULL)
}

# This function prints the columns used by a hypothesis table.
.hypothesis_header <- function(object) {
  if (!.has_boot_summary(object)) {
    return(sprintf("  %-14s %9s %8s  %-18s %7s %9s\n",
                   "Parameter", "Estimate", "SE",
                   sprintf("%.0f%% CI", 100 * object$level), "z", "p-value"))
  }
  type <- switch(object$boot$type,
                 perc = "percentile", norm = "normal", basic = "basic")
  sprintf("  %-14s %11s %8s  %-25s %12s\n",
          "Parameter", "Estimate", "Boot SE",
          sprintf("Boot %.0f%% CI (%s)", 100 * object$level, type),
          "Boot p-value")
}

# This function checks whether bootstrap summaries are available.
.has_boot_summary <- function(object) {
  !is.null(object$boot) && isTRUE(object$boot$summary_available)
}

# This function prints one result from a hypothesis table.
.hypothesis_row <- function(parameter, row, bootstrap) {
  if (bootstrap) {
    return(sprintf("  %-14s %11s %8s  %-25s %12s\n",
                   parameter, .fmt_n(row$estimate, 4L, 11L),
                   .fmt_n(row$se, 4L, 8L),
                   .fmt_ci(row$conf.low, row$conf.high),
                   .fmt_p(row$p.value)))
  }
  sprintf("  %-14s %9s %8s  %-18s %7s %9s\n",
          parameter, .fmt_n(row$estimate, 4L, 9L),
          .fmt_n(row$se, 4L, 8L),
          .fmt_ci(row$conf.low, row$conf.high),
          if (is.na(row$statistic)) "NA" else
            formatC(row$statistic, format = "f", digits = 2, width = 7),
          .fmt_p(row$p.value))
}
