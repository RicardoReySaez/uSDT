# plan.R
# Choose the number of subjects and trials of a uSDT study, and find the
# smallest effect a design detects
# Author: Ricardo Rey-Sáez
# Last modified: 28-10-2026

# Public functions

#' Plan the subjects and trials of a uSDT study
#'
#' Finds the design, the number of subjects and of trials per task, that
#' reaches a given power for the tests of the uSDT hypotheses at the lowest
#' cost, the design with the highest power within a budget, or the design
#' that estimates a quantity with a given precision. What is free is given as
#' a range and what is fixed as a single number, so the same function finds
#' the subjects for given trials, the trials for given subjects, or both.
#'
#' @inheritParams usdt_power
#' @param population A `usdt_population` object from [usdt_population()]. With
#'   several scenarios, one plan for each.
#' @param power Power to reach (default 0.80).
#' @param subjects Number of subjects: one number fixes it, and a range
#'   `c(min, max)` lets the plan choose it (default `c(10, 1000)`).
#' @param trials Number of trials per subject in each task, signal and noise
#'   trials together: one number fixes it, and a range `c(min, max)` lets the
#'   plan choose it. One value or range for both tasks, which keeps them equal,
#'   or one for each as `c(direct = , indirect = )` or
#'   `list(direct = , indirect = )`.
#' @param cost Cost of a design, which the plan minimises: a named vector with
#'   the cost of each `subject` (a session, recruitment) and of each `direct`
#'   and `indirect` trial, in any unit (missing names cost 0), so that a design
#'   costs `subjects * (subject + direct * trials_direct + indirect *
#'   trials_indirect)`; or a function of `subjects`, `direct` and `indirect`
#'   that returns the cost of each design. Needed when more than one of the
#'   subjects and the trials is free; with only one free, the plan minimises
#'   it.
#' @param budget Largest cost allowed. With a budget, the plan finds the
#'   design with the highest power that costs no more, and `power` is not
#'   used.
#' @param precision Named vector of the widths that the confidence intervals
#'   of usdt_tests() (level `1 - alpha`) should not exceed, by quantity: `diff`,
#'   `rho`, `intercept` or `slope`, e.g. `c(intercept = 0.2)`. With it, the
#'   plan finds the cheapest design whose expected interval is that narrow,
#'   and `hypothesis` and `power` are not used.
#' @param combine `"each"` plans a design for each hypothesis; `"all"` plans
#'   one design in which every hypothesis rejects H0 in the same study with
#'   the given power.
#'
#' @return An object of class `usdt_plan`: a list with `$design`, a data frame
#'   with one row per scenario and hypothesis (or quantity, for a precision
#'   target) and the columns `scenario`, `hypothesis`, `subjects`,
#'   `trials_direct`, `trials_indirect`, `power` (approximate, as in
#'   [usdt_power()]), `estimable`, `width` (for a precision target), `cost`
#'   and `limit` (the free argument whose limit the design reaches, if any);
#'   `$frontier`, the smallest number of subjects that reaches the target
#'   with each number of trials, from the closed form; and `$settings`. It
#'   prints one report per design.
#'
#' @details
#' # How the design is found
#'
#' The search follows the two approximations of [usdt_power()]. With each
#' number of trials in the ranges, on a grid of at most 400 values per task,
#' the closed form gives the smallest number of subjects that reaches the
#' power, and with it the cost of each design. The cheapest designs are then
#' checked with the draws of the sample mean and covariance, which see the
#' fits on the boundary, and the number of subjects of each is raised or
#' lowered until those draws reach the power; the search stops when the next
#' design costs more in the closed form than the best design checked so far.
#' A design whose subjects or trials reach the limit of their range is marked:
#' the range may be too narrow.
#'
#' The design chosen this way rests on the approximation. Simulate it with
#' `usdt_power(plan, method = "simulation")` before relying on it, as the
#' approximation is somewhat optimistic for H2 and H3 when the direct task is
#' unreliable.
#'
#' # Costs
#'
#' Many designs reach the same power, trading subjects for trials, so the
#' choice among them depends on what each costs. With linear costs and the
#' difference between means (H1), the cheapest number of trials of each task
#' is close to \eqn{\sqrt{(c_s / c_t) (a_t / \tau^2)}}, where \eqn{c_s} is the
#' cost of a subject, \eqn{c_t} that of a trial of the task, \eqn{a_t} the
#' error variance of a single trial (about \eqn{2\pi} near \eqn{d' = 0}) and
#' \eqn{\tau^2} the variance of the difference between the two \eqn{d'}
#' across subjects, whatever the power.
#'
#' @seealso [usdt_power()], [usdt_mdes()], [usdt_population()]
#'
#' @examples
#' pop <- usdt_population(dprime    = c(direct = 0.15, indirect = 0.30),
#'                        sd_dprime = c(direct = 0.35, indirect = 0.30),
#'                        cor_dprime = 0.3, criterion = 0.1)
#'
#' # Subjects for given trials
#' usdt_plan(pop, hypothesis = "diff = 0",
#'           trials = c(direct = 64, indirect = 200), dichotomize = "indirect")
#'
#' # Subjects and trials together, with the session time of each in seconds
#' plan <- usdt_plan(pop, hypothesis = "intercept > 0", power = 0.9,
#'                   subjects = c(20, 300),
#'                   trials = list(direct = c(32, 160), indirect = c(100, 300)),
#'                   cost = c(subject = 600, direct = 4, indirect = 1.5),
#'                   dichotomize = "indirect")
#' plan
#'
#' # The most powerful design within 30 hours of sessions
#' usdt_plan(pop, hypothesis = "rho = 0", subjects = c(20, 400),
#'           trials = list(direct = c(32, 256), indirect = c(100, 800)),
#'           cost = c(subject = 600, direct = 4, indirect = 1.5),
#'           budget = 30 * 3600, dichotomize = "indirect")
#'
#' # The intercept within an interval 0.2 wide
#' usdt_plan(pop, precision = c(intercept = 0.2),
#'           trials = c(direct = 64, indirect = 200), dichotomize = "indirect")
#'
#' \donttest{
#' # Check the chosen design by simulation
#' usdt_power(plan, method = "simulation", nsim = 50, seed = 1)
#' }
#'
#' @export
usdt_plan <- function(population,
                      hypothesis  = c("diff = 0", "rho = 0", "intercept = 0"),
                      power       = 0.80,
                      subjects    = c(10, 1000),
                      trials,
                      cost        = NULL,
                      budget      = NULL,
                      precision   = NULL,
                      combine     = c("each", "all"),
                      alpha       = 0.05,
                      signal      = 0.5,
                      dichotomize = c("none", "indirect", "direct", "both")) {

  # The function checks the population, the targets and the design space.
  .check_plan_population(population)
  if (missing(trials)) {
    .usdt_stop("`trials` is required: a number of trials per task, or a ",
               "range c(min, max) for the plan to choose.")
  }
  combine <- match.arg(combine)
  .check_scalar_number(alpha, "alpha", lower = 0, upper = 0.5,
                       open_lower = TRUE, open_upper = TRUE)
  split <- .check_dichotomize(dichotomize)
  space <- .plan_space(subjects, trials, signal)
  mode <- if (!is.null(precision)) "precision" else if (!is.null(budget))
    "budget" else "power"
  if (mode == "power") {
    .check_scalar_number(power, "power", lower = 0, upper = 1,
                         open_lower = TRUE, open_upper = TRUE)
  }
  if (mode == "precision") {
    precision <- .check_precision(precision)
    tests <- NULL
  } else {
    tests <- .check_hypotheses(hypothesis)
  }
  free <- space$free
  if (is.null(cost)) {
    if (mode == "budget") .usdt_stop("a `budget` needs the `cost` of a design.")
    if (sum(free) > 1L) {
      .usdt_stop("more than one of the subjects and the trials is free, and ",
                 "many of their combinations reach the target. Give the ",
                 "`cost` of a design to choose among them, e.g. cost = ",
                 "c(subject = 600, direct = 4, indirect = 1.5) in seconds of ",
                 "session, or cost = c(subject = 1) to use as few subjects as ",
                 "the trials allow.")
    }
    if (!any(free)) {
      .usdt_stop("nothing is free: give a range for the subjects or the ",
                 "trials, or use usdt_power() for a fixed design.")
    }
  }
  cost_fun <- .plan_cost(cost, free)
  if (mode == "budget") {
    .check_scalar_number(budget, "budget", lower = 0, open_lower = TRUE)
  }

  # Every scenario has its own plan, searched on the same grid.
  base <- .moment_base()
  designs <- list()
  frontiers <- list()
  for (s in seq_len(nrow(population))) {
    pop <- as.list(as.data.frame(population)[s, ])
    found <- .plan_scenario(pop, mode, tests, precision, power, budget,
                            combine, alpha, split, space, cost_fun, base)
    designs[[s]] <- cbind(scenario = pop$scenario, found$design)

    # A budget has no frontier: it fixes the cost, not the power.
    if (nrow(found$frontier)) {
      frontiers[[s]] <- cbind(scenario = pop$scenario, found$frontier)
    }
  }
  design <- do.call(rbind, designs)
  frontier <- if (length(frontiers)) do.call(rbind, frontiers) else
    data.frame()
  rownames(design) <- NULL
  if (nrow(frontier)) rownames(frontier) <- NULL

  structure(list(design = design, frontier = frontier,
                 settings = list(population = population, tests = tests,
                                 hypothesis = if (!is.null(tests))
                                   vapply(tests, `[[`, "", "text"),
                                 mode = mode, power = power, budget = budget,
                                 precision = precision, combine = combine,
                                 alpha = alpha, signal = signal, split = split,
                                 space = space, cost = cost)),
            class = "usdt_plan")
}

#' @param x A `usdt_plan` object.
#' @param digits Number of decimals printed.
#' @param ... Ignored.
#' @rdname usdt_plan
#' @export
print.usdt_plan <- function(x, digits = 3L, ...) {
  st <- x$settings
  labels <- attr(st$population, "labels")
  chars <- .usdt_chars()
  target <- switch(st$mode,
    power = paste0("power ", format(st$power)),
    budget = paste0("budget ", format(st$budget)),
    precision = "precision")
  scenarios <- length(unique(x$design$scenario))
  cat(.rule("Plan", paste0("approximate ", chars$dot, " ", target,
                           if (scenarios > 1L) paste0(" ", chars$dot, " ",
                                                      scenarios, " scenarios"))),
      "\n\n")

  # The design space says what the plan could choose.
  sp <- st$space
  range_text <- function(r) {
    if (length(r) == 1L) paste0(r, " (fixed)") else
      paste0(r[1L], " to ", r[2L])
  }
  cat("  Subjects     ", range_text(sp$subjects), "\n", sep = "")
  if (sp$linked) {
    cat("  Trials       ", range_text(sp$trials$direct), " in each task\n",
        sep = "")
  } else {
    cat("  Trials       ", labels[["direct"]], " ",
        range_text(sp$trials$direct), "; ", labels[["indirect"]], " ",
        range_text(sp$trials$indirect), "\n", sep = "")
  }
  if (!is.null(st$cost)) {
    cat("  Cost         ", .cost_text(st$cost, labels), "\n", sep = "")
  }
  cat("\n")

  number <- function(v) ifelse(is.na(v), "-", formatC(v, format = "f",
                                                         digits = digits))
  d <- x$design
  head <- c(if (scenarios > 1L) "Scenario", "Target", "Subjects",
            labels[["direct"]], labels[["indirect"]],
            if (st$mode == "precision") "Width" else "Power", "Estimable",
            "Cost")
  rows <- lapply(seq_len(nrow(d)), function(i) {
    r <- d[i, ]
    c(if (scenarios > 1L) as.character(r$scenario), r$hypothesis,
      if (is.na(r$subjects)) "-" else as.character(r$subjects),
      if (is.na(r$trials_direct)) "-" else as.character(r$trials_direct),
      if (is.na(r$trials_indirect)) "-" else as.character(r$trials_indirect),
      number(if (st$mode == "precision") r$width else r$power),
      number(r$estimable),
      if (is.na(r$cost)) "-" else formatC(r$cost, format = "fg", digits = 6L,
                                          big.mark = ","))
  })
  table <- rbind(head, do.call(rbind, rows))
  widths <- apply(table, 2L, function(v) max(nchar(v)))
  widths[seq_len(1L + (scenarios > 1L))] <-
    -widths[seq_len(1L + (scenarios > 1L))]
  for (i in seq_len(nrow(table))) {
    cat("  ", paste(.pad(table[i, ], widths), collapse = "  "), "\n", sep = "")
  }

  # The notes flag designs the ranges may have cut short.
  cat("\n")
  missing <- d$hypothesis[is.na(d$subjects)]
  if (length(missing)) {
    cat("  No design in the ranges reaches the target for ",
        paste(unique(missing), collapse = ", "), ".\n", sep = "")
  }
  limits <- d[!is.na(d$limit), , drop = FALSE]
  for (i in seq_len(nrow(limits))) {
    cat("  ! ", limits$hypothesis[i], ": the design is at the limit of ",
        limits$limit[i], "; a wider range may give a cheaper one.\n", sep = "")
  }
  if (st$mode != "precision") {
    cat("  Approximate power. Check the design with usdt_power(plan, method",
        "= \"simulation\").\n")
  } else {
    cat("  Expected width of the ", format(100 * (1 - st$alpha)), "% interval ",
        "of usdt_tests() at the population values.\n", sep = "")
  }
  invisible(x)
}

#' Smallest effect a uSDT study detects
#'
#' Finds, for each hypothesis, the smallest effect that a design detects with
#' the given power: the population value of the tested quantity, moved away
#' from H0 with the other values of the population fixed, at which the
#' approximate power of [usdt_power()] reaches it. For an equivalence test
#' (`in`), it finds the narrowest equivalence region, around the same centre,
#' that the design can establish.
#'
#' @inheritParams usdt_power
#' @param population A `usdt_population` object, or a `usdt_plan` from
#'   [usdt_plan()], whose designs, hypotheses and settings are then used.
#' @param power Power to reach (default 0.80).
#' @param vary Optional name of the population value to move instead of the
#'   tested quantity: `"mu_D"`, `"mu_I"`, `"sigma_D"`, `"sigma_I"`, `"rho"`,
#'   `"intercept"`, `"criterion_D"`, `"criterion_I"`, `"sd_criterion_D"` or
#'   `"sd_criterion_I"`. Needs `range`.
#' @param range With `vary`, the interval of values to search, `c(from, to)`;
#'   the result is the value nearest to `from` that reaches the power.
#'
#' @return An object of class `usdt_mdes`: a data frame with one row per
#'   scenario, design and hypothesis, and the columns `scenario`, `subjects`,
#'   `trials_direct`, `trials_indirect`, `hypothesis`, `varies` (the value
#'   moved), `population` (its value in the population), `power_population`
#'   (the power there) and `mdes` (the value at which the power is reached, or
#'   `NA` if it is not reached in the range searched). For an equivalence
#'   test, `varies` is `"margin"` and `mdes` the half-width of the narrowest
#'   region. The attribute `curve` holds the power over the values searched.
#'
#' @details
#' Moving a quantity keeps the other values of the population fixed as its
#' parameterization says (see [usdt_population()]): the difference moves the
#' indirect mean with the direct mean fixed; the intercept moves the indirect
#' mean with the means and SDs and the correlation fixed; the correlation
#' (and the slope, through it) keeps both means under the correlation
#' parameterization and the intercept under the latent regression. The search
#' runs from the value of H0 (or the nearest bound) outwards, in the
#' direction of H1, and the power along the way comes from the same draws of
#' the sample mean and covariance as in [usdt_power()], so it can level off
#' below the target when fits on the boundary are frequent; the result is
#' then `NA`.
#'
#' @seealso [usdt_power()], [usdt_plan()]
#'
#' @examples
#' pop <- usdt_population(dprime    = c(direct = 0.15, indirect = 0.30),
#'                        sd_dprime = c(direct = 0.35, indirect = 0.30),
#'                        cor_dprime = 0.3, criterion = 0.1)
#'
#' usdt_mdes(pop, subjects = 80, trials = c(direct = 64, indirect = 200),
#'           dichotomize = "indirect")
#'
#' # The narrowest equivalence margin for the difference
#' usdt_mdes(update(pop, dprime = c(direct = 0.15, indirect = 0.15)),
#'           subjects = 80, trials = c(direct = 64, indirect = 200),
#'           hypothesis = "diff in [-0.2, 0.2]", dichotomize = "indirect")
#'
#' # The SD of the direct d' the intercept test needs
#' usdt_mdes(pop, subjects = 80, trials = c(direct = 64, indirect = 200),
#'           hypothesis = "intercept > 0", vary = "sigma_D",
#'           range = c(0.05, 0.6), dichotomize = "indirect")
#'
#' @export
usdt_mdes <- function(population, subjects, trials,
                      hypothesis  = c("diff = 0", "rho = 0", "intercept = 0"),
                      power       = 0.80,
                      alpha       = 0.05,
                      signal      = 0.5,
                      dichotomize = c("none", "indirect", "direct", "both"),
                      vary        = NULL,
                      range       = NULL) {

  # A plan supplies its designs and settings, and its power unless another
  # is given.
  if (inherits(population, "usdt_plan")) {
    if (!missing(subjects) || !missing(trials)) {
      .usdt_stop("a plan gives the subjects and trials; leave them out.")
    }
    if (missing(power) && population$settings$mode == "power") {
      power <- population$settings$power
    }
    .check_vary(vary, range)
    return(.mdes_plan(population, power, vary, range))
  }
  .check_plan_population(population)
  subjects <- .design_subjects(subjects)
  design <- .design_trials(trials, signal)
  tests <- .check_hypotheses(hypothesis)
  .check_scalar_number(power, "power", lower = 0, upper = 1,
                       open_lower = TRUE, open_upper = TRUE)
  .check_scalar_number(alpha, "alpha", lower = 0, upper = 0.5,
                       open_lower = TRUE, open_upper = TRUE)
  split <- .check_dichotomize(dichotomize)
  .check_vary(vary, range)

  base <- .moment_base()
  rows <- list()
  curves <- list()
  for (s in seq_len(nrow(population))) {
    pop <- as.list(as.data.frame(population)[s, ])
    for (t in seq_len(nrow(design$trials))) {
      cells <- list(signal = unlist(design$signal[t, ]),
                    noise = unlist(design$noise[t, ]))
      for (N in subjects) {
        for (h in tests) {
          res <- .mdes_one(pop, h, N, cells, split, alpha, power, base,
                           attr(population, "parameterization"), vary, range)
          trials_t <- unlist(design$trials[t, ])
          rows[[length(rows) + 1L]] <- data.frame(
            scenario = pop$scenario, subjects = N,
            trials_direct = trials_t[["direct"]],
            trials_indirect = trials_t[["indirect"]],
            hypothesis = .bf_label(.bf_name(h$quantity), h),
            quantity = h$quantity, varies = res$varies,
            population = res$current, power_population = res$power_current,
            mdes = res$mdes, stringsAsFactors = FALSE)
          curves[[length(curves) + 1L]] <- res$curve
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  structure(out, class = c("usdt_mdes", "data.frame"), power = power,
            alpha = alpha, curve = curves, labels = attr(population, "labels"),
            split = split)
}

#' @param x A `usdt_mdes` object.
#' @rdname usdt_mdes
#' @export
print.usdt_mdes <- function(x, digits = 3L, ...) {
  labels <- attr(x, "labels")
  chars <- .usdt_chars()
  cat(.rule("Smallest detectable effect", paste0(
    "approximate ", chars$dot, " power ", format(attr(x, "power")))), "\n\n")
  number <- function(v) ifelse(is.na(v), "NA", formatC(v, format = "f",
                                                         digits = digits))
  number_h <- c(diff = "H1", rho = "H2", intercept = "H3", slope = "H3")
  several <- nrow(unique(x[c("scenario", "subjects", "trials_direct",
                             "trials_indirect")])) > 1L
  head <- c(if (several) c("Scenario", "Subjects", labels[["direct"]],
                           labels[["indirect"]]),
            "Hypothesis", "Varies", "Population", "Power there", "Smallest")
  rows <- lapply(seq_len(nrow(x)), function(i) {
    r <- x[i, ]
    c(if (several) c(as.character(r$scenario), as.character(r$subjects),
                     as.character(r$trials_direct),
                     as.character(r$trials_indirect)),
      paste(number_h[[r$quantity]], r$hypothesis), r$varies,
      number(r$population), number(r$power_population), number(r$mdes))
  })
  if (!several) {
    r <- x[1L, ]
    cat("  Design       ", r$subjects, " subjects; trials per subject: ",
        r$trials_direct, " ", labels[["direct"]], ", ", r$trials_indirect, " ",
        labels[["indirect"]], "\n\n", sep = "")
  }
  table <- rbind(head, do.call(rbind, rows))
  widths <- apply(table, 2L, function(v) max(nchar(v)))
  left <- match(c("Hypothesis", "Varies"), head)
  widths[left] <- -widths[left]
  for (i in seq_len(nrow(table))) {
    cat("  ", paste(.pad(table[i, ], widths), collapse = "  "), "\n", sep = "")
  }
  cat("\n  Approximate power, as in usdt_power(). NA: the power is not reached",
      "\n  in the range searched.\n")
  if (any(x$varies == "margin")) {
    cat("  For an equivalence test, the smallest value is the half-width of",
        "the\n  narrowest region the design can establish.\n")
  }
  invisible(x)
}

# Internal functions

# This function checks a population for the approximate engines.
.check_plan_population <- function(population) {
  if (!inherits(population, "usdt_population")) {
    .usdt_stop("`population` must come from usdt_population(), not a plain ",
               class(population)[1L], ".")
  }
  if (any(c(population$sd_signal_D, population$sd_signal_I) != 1)) {
    .usdt_stop("planning uses the power approximation, which assumes the ",
               "equal-variance model that hsdt() fits, so `sd_signal` must ",
               "be 1.")
  }
  invisible(TRUE)
}

# This function reads the hypotheses.
.check_hypotheses <- function(hypothesis) {
  if (!is.character(hypothesis) || !length(hypothesis) || anyNA(hypothesis)) {
    .usdt_stop("`hypothesis` must be one or more strings such as ",
               "\"intercept > 0\" or \"diff in [-0.2, 0.2]\".")
  }
  lapply(hypothesis, .parse_hypothesis)
}

# This function checks a precision target.
.check_precision <- function(precision) {
  quantities <- c("diff", "rho", "intercept", "slope")
  if (!is.numeric(precision) || !length(precision) || is.null(names(precision)) ||
      anyNA(precision) || any(precision <= 0) ||
      !all(names(precision) %in% quantities) || anyDuplicated(names(precision))) {
    .usdt_stop("`precision` must give a positive interval width for one or ",
               "more of ", paste0("`", quantities, "`", collapse = ", "),
               ", e.g. c(intercept = 0.2).")
  }
  precision
}

# This function checks the value usdt_mdes() moves instead of the quantity.
.mdes_vary <- c("mu_D", "mu_I", "sigma_D", "sigma_I", "rho", "intercept",
                "criterion_D", "criterion_I", "sd_criterion_D",
                "sd_criterion_I")
.check_vary <- function(vary, range) {
  if (is.null(vary)) {
    if (!is.null(range)) .usdt_stop("`range` goes with `vary`.")
    return(invisible(TRUE))
  }
  if (!is.character(vary) || length(vary) != 1L || !vary %in% .mdes_vary) {
    .usdt_stop("`vary` must be one of ", paste0("`", .mdes_vary, "`",
                                                collapse = ", "), ".")
  }
  if (!is.numeric(range) || length(range) != 2L || anyNA(range) ||
      range[1L] == range[2L]) {
    .usdt_stop("`vary` needs `range = c(from, to)`, two different values.")
  }
  invisible(TRUE)
}

# This function reads the design space: the subjects and the trials of each
# task, fixed (one number) or free (a range), and the grid of trial counts
# with whole numbers of signal and noise trials. A range shared by both tasks
# keeps them equal.
.plan_space <- function(subjects, trials, signal) {
  range_of <- function(x, what, lower) {
    if (is.list(x)) x <- unlist(x, use.names = FALSE)
    if (!is.numeric(x) || !length(x) %in% 1:2 || anyNA(x) ||
        any(!is.finite(x)) || any(abs(x - round(x)) > 1e-8) || any(x < lower)) {
      .usdt_stop("`", what, "` must be one whole number of at least ", lower,
                 ", or a range c(min, max) of them.")
    }
    x <- sort(as.integer(round(x)))
    if (length(x) == 2L && x[1L] == x[2L]) x <- x[1L]
    x
  }
  subjects <- range_of(subjects, "subjects", 3L)
  linked <- is.null(names(trials))
  t <- if (linked) list(direct = trials, indirect = trials) else
    .per_task(trials, "trials")
  ranges <- list(direct = range_of(t$direct, if (linked) "trials" else
                   "trials$direct", 2L),
                 indirect = range_of(t$indirect, if (linked) "trials" else
                   "trials$indirect", 2L))

  proportion <- .per_task(signal, "signal")
  for (k in c("direct", "indirect")) {
    p <- proportion[[k]]
    if (!is.numeric(p) || length(p) != 1L || is.na(p) || p <= 0 || p >= 1) {
      .usdt_stop("`signal` must be one proportion inside (0, 1) for each ",
                 "task.")
    }
  }
  proportion <- unlist(proportion)
  step <- vapply(proportion, .trial_step, 0L)

  axis <- function(r, k_step, what) {
    if (length(r) == 1L) {
      if (r %% k_step != 0L) {
        .usdt_stop(r, " ", what, " trials do not split into whole numbers ",
                   "of signal and noise trials with this `signal`.")
      }
      return(r)
    }
    values <- seq.int(k_step * ceiling(r[1L] / k_step),
                      k_step * floor(r[2L] / k_step), by = k_step)
    if (!length(values)) {
      .usdt_stop("the range of ", what, " trials holds no number that splits ",
                 "into whole numbers of signal and noise trials.")
    }
    if (length(values) > 400L) {
      values <- values[unique(round(seq(1, length(values), length.out = 400L)))]
    }
    as.integer(values)
  }
  if (linked) {
    common <- step[["direct"]] * step[["indirect"]] /
      .gcd(step[["direct"]], step[["indirect"]])
    v <- axis(ranges$direct, common, "the")
    grid <- data.frame(direct = v, indirect = v)
  } else {
    grid <- expand.grid(direct = axis(ranges$direct, step[["direct"]],
                                      "direct"),
                        indirect = axis(ranges$indirect, step[["indirect"]],
                                        "indirect"),
                        KEEP.OUT.ATTRS = FALSE)
  }
  free <- c(subjects = length(subjects) == 2L,
            direct = length(ranges$direct) == 2L,
            indirect = !linked && length(ranges$indirect) == 2L)
  list(subjects = subjects, trials = ranges, linked = linked, grid = grid,
       proportion = proportion, free = free)
}

# The smallest number of trials whose share `p` of signal trials is whole.
.trial_step <- function(p) {
  for (k in 1:1000) {
    if (abs(k * p - round(k * p)) < 1e-8) return(as.integer(k))
  }
  .usdt_stop("`signal` = ", p, " needs more than 1000 trials to give whole ",
             "numbers of signal and noise trials.")
}

# The greatest common divisor of two whole numbers.
.gcd <- function(a, b) if (b == 0) a else .gcd(b, a %% b)

# This function turns the cost argument into a function of the subjects and
# the trials of each task. Without a cost, the plan minimises its single free
# value.
.plan_cost <- function(cost, free) {
  if (is.null(cost)) {
    if (free[["subjects"]]) return(function(N, TD, TI) N)
    if (free[["direct"]]) return(function(N, TD, TI) TD)
    return(function(N, TD, TI) TI)
  }
  if (is.function(cost)) {
    return(function(N, TD, TI) {
      out <- cost(subjects = N, direct = TD, indirect = TI)
      if (!is.numeric(out) || length(out) != length(N) || anyNA(out)) {
        .usdt_stop("`cost` must return one number per design.")
      }
      out
    })
  }
  parts <- c("subject", "direct", "indirect")
  if (!is.numeric(cost) || is.null(names(cost)) || !all(names(cost) %in% parts) ||
      anyNA(cost) || any(cost < 0) || all(cost == 0)) {
    .usdt_stop("`cost` must be a function or a named vector of costs, e.g. ",
               "c(subject = 600, direct = 4, indirect = 1.5).")
  }
  unit <- stats::setNames(rep(0, 3L), parts)
  unit[names(cost)] <- cost
  function(N, TD, TI) {
    N * (unit[["subject"]] + unit[["direct"]] * TD + unit[["indirect"]] * TI)
  }
}

# This function describes a cost argument.
.cost_text <- function(cost, labels) {
  if (is.function(cost)) return("a function of the subjects and trials")
  parts <- c(subject = "subject", direct = paste(labels[["direct"]], "trial"),
             indirect = paste(labels[["indirect"]], "trial"))
  values <- vapply(cost, format, "")
  paste(paste(values, "per", parts[names(cost)]), collapse = " + ")
}

# This function gives the error variance of each task for every number of
# trials of the grid.
.plan_error_variance <- function(pop, grid, proportion, split) {
  one <- function(k, values) {
    vapply(values, function(t) {
      cells <- list(signal = c(direct = t * proportion[["direct"]],
                               indirect = t * proportion[["indirect"]]),
                    noise = c(direct = t * (1 - proportion[["direct"]]),
                              indirect = t * (1 - proportion[["indirect"]])))
      .error_variance(pop, cells, split)[[k]]
    }, 0)
  }
  direct <- unique(grid$direct)
  indirect <- unique(grid$indirect)
  cbind(direct = one("direct", direct)[match(grid$direct, direct)],
        indirect = one("indirect", indirect)[match(grid$indirect, indirect)])
}

# This function gives the closed-form standard errors of every test with one
# subject, for each design of the grid. With N subjects they shrink by
# sqrt(N).
.plan_se1 <- function(pop, v) {
  n <- nrow(v)
  s2D <- pop$sigma_D^2; s2I <- pop$sigma_I^2
  sDI <- pop$rho * pop$sigma_D * pop$sigma_I
  truth <- matrix(c(pop$mu_D, pop$mu_I, s2D, s2I, sDI), n, 5L, byrow = TRUE,
                  dimnames = list(NULL, c("gamma_D", "gamma_I", "s2_D", "s2_I",
                                          "s_DI")))
  .moment_tests(truth, s2D + v[, "direct"], s2I + v[, "indirect"], sDI, 1)
}

# This function gives the closed-form power of one hypothesis for each design,
# with N subjects each.
.plan_closed_power <- function(h, pop, se1, N, alpha) {
  test <- .power_test(h)
  centre <- switch(test,
    covariance = pop$rho * pop$sigma_D * pop$sigma_I,
    fisher = atanh(pop$rho),
    c(diff = pop$diff, intercept = pop$intercept, slope = pop$slope)[[test]])
  .power_normal(h$op, .power_bounds(h), centre, se1[[test]]$se / sqrt(N),
                alpha)
}

# This function gives the expected width of the (1 - alpha) interval of a
# quantity in usdt_tests(), for each design with N subjects: Wald, and
# Fisher-z for the correlation.
.plan_width <- function(q, pop, se1, N, alpha) {
  z <- stats::qnorm(1 - alpha / 2)
  if (q == "rho") {
    se_z <- se1$rho$se / sqrt(N) / (1 - pop$rho^2)
    return(tanh(atanh(pop$rho) + z * se_z) - tanh(atanh(pop$rho) - z * se_z))
  }
  2 * z * se1[[q]]$se / sqrt(N)
}

# This function finds, for each design, the smallest whole N in [lo, hi] at
# which `f` (a function of a vector of N, one per design) reaches `target`,
# by a bisection that runs on every design at once. Designs that do not reach
# it at hi get NA.
.min_subjects <- function(f, lo, hi, n, target) {
  a <- rep(lo, n)
  b <- rep(hi, n)
  out <- rep(NA_real_, n)
  reach <- f(b) >= target
  at_lo <- f(a) >= target
  out[at_lo] <- lo
  todo <- reach & !at_lo
  while (any(todo & b - a > 1)) {
    m <- floor((a + b) / 2)
    up <- f(m) >= target
    b <- ifelse(todo & up, m, b)
    a <- ifelse(todo & !up, m, a)
  }
  out[todo] <- b[todo]
  out
}

# This function finds the smallest whole N in [lo, hi] at which the scalar
# function `g` reaches `target`, starting from a guess: it brackets the answer
# by steps from the guess and then bisects, so it needs few evaluations when
# the guess is close.
.refine_subjects <- function(g, guess, lo, hi, target) {
  if (lo == hi) return(if (g(lo) >= target) lo else NA_real_)
  s <- min(max(round(guess), lo), hi)
  if (g(s) >= target) {
    b <- s
    a <- s
    repeat {
      if (a <= lo) return(lo)
      a <- max(lo, floor(a - max(1, 0.15 * a)))
      if (g(a) < target) break
      b <- a
    }
  } else {
    a <- s
    b <- s
    repeat {
      if (b >= hi) return(NA_real_)
      b <- min(hi, ceiling(b + max(1, 0.25 * b)))
      if (g(b) >= target) break
      a <- b
    }
  }
  while (b - a > 1) {
    m <- floor((a + b) / 2)
    if (g(m) >= target) b <- m else a <- m
  }
  b
}

# This function plans one scenario: it returns the chosen design for each
# target and the closed-form frontier.
.plan_scenario <- function(pop, mode, tests, precision, power, budget, combine,
                           alpha, split, space, cost_fun, base) {
  grid <- space$grid
  v <- .plan_error_variance(pop, grid, space$proportion, split)
  se1 <- .plan_se1(pop, v)
  lo <- min(space$subjects)
  hi <- max(space$subjects)
  n <- nrow(grid)

  # The moment engine checks a design; its chi variables are kept per N. The
  # search compares designs with the first 5000 draws, and the chosen design
  # is settled with all of them. A precision target asks the engine only how
  # often the correlation is estimable.
  engine_tests <- tests %||% list(.parse_hypothesis("rho = 0"))
  engine <- function(draws) {
    chi <- new.env(parent = emptyenv())
    function(i, N, which) {
      key <- as.character(N)
      if (is.null(chi[[key]])) chi[[key]] <- .moment_chi(draws, N)
      res <- .power_approximate(pop, N, v[i, ], engine_tests, alpha, draws,
                                chi[[key]])
      if (identical(which, "all")) c(power = res$joint,
                                     estimable = res$joint_estimable) else
        c(power = res$rows$power[which], estimable = res$rows$estimable[which])
    }
  }
  moment <- engine(base)
  attr(moment, "fast") <- engine(.moment_subset(base))

  targets <- if (mode == "precision") names(precision) else
    if (combine == "all") "all" else seq_along(tests)
  label <- function(j) {
    if (mode == "precision") return(paste0(j, " width <= ",
                                           format(precision[[j]])))
    if (identical(j, "all")) return("all at once")
    .bf_label(.bf_name(tests[[j]]$quantity), tests[[j]])
  }

  designs <- list()
  frontier <- list()
  for (j in targets) {
    # The closed form gives each design of the grid its smallest N.
    closed <- function(N) {
      if (mode == "precision") {
        return(-.plan_width(j, pop, se1, N, alpha))
      }
      if (identical(j, "all")) {
        return(do.call(pmin, lapply(tests, .plan_closed_power, pop = pop,
                                    se1 = se1, N = N, alpha = alpha)))
      }
      .plan_closed_power(tests[[j]], pop, se1, N, alpha)
    }
    goal <- if (mode == "precision") -precision[[j]] else power

    if (mode == "budget") {
      found <- .plan_budget(closed, moment, j, space, cost_fun, budget, lo, hi)
    } else {
      N0 <- .min_subjects(closed, lo, hi, n, goal)
      cost0 <- ifelse(is.na(N0), NA_real_, cost_fun(N0, grid$direct,
                                                     grid$indirect))
      frontier[[length(frontier) + 1L]] <- data.frame(
        hypothesis = label(j), trials_direct = grid$direct,
        trials_indirect = grid$indirect, subjects = N0, cost = cost0,
        stringsAsFactors = FALSE)
      found <- if (mode == "precision") {
        .plan_precision(N0, cost0, moment, grid, j, pop, se1, alpha)
      } else {
        .plan_power(N0, cost0, moment, j, space, cost_fun, lo, hi, goal)
      }
    }
    designs[[length(designs) + 1L]] <- cbind(
      data.frame(hypothesis = label(j), stringsAsFactors = FALSE), found,
      limit = .plan_limit(found, space))
  }
  list(design = do.call(rbind, designs),
       frontier = if (length(frontier)) do.call(rbind, frontier) else
         data.frame())
}

# This function searches the grid of trials for the design with the lowest
# score, checking few designs with the moment engine. `bound` gives every
# design a lower bound of its score from the closed form (NA when it cannot
# reach the target), which is optimistic, so a design whose bound is not
# below the best score so far is skipped. `score(i)` checks design i. The
# search tries the five designs with the lowest bound, then a lattice of at
# most 8 values per free axis, and then moves from the best design to better
# neighbours on the full grid, halving the step until it is one grid value;
# it checks at most `max_checks` designs, and returns those it checked with a
# finite score, best first. The lattice matters because the closed form is
# most optimistic where boundary fits are frequent, so the designs it ranks
# first need not be the best.
.plan_search <- function(bound, score, space, max_checks = 150L) {
  grid <- space$grid
  axes <- if (space$linked) list(unique(grid$direct)) else
    list(unique(grid$direct), unique(grid$indirect))
  dims <- lengths(axes)
  row_of <- function(pos) {
    if (length(dims) == 1L) pos[1L] else pos[1L] + (pos[2L] - 1L) * dims[1L]
  }
  pos_of <- function(i) {
    if (length(dims) == 1L) i else c((i - 1L) %% dims[1L] + 1L,
                                     (i - 1L) %/% dims[1L] + 1L)
  }

  scores <- rep(NA_real_, nrow(grid))
  best <- NA_integer_
  checks <- 0L
  try_row <- function(i) {
    if (!is.na(scores[i]) || is.na(bound[i]) || checks >= max_checks) return()
    if (!is.na(best) && bound[i] >= scores[best]) return()
    checks <<- checks + 1L
    scores[i] <<- score(i)
    if (is.na(best) || scores[i] < scores[best]) best <<- i
  }

  for (i in utils::head(order(bound, na.last = NA), 5L)) try_row(i)
  lattice <- as.matrix(expand.grid(lapply(dims, function(n) {
    unique(round(seq(1, n, length.out = min(n, 8L))))
  })))
  rows <- apply(lattice, 1L, row_of)
  for (i in rows[order(bound[rows], na.last = NA)]) try_row(i)

  if (!is.na(best) && is.finite(scores[best])) {
    step <- pmax(1L, ceiling(dims / 16))
    offsets <- as.matrix(expand.grid(rep(list(-1:1), length(dims))))
    offsets <- offsets[rowSums(abs(offsets)) > 0, , drop = FALSE]
    repeat {
      from <- best
      centre <- pos_of(best)
      for (k in seq_len(nrow(offsets))) {
        pos <- centre + offsets[k, ] * step
        if (all(pos >= 1L & pos <= dims)) try_row(row_of(pos))
      }
      if (best == from) {
        if (all(step == 1L) || checks >= max_checks) break
        step <- pmax(1L, step %/% 2L)
      }
    }
  }

  # The designs checked, best first, so the caller can settle the close
  # ones with more draws.
  checked <- which(is.finite(scores))
  checked[order(scores[checked])]
}

# This function finds the cheapest design that reaches the power: each design
# checked gets the smallest N that the moment engine says reaches it.
.plan_power <- function(N0, cost0, moment, j, space, cost_fun, lo, hi, goal) {
  grid <- space$grid
  fast <- attr(moment, "fast")
  subjects <- rep(NA_real_, nrow(grid))
  score <- function(i) {
    N <- .refine_subjects(function(N) fast(i, N, j)[["power"]], N0[i], lo,
                          hi, goal)
    subjects[i] <<- N
    if (is.na(N)) Inf else cost_fun(N, grid$direct[i], grid$indirect[i])
  }
  checked <- .plan_search(cost0, score, space)
  if (!length(checked)) return(.plan_none())

  # All the draws settle the subjects of the five best designs, which differ
  # by less than the Monte Carlo error of the search, and the cheapest wins.
  best <- NULL
  for (i in utils::head(checked, 5L)) {
    N <- .refine_subjects(function(N) moment(i, N, j)[["power"]], subjects[i],
                          lo, hi, goal)
    if (is.na(N)) next
    cst <- cost_fun(N, grid$direct[i], grid$indirect[i])
    if (is.null(best) || cst < best$cost) {
      m <- moment(i, N, j)
      best <- data.frame(subjects = N, trials_direct = grid$direct[i],
                         trials_indirect = grid$indirect[i],
                         power = m[["power"]], estimable = m[["estimable"]],
                         width = NA_real_, cost = cst)
    }
  }
  best %||% .plan_none()
}

# This function picks the cheapest design whose expected interval is narrow
# enough, and reports how often it is estimable.
.plan_precision <- function(N0, cost0, moment, grid, q, pop, se1, alpha) {
  if (all(is.na(cost0))) return(.plan_none())
  i <- which.min(cost0)
  m <- moment(i, N0[i], 1L)
  data.frame(subjects = N0[i], trials_direct = grid$direct[i],
             trials_indirect = grid$indirect[i], power = NA_real_,
             estimable = if (q == "diff") 1 else m[["estimable"]],
             width = .plan_width(q, pop, lapply(se1, function(e)
               list(se = e$se[i])), N0[i], alpha),
             cost = cost0[i])
}

# This function finds the design with the highest power within the budget:
# with each number of trials, the largest N it pays for. The closed-form
# power there bounds the power of the moment engine from above, so its
# negative is the bound of the search.
.plan_budget <- function(closed, moment, j, space, cost_fun, budget, lo, hi) {
  grid <- space$grid
  n <- nrow(grid)
  affordable <- function(N) -cost_fun(N, grid$direct, grid$indirect)
  N <- .min_subjects(function(N) affordable(hi - N + lo), lo, hi, n, -budget)
  N <- ifelse(is.na(N), NA_real_, hi - N + lo)
  p0 <- ifelse(is.na(N), NA_real_, closed(ifelse(is.na(N), lo, N)))
  fast <- attr(moment, "fast")
  checked <- .plan_search(-p0, function(i) -fast(i, N[i], j)[["power"]],
                          space)
  if (!length(checked)) return(.plan_none())

  # All the draws decide among the five most powerful designs found.
  power_full <- vapply(utils::head(checked, 5L), function(i) {
    moment(i, N[i], j)[["power"]]
  }, 0)
  i <- utils::head(checked, 5L)[which.max(power_full)]
  m <- moment(i, N[i], j)
  data.frame(subjects = N[i], trials_direct = grid$direct[i],
             trials_indirect = grid$indirect[i], power = m[["power"]],
             estimable = m[["estimable"]], width = NA_real_,
             cost = cost_fun(N[i], grid$direct[i], grid$indirect[i]))
}

# The row of a target no design reaches.
.plan_none <- function() {
  data.frame(subjects = NA_real_, trials_direct = NA_real_,
             trials_indirect = NA_real_, power = NA_real_,
             estimable = NA_real_, width = NA_real_, cost = NA_real_)
}

# This function names the free argument whose limit a design reaches.
.plan_limit <- function(found, space) {
  if (is.na(found$subjects)) return(NA_character_)
  at <- character(0)
  if (space$free[["subjects"]] && found$subjects %in% space$subjects) {
    at <- c(at, "`subjects`")
  }
  for (k in c("direct", "indirect")) {
    r <- space$trials[[k]]
    if (length(r) == 2L && (space$free[[k]] || (k == "indirect" && space$linked)) &&
        found[[paste0("trials_", k)]] %in% c(min(space$grid[[k]]),
                                             max(space$grid[[k]]))) {
      at <- c(at, if (space$linked) "`trials`" else paste0("`trials$", k, "`"))
    }
  }
  at <- unique(at)
  if (!length(at)) NA_character_ else paste(at, collapse = " and ")
}

# This function moves one population value of a scenario and derives the
# others again, keeping fixed what the parameterization keeps fixed.
.scenario_with <- function(pop, what, value, parameterization) {
  regression <- identical(parameterization, "regression")
  slope <- function() pop$rho * pop$sigma_I / pop$sigma_D
  switch(what,
    diff = pop$mu_I <- pop$mu_D + value,
    intercept = pop$mu_I <- value + slope() * pop$mu_D,
    slope = ,
    rho = {
      rho <- if (what == "slope") value * pop$sigma_D / pop$sigma_I else value
      if (regression) {
        pop$mu_I <- pop$intercept + rho * pop$sigma_I / pop$sigma_D * pop$mu_D
      }
      pop$rho <- rho
    },
    {
      if (what %in% c("sigma_D", "sigma_I") && regression) {
        pop[[what]] <- value
        pop$mu_I <- pop$intercept + slope() * pop$mu_D
      } else if (what == "mu_D" && regression) {
        pop$mu_D <- value
        pop$mu_I <- pop$intercept + slope() * pop$mu_D
      } else {
        pop[[what]] <- value
      }
    })
  q <- .usdt_quantities(.bivariate_primitives(cbind(
    mu_D = pop$mu_D, mu_I = pop$mu_I, sigma_D = pop$sigma_D,
    sigma_I = pop$sigma_I, rho = pop$rho)))
  pop$diff <- q[, "diff"][[1L]]
  pop$intercept <- q[, "intercept"][[1L]]
  pop$slope <- q[, "slope"][[1L]]
  pop
}

# This function finds the smallest effect of one hypothesis for one design.
# The power along the way is that of the moment engine, with the error
# variance of each moved population.
.mdes_one <- function(pop, h, N, cells, split, alpha, power, base,
                      parameterization, vary, range) {
  # The search runs on the first 5000 draws, and the power of the population
  # itself on all of them, as in usdt_power().
  engine <- function(draws) {
    chi <- .moment_chi(draws, N)
    function(p, hyp = h) {
      if (!is.finite(p$sigma_D) || p$sigma_D <= 0 || !is.finite(p$sigma_I) ||
          p$sigma_I <= 0 || abs(p$rho) >= 1) return(0)
      v <- .error_variance(p, cells, split)
      .power_approximate(p, N, v, list(hyp), alpha, draws, chi)$rows$power
    }
  }
  at <- engine(.moment_subset(base))
  values <- c(diff = pop$diff, rho = pop$rho, intercept = pop$intercept,
              slope = pop$slope)
  current <- if (is.null(vary)) values[[h$quantity]] else pop[[vary]]
  power_current <- engine(base)(pop)

  # The search runs from H0 outwards, or along the given range.
  if (!is.null(vary)) {
    varies <- vary
    from <- range[1L]
    to <- range[2L]
    f <- function(x) at(.scenario_with(pop, vary, x, parameterization))
  } else if (h$op == "in") {
    varies <- "margin"
    centre <- mean(h$value)
    from <- abs(values[[h$quantity]] - centre)
    to <- from + if (h$quantity == "rho") 2 else 3
    f <- function(x) {
      bounds <- centre + c(-1, 1) * x
      if (h$quantity == "rho") bounds <- pmin(pmax(bounds, -0.999), 0.999)
      at(pop, utils::modifyList(h, list(value = bounds)))
    }
  } else {
    varies <- h$quantity
    a <- h$value[1L]
    b <- h$value[length(h$value)]
    direction <- switch(h$op, ">" = 1, "<" = -1,
      out = if (values[[h$quantity]] < a) -1 else 1,
      if (values[[h$quantity]] < a) -1 else 1)
    from <- if (h$op == "out" && direction > 0) b else a
    limit <- if (h$quantity == "rho") 0.99 else
      if (h$quantity == "slope") 0.99 * pop$sigma_I / pop$sigma_D else NA
    to <- if (is.na(limit)) from + direction * 3 else direction * limit
    f <- function(x) at(.scenario_with(pop, h$quantity, x, parameterization))
  }

  # The curve of power over the values searched, and the first crossing.
  grid <- seq(from, to, length.out = 25L)
  curve <- vapply(grid, f, 0)
  mdes <- NA_real_
  cross <- which(curve >= power)[1L]
  if (!is.na(cross)) {
    mdes <- if (cross == 1L) grid[1L] else stats::uniroot(
      function(x) f(x) - power, grid[c(cross - 1L, cross)],
      tol = 1e-4 * max(1, abs(to - from)))$root
  }
  list(varies = varies, current = current, power_current = power_current,
       mdes = mdes, curve = data.frame(value = grid, power = curve))
}

# This function finds the smallest effects at the designs of a plan.
.mdes_plan <- function(plan, power, vary, range) {
  st <- plan$settings
  if (st$mode == "precision") {
    .usdt_stop("a precision plan has no hypotheses; give usdt_mdes() the ",
               "population, subjects and trials.")
  }
  d <- plan$design[!is.na(plan$design$subjects), , drop = FALSE]
  if (!nrow(d)) .usdt_stop("no design of the plan reaches its target.")
  out <- lapply(seq_len(nrow(d)), function(i) {
    hyp <- if (d$hypothesis[i] == "all") st$hypothesis else
      st$hypothesis[match(d$hypothesis[i], vapply(st$tests, function(h)
        .bf_label(.bf_name(h$quantity), h), ""))]
    usdt_mdes(.population_rows(st$population, d$scenario[i]),
              subjects = d$subjects[i],
              trials = c(direct = d$trials_direct[i],
                         indirect = d$trials_indirect[i]),
              hypothesis = hyp, power = power, alpha = st$alpha,
              signal = st$signal, dichotomize = st$split, vary = vary,
              range = range)
  })
  .bind_results(out, "usdt_mdes")
}

# This function keeps some scenarios of a population as a population.
.population_rows <- function(population, scenarios) {
  keep <- as.data.frame(population)$scenario %in% scenarios
  out <- as.data.frame(population)[keep, , drop = FALSE]
  rownames(out) <- NULL
  structure(out, class = c("usdt_population", "data.frame"),
            parameterization = attr(population, "parameterization"),
            labels = attr(population, "labels"))
}

# This function binds results of usdt_power() or usdt_mdes() for several
# designs into one object: the rows, and the attributes that hold one element
# per design.
.bind_results <- function(parts, class) {
  first <- parts[[1L]]
  out <- do.call(rbind, lapply(parts, as.data.frame))
  rownames(out) <- NULL
  keep <- attributes(first)
  keep$names <- keep$row.names <- NULL
  for (a in intersect(c("replicates", "curve"), names(keep))) {
    keep[[a]] <- do.call(c, lapply(parts, attr, a))
  }
  if ("designs" %in% names(keep)) {
    keep$designs <- do.call(rbind, lapply(parts, attr, "designs"))
  }
  if ("population" %in% names(keep)) {
    pops <- do.call(rbind, lapply(parts, attr, "population"))
    keep$population <- pops[!duplicated(pops$scenario), , drop = FALSE]
  }
  attributes(out) <- c(attributes(out)[c("names", "row.names")], keep)
  class(out) <- c(class, "data.frame")
  out
}
