# data.R
# This script prepares direct and indirect task data for uSDT models.
# Author: Ricardo Rey-Sáez
# Last modified: 08-09-2026

# Public functions

#' Prepare the data of a direct and an indirect task
#'
#' Both functions build the same object. Use `usdt_data_tasks()` when each task
#' has its own data frame, and `usdt_data_long()` when a single data frame
#' holds both tasks together with a column that identifies them.
#'
#' The functions count the responses of every subject in each condition, check
#' that the two tasks describe the same subjects, and record how each column
#' was read. Printing the result shows all of that, so the coding can be
#' checked before the model runs.
#'
#' @param direct,indirect The data frame of each task, for `usdt_data_tasks()`.
#'   `usdt_data_long()` ignores them and uses `task_levels` instead.
#' @param data A single data frame holding both tasks, for `usdt_data_long()`.
#' @param task_col Name of the column that identifies the task, for
#'   `usdt_data_long()`.
#' @param task_levels Which value of `task_col` belongs to each task, as
#'   `c(direct = "D", indirect = "I")`.
#' @param subject_col Name of the column that identifies the subject.
#' @param condition_col Name of the column that holds the signal and noise
#'   condition. It is not needed when the data arrive as an SDT table through
#'   `sdt_cols`.
#' @param condition_levels Which value of `condition_col` plays each role, as
#'   `c(signal = "old", noise = "new")`. The function guesses them and reports
#'   its choice when they are missing.
#' @param response_col Name of the column that holds the response. It can be a
#'   binary response, or a continuous measure such as response times when the
#'   task is dichotomized.
#' @param response_levels Which value of `response_col` counts as a signal
#'   response, as `c(signal = 1, noise = 0)`. A task that is dichotomized takes
#'   the side of the median instead, as
#'   `c(signal = "faster", noise = "slower")`. The function guesses them and
#'   reports its choice when they are missing.
#' @param successes_col,trials_col Names of the columns that hold data already
#'   summed up, the number or proportion of signal responses and the number of
#'   trials. Give these instead of `response_col`.
#' @param successes_type Format of `successes_col`. Use `"counts"` for counts
#'   and `"proportions"` for proportions. `"auto"` recognises clear cases and
#'   asks for an explicit choice when every value is zero or one.
#' @param sdt_cols Names of the columns of an SDT table, as
#'   `c(hit = "H", miss = "M", fa = "FA", cr = "CR")`. Give these instead of
#'   `condition_col` and `response_col`.
#' @param dichotomize Which tasks need the median split of [meyen_split()].
#'   Name them with `"none"`, `"direct"`, `"indirect"` or `"both"`, or give one
#'   logical value per task, as `list(direct = FALSE, indirect = TRUE)`.
#' @param ties What to do with trials that fall exactly on the median. See
#'   [meyen_split()].
#' @param coding How the condition enters the model. See Details.
#' @param labels Display names for the two tasks. They only affect printed
#'   output.
#'
#' @return An object of class `usdt_data`. Its `agg` element is the data frame
#'   of counts that the model uses, with one row per subject, task and
#'   condition. Its `meta` element records how every column was read, which
#'   tasks were split at the median, and the descriptive summaries shown when
#'   the object is printed. Pass the object to [hsdt()].
#'
#' @details
#' # One value or one per task
#'
#' The two tasks rarely come from the same experimental design, so every
#' argument that names a column, a level or a format accepts two forms. Give
#' one value and both tasks use it. Give one value per task and each task is
#' read on its own.
#'
#' ```
#' subject_col      = "subj"                       # both tasks
#' condition_col    = list(direct   = "condition",
#'                         indirect = "cue")       # one per task
#' condition_levels = list(
#'   direct   = c(signal = "old",  noise = "new"),
#'   indirect = c(signal = "cued", noise = "uncued"))
#' dichotomize      = list(direct = FALSE, indirect = TRUE)
#' ```
#'
#' This applies to `subject_col`, `condition_col`, `condition_levels`,
#' `response_col`, `response_levels`, `successes_col`, `trials_col`,
#' `successes_type`, `sdt_cols`, `dichotomize` and `ties`. The two tasks may
#' therefore use different columns, different values inside those columns, and
#' even different formats, with one task given trial by trial and the other as
#' an SDT table. Use `list()` rather than `c()` when the value of a task is
#' itself a vector, as happens with the `*_levels` and `sdt_cols` arguments.
#'
#' Only `coding` works differently. It describes the model itself, so it always
#' applies to both tasks at once.
#'
#' # Condition coding
#'
#' The two codings answer different questions and give different intercepts.
#' Under `"deviation"` the condition takes the values -0.5 and +0.5, and the
#' intercept is `-c`, the criterion measured from the midpoint between the
#' signal and noise distributions. Under `"treatment"` the condition takes the
#' values 0 and 1, and the intercept is `z(FAR)`, the criterion measured from
#' the noise distribution. The two intercepts are related by
#' `intercept_treatment = intercept_deviation - d'/2`.
#'
#' The choice matters for a task that was split at the median. The split leaves
#' each subject with half signal responses, so with balanced conditions the
#' deviation intercept is exactly zero and needs no estimation. The treatment
#' intercept equals `-d'/2` instead, which is not zero and has to be estimated.
#' The object records this in `criterion_zero`, and [hsdt()] uses it to decide.
#'
#' @seealso [meyen_split()], [hsdt()], [sdt_moments()]
#'
#' @examples
#' set.seed(1)
#' df <- usdt_simulate(n_subj = 30, n_trials = 80)
#'
#' # Both tasks share every column name and every level here.
#' d  <- usdt_data_long(df, task_col = "task",
#'                      task_levels   = c(direct = "D", indirect = "I"),
#'                      subject_col   = "subj",
#'                      condition_col = "cond",
#'                      condition_levels = c(signal = 1, noise = 0),
#'                      response_col  = "response",
#'                      response_levels  = c(signal = 1, noise = 0))
#' d
#'
#' # When they do not, each task gets its own column and its own levels.
#' aware <- df[df$task == "D", ]
#' cuing <- df[df$task == "I", ]
#' names(aware)[names(aware) == "cond"] <- "seen"
#' names(cuing)[names(cuing) == "cond"] <- "cue"
#' aware$seen <- ifelse(aware$seen == 1, "old",  "new")
#' cuing$cue  <- ifelse(cuing$cue  == 1, "cued", "uncued")
#'
#' usdt_data_tasks(
#'   direct = aware, indirect = cuing,
#'   subject_col      = "subj",
#'   condition_col    = list(direct = "seen", indirect = "cue"),
#'   condition_levels = list(direct   = c(signal = "old",  noise = "new"),
#'                           indirect = c(signal = "cued", noise = "uncued")),
#'   response_col     = "response",
#'   response_levels  = c(signal = 1, noise = 0))
#'
#' @name usdt_data
NULL

#' @rdname usdt_data
#' @export
usdt_data_tasks <- function(direct, indirect,
                            subject_col,
                            condition_col    = NULL, condition_levels = NULL,
                            response_col     = NULL, response_levels  = NULL,
                            successes_col    = NULL, trials_col       = NULL,
                            successes_type   = c("auto", "counts", "proportions"),
                            sdt_cols         = NULL,
                            dichotomize      = c("none", "indirect", "direct", "both"),
                            ties             = c("noise", "random"),
                            coding           = c("deviation", "treatment"),
                            labels           = c(direct = "Direct", indirect = "Indirect")) {

  # Both task inputs must contain data.
  .check_df(direct,   "direct")
  .check_df(indirect, "indirect")

  # The shared builder prepares both tasks.
  .usdt_build(parts       = list(direct = direct, indirect = indirect),
              input       = "2 data frames",
              entry       = "usdt_data_tasks()",
              subject_col = subject_col,
              condition_col = condition_col, condition_levels = condition_levels,
              response_col  = response_col,  response_levels  = response_levels,
              successes_col = successes_col, trials_col       = trials_col,
              successes_type = successes_type,
              sdt_cols      = sdt_cols,
              dichotomize   = dichotomize,
              ties          = ties,
              coding        = match.arg(coding),
              labels        = labels)
}

#' @rdname usdt_data
#' @export
usdt_data_long <- function(data, task_col, task_levels,
                           subject_col,
                           condition_col    = NULL, condition_levels = NULL,
                           response_col     = NULL, response_levels  = NULL,
                           successes_col    = NULL, trials_col       = NULL,
                           successes_type   = c("auto", "counts", "proportions"),
                           sdt_cols         = NULL,
                           dichotomize      = c("none", "indirect", "direct", "both"),
                           ties             = c("noise", "random"),
                           coding           = c("deviation", "treatment"),
                           labels           = NULL) {

  # The input must contain a valid task column.
  .check_df(data, "data")
  .check_string(task_col, "task_col")
  .check_cols(data, task_col, "task_col", "`data`")
  task_levels <- .check_task_levels(task_levels, data[[task_col]], task_col)

  # The task labels separate the direct and indirect data.
  tv <- as.character(data[[task_col]])
  parts <- list(direct   = data[tv == as.character(task_levels[["direct"]]),   , drop = FALSE],
                indirect = data[tv == as.character(task_levels[["indirect"]]), , drop = FALSE])
  for (k in names(parts)) {
    if (nrow(parts[[k]]) == 0L) {
      .usdt_stop("no rows of `", task_col, "` equal `", task_levels[[k]],
                 "`, the value given for the ", k, " task.")
    }
  }

  # The shared builder prepares both tasks.
  .usdt_build(parts       = parts,
              input       = "1 long data frame",
              entry       = "usdt_data_long()",
              subject_col = subject_col,
              condition_col = condition_col, condition_levels = condition_levels,
              response_col  = response_col,  response_levels  = response_levels,
              successes_col = successes_col, trials_col       = trials_col,
              successes_type = successes_type,
              sdt_cols      = sdt_cols,
              dichotomize   = dichotomize,
              ties          = ties,
              coding        = match.arg(coding),
              labels        = if (is.null(labels))
                c(direct = "Direct", indirect = "Indirect") else labels)
}

# Data builder

# This function prepares and aggregates both task inputs.
.usdt_build <- function(parts, input, entry, subject_col,
                        condition_col, condition_levels,
                        response_col,  response_levels,
                        successes_col, trials_col, successes_type, sdt_cols,
                        dichotomize, ties, coding, labels) {

  labels <- .check_labels(labels)

  # Each argument receives one value for each task.
  arg <- list(
    subject   = .per_task(subject_col,      "subject_col"),
    condition = .per_task(condition_col,    "condition_col",    allow_null = TRUE),
    cond_lev  = .per_task(condition_levels, "condition_levels", allow_null = TRUE),
    response  = .per_task(response_col,     "response_col",     allow_null = TRUE),
    resp_lev  = .per_task(response_levels,  "response_levels",  allow_null = TRUE),
    successes = .per_task(successes_col,    "successes_col",    allow_null = TRUE),
    trials    = .per_task(trials_col,       "trials_col",       allow_null = TRUE),
    sdt       = .per_task(sdt_cols,         "sdt_cols",         allow_null = TRUE),
    stype     = .per_task_choice(successes_type, "successes_type",
                                 c("auto", "counts", "proportions")),
    ties      = .per_task_choice(ties, "ties", c("noise", "random"))
  )
  dic <- .check_dichotomize(dichotomize)

  # Each task becomes a table of response counts.
  cells <- list(); info <- list()
  for (k in c("direct", "indirect")) {
    res <- .prepare_task(
      data      = parts[[k]],          task     = k,
      subject   = arg$subject[[k]],    condition = arg$condition[[k]],
      cond_lev  = arg$cond_lev[[k]],   response  = arg$response[[k]],
      resp_lev  = arg$resp_lev[[k]],   successes = arg$successes[[k]],
      trials    = arg$trials[[k]],     sdt       = arg$sdt[[k]],
      successes_type = arg$stype[[k]],
      dichotomize = dic[[k]],          ties      = arg$ties[[k]]
    )
    cells[[k]] <- res$cells
    info[[k]]  <- res$info
  }

  # The code joins both tasks and creates the model columns.
  agg <- rbind(
    data.frame(subj = cells$direct$subj,   task = "D", sig = cells$direct$sig,
               y = cells$direct$y,   n = cells$direct$n,   stringsAsFactors = FALSE),
    data.frame(subj = cells$indirect$subj, task = "I", sig = cells$indirect$sig,
               y = cells$indirect$y, n = cells$indirect$n, stringsAsFactors = FALSE)
  )
  agg$task <- factor(agg$task, levels = c("D", "I"))
  agg$cond <- if (coding == "deviation") ifelse(agg$sig, .5, -.5) else
    as.numeric(agg$sig)
  agg$c_D <- as.integer(agg$task == "D")
  agg$c_I <- as.integer(agg$task == "I")
  agg$d_D <- agg$cond * (agg$task == "D")
  agg$d_I <- agg$cond * (agg$task == "I")
  agg <- agg[order(agg$subj, agg$task, agg$sig), ]
  rownames(agg) <- NULL

  # The code checks which subjects completed both tasks.
  sD <- unique(agg$subj[agg$task == "D"]); sI <- unique(agg$subj[agg$task == "I"])
  both <- intersect(sD, sI); only <- setdiff(union(sD, sI), both)
  if (!length(both)) {
    .usdt_stop("no subject appears in both tasks. The latent correlation ",
               "is not identified.\n  Direct task ids: ",
               paste(utils::head(sD, 5L), collapse = ", "), ", ...\n",
               "  Indirect task ids: ",
               paste(utils::head(sI, 5L), collapse = ", "), ", ...\n",
               "  Do the two data frames use the same subject labels?")
  }
  if (length(only)) {
    .usdt_warn(length(only), " subject", if (length(only) > 1L) "s" else "",
               " appear in only one task. They contribute to the group means ",
               "but not to the latent correlation.")
  }
  if (length(both) < 20L) {
    .usdt_warn("only ", length(both), " subjects in both tasks. The latent ",
               "correlation and regression will be very imprecise.")
  }

  # The code checks the criterion left by each Meyen split.
  chk   <- .criterion_check(agg, coding)
  czero <- stats::setNames(logical(2), c("direct", "indirect"))
  for (k in c("direct", "indirect")) {
    czero[[k]] <- info[[k]]$dichotomized && chk[[k]]$negligible
    if (info[[k]]$dichotomized && !chk[[k]]$negligible) {
      .usdt_msg("task ", labels[[k]], " was Meyen-split but its implied ",
                "criterion averages ", signif(chk[[k]]$mean_abs, 3), ", above the ",
                "tolerance of ", chk[[k]]$tol, ", so it will be estimated.",
                if (coding == "treatment")
                  " Under treatment coding that criterion is -d'/2, never zero; use `coding = \"deviation\"` to drop it."
                else "")
    }
  }

  # The result stores the data and its preparation details.
  structure(list(
    agg  = agg,
    meta = list(
      input = input, entry = entry, coding = coding,
      ties = unlist(arg$ties)[c("direct", "indirect")],
      labels = labels, tasks = info, criterion_zero = czero, criterion = chk,
      n_subj = length(union(sD, sI)), n_both = length(both), n_only = length(only),
      n_trials = sum(agg$n), n_rows = nrow(agg),
      descriptives = .descriptives(agg, coding)
    )
  ), class = "usdt_data")
}

# This function prepares one task for the model.
.prepare_task <- function(data, task, subject, condition, cond_lev,
                          response, resp_lev, successes, trials, sdt,
                          successes_type, dichotomize, ties) {

  where <- paste0("the ", task, " task")
  .check_string(subject, "subject_col")
  .check_cols(data, subject, "subject_col", where)
  .check_ids(data[[subject]], subject, where)
  subj <- as.character(data[[subject]])

  # The supplied columns reveal the input format.
  gran <- if (!is.null(sdt)) "sdt" else
    if (!is.null(successes) && !is.null(trials)) "counts" else
      if (!is.null(response)) "trials" else
        .usdt_stop("for ", where, " supply one of `response_col` (trial level), ",
                   "`successes_col` + `trials_col` (counts or proportions), or ",
                   "`sdt_cols` (an SDT table).")
  if (dichotomize && gran != "trials") {
    .usdt_stop("`dichotomize` was requested for ", where, ", but its data are ",
               "aggregated (", gran, "). The Meyen median split needs one row ",
               "per trial.")
  }

  # Every input format becomes signal and noise cells.
  out <- switch(gran,
    sdt    = .cells_from_sdt(data, subj, sdt, where),
    counts = .cells_from_counts(data, subj, condition, cond_lev, successes,
                                trials, successes_type, where),
    trials = .cells_from_trials(data, subj, condition, cond_lev, response,
                                resp_lev, dichotomize, ties, where)
  )

  # The result records the subject column.
  out$info$subject <- subject

  # Every subject must have trials in both conditions.
  if (any(out$cells$n == 0L)) {
    .usdt_stop("some cells of ", where, " contain no trials.")
  }
  keys <- table(out$cells$subj)
  if (any(keys != 2L)) {
    bad <- names(keys)[keys != 2L]
    .usdt_stop(length(bad), " subject", if (length(bad) > 1L) "s" else "",
               " of ", where, " lack one of the two conditions (e.g. `",
               bad[1L], "`). Both signal and noise trials are required.")
  }
  out
}

# This function prepares trial-level data.
.cells_from_trials <- function(data, subj, condition, cond_lev, response,
                               resp_lev, dichotomize, ties, where) {

  # The input must contain condition and response columns.
  .check_string(condition, "condition_col"); .check_string(response, "response_col")
  .check_cols(data, c(condition, response), "condition_col/response_col", where)
  cv <- data[[condition]]; rv <- data[[response]]

  # The condition values receive signal and noise roles.
  cond_lev <- .check_levels(cond_lev, "condition_levels")
  if (is.null(cond_lev)) {
    cond_lev <- .guess_levels(cv, "condition_levels", paste0("`", condition, "` in ", where))
  }
  sig <- .role_match(cv, cond_lev, condition, where)

  # The response becomes a binary signal choice.
  if (dichotomize) {
    .check_numeric(rv, response, where, allow_na = TRUE)
    side <- .split_side(resp_lev, response, where)
    resp <- meyen_split(rv, by = subj, signal = side, ties = ties)
    dic  <- list(side = side, p = tapply(resp, subj, mean, na.rm = TRUE),
                 odd = sum(table(subj) %% 2L == 1L))
    resp_lev <- stats::setNames(
      c(side, if (side == "faster") "slower" else "faster"),
      c("signal", "noise"))
  } else {
    resp_lev <- .check_levels(resp_lev, "response_levels")
    if (is.null(resp_lev)) {
      resp_lev <- .guess_levels(rv, "response_levels", paste0("`", response, "` in ", where))
    }
    resp <- as.integer(.role_match(rv, resp_lev, response, where))
    dic  <- NULL
  }

  # Complete trials become two condition rows for each subject.
  ok <- !is.na(resp) & !is.na(sig)
  if (any(!ok)) {
    .usdt_warn(sum(!ok), " trials of ", where, " dropped for missing values.")
  }
  cells <- .agg_counts(resp[ok], by = list(subj = subj[ok], sig = sig[ok]))
  list(cells = cells,
       info  = list(granularity = "trial level", dichotomized = dichotomize,
                    condition_col = condition, response_col = response,
                    condition_levels = cond_lev, response_levels = resp_lev,
                    n_trials = sum(ok), dic = dic))
}

# This function prepares counts or proportions.
.cells_from_counts <- function(data, subj, condition, cond_lev, successes,
                               trials, successes_type, where) {

  # The input columns and condition roles are checked first.
  .check_string(condition, "condition_col")
  .check_string(successes, "successes_col")
  .check_string(trials, "trials_col")
  .check_cols(data, c(condition, successes, trials),
              "condition_col/successes_col/trials_col", where)
  cond_lev <- .check_levels(cond_lev, "condition_levels")
  if (is.null(cond_lev)) {
    cond_lev <- .guess_levels(data[[condition]], "condition_levels",
                              paste0("`", condition, "` in ", where))
  }
  sig <- .role_match(data[[condition]], cond_lev, condition, where)
  s <- data[[successes]]
  n <- data[[trials]]
  .check_numeric(n, trials, where, positive = TRUE, whole = TRUE)
  .check_numeric(s, successes, where, nonnegative = TRUE)

  # The input format must be clear before values become counts.
  is_prop <- successes_type == "proportions" ||
    (successes_type == "auto" && any(abs(s - round(s)) > 1e-8))
  ambiguous <- successes_type == "auto" && all(s %in% c(0, 1)) && any(n > 1)
  if (ambiguous) {
    .usdt_stop("`", successes, "` in ", where, " contains only zero and one, ",
               "so its format is ambiguous. Set `successes_type` to `counts` ",
               "or `proportions`.")
  }
  if (is_prop) {
    if (any(s > 1)) {
      .usdt_stop("`", successes, "` in ", where,
                 " contains proportions outside the range from zero to one.")
    }
    y <- round(s * n)
    if (any(abs(s * n - y) > 1e-8)) {
      .usdt_stop("some values of `", successes, " * ", trials, "` in ", where,
                 " are not whole counts.")
    }
  } else {
    .check_numeric(s, successes, where, nonnegative = TRUE, whole = TRUE)
    y <- s
    if (any(y > n)) {
      .usdt_stop("in ", where, ", `", successes, "` exceeds `", trials,
                 "` in ", sum(y > n), " rows.")
    }
  }

  # Extra rows are combined within each subject and condition.
  cells <- .agg_sums(y, n, list(subj = subj, sig = sig))
  list(cells = cells,
       info  = list(granularity = if (is_prop) "proportions" else "counts",
                    dichotomized = FALSE, condition_col = condition,
                    response_col = successes, condition_levels = cond_lev,
                    response_levels = NULL, n_trials = sum(cells$n), dic = NULL))
}

# This function prepares a table with the four response counts.
.cells_from_sdt <- function(data, subj, sdt, where) {

  # The table must identify all four response counts.
  if (!setequal(names(sdt), c("hit", "miss", "fa", "cr"))) {
    .usdt_stop("`sdt_cols` must be named `hit`, `miss`, `fa` and `cr`, e.g.\n",
               '  sdt_cols = c(hit = "H", miss = "M", fa = "FA", cr = "CR")')
  }
  if (anyDuplicated(unname(sdt))) {
    .usdt_stop("`sdt_cols` must use a different column for each response count.")
  }
  .check_cols(data, unlist(sdt), "sdt_cols", where)
  g <- lapply(sdt[c("hit", "miss", "fa", "cr")], function(cl) {
    .check_numeric(data[[cl]], cl, where, nonnegative = TRUE, whole = TRUE)
    data[[cl]]
  })

  # Hits form the signal count and false alarms form the noise count.
  cells <- .agg_sums(c(g$hit, g$fa),
                     c(g$hit + g$miss, g$fa + g$cr),
                     list(subj = rep(subj, 2L),
                          sig  = rep(c(TRUE, FALSE), each = length(subj))))
  list(cells = cells,
       info  = list(granularity = "SDT table", dichotomized = FALSE,
                    condition_col = NA_character_,
                    response_col = paste(unlist(sdt), collapse = "/"),
                    condition_levels = NULL, response_levels = NULL,
                    n_trials = sum(cells$n), dic = NULL))
}

# Internal helpers

# This function maps raw values to signal and noise roles.
.role_match <- function(v, lev, col, where) {
  vs <- as.character(v)
  out <- ifelse(vs == as.character(lev[["signal"]]), TRUE,
                ifelse(vs == as.character(lev[["noise"]]), FALSE, NA))
  bad <- is.na(out) & !is.na(vs)
  if (any(bad)) {
    extra <- unique(vs[bad])
    .usdt_stop("`", col, "` in ", where, " has ", sum(bad), " values that are ",
               "neither signal (`", lev[["signal"]], "`) nor noise (`",
               lev[["noise"]], "`): ",
               paste0("`", utils::head(extra, 5L), "`", collapse = ", "),
               if (length(extra) > 5L) ", ..." else "", ".")
  }
  out
}

# This function identifies which side of the median means signal.
.split_side <- function(resp_lev, col, where) {
  if (is.null(resp_lev)) {
    .usdt_msg("`", col, "` in ", where, " is being median-split; taking ",
              "faster responses as signal. Set `response_levels = ",
              'c(signal = "slower", noise = "faster")` to reverse it.')
    return("faster")
  }
  resp_lev <- .check_levels(resp_lev, "response_levels")
  side <- as.character(resp_lev[["signal"]])
  if (!side %in% c("faster", "slower")) {
    .usdt_stop("`response_levels` for a dichotomized task must use the side of ",
               "the median, not a data value: c(signal = \"faster\", noise = ",
               "\"slower\") or the reverse. Got `", side, "`.")
  }
  side
}

# This function combines existing counts within each group.
.agg_sums <- function(y, n, by) {
  key <- interaction(by, drop = TRUE, sep = "\r")
  ys <- rowsum(as.numeric(y), key, reorder = TRUE)
  ns <- rowsum(as.numeric(n), key, reorder = TRUE)
  parts <- do.call(rbind, strsplit(rownames(ys), "\r", fixed = TRUE))
  out <- as.data.frame(parts, stringsAsFactors = FALSE)
  names(out) <- names(by)
  out$sig <- as.logical(out$sig)
  out$y <- as.integer(round(ys[, 1L]))
  out$n <- as.integer(round(ns[, 1L]))
  rownames(out) <- NULL
  out
}

# This function checks whether each task has a negligible criterion.
.criterion_check <- function(agg, coding, tol = 0.02) {
  out <- list()
  for (k in c("direct", "indirect")) {
    crit <- .task_moments(agg, if (k == "direct") "D" else "I",
                          coding)$criterion
    mean_abs <- mean(abs(crit))
    out[[k]] <- list(mean_abs = mean_abs, max_abs = max(abs(crit)),
                     sd = stats::sd(crit), negligible = mean_abs < tol,
                     tol = tol)
  }
  out
}

# This function calculates simple summaries for each subject.
.descriptives <- function(agg, coding) {
  out <- list()
  for (k in c("D", "I")) {
    m <- .task_moments(agg, k, coding)
    out[[k]] <- list(
      trials = stats::median(agg$n[agg$task == k]),
      hr = m$hr, far = m$far, dprime = m$dprime, edge = sum(m$corrected)
    )
  }
  out
}

# This function checks the direct and indirect task labels.
.check_task_levels <- function(x, v, col) {
  if (is.null(x) || length(x) != 2L || !setequal(names(x), c("direct", "indirect"))) {
    u <- unique(as.character(v))
    .usdt_stop("`task_levels` must be named `direct` and `indirect`, e.g.\n",
               '  task_levels = c(direct = "D", indirect = "I")\n',
               "  Values found in `", col, "`: ",
               paste0("`", utils::head(u, 6L), "`", collapse = ", "))
  }
  if (as.character(x[["direct"]]) == as.character(x[["indirect"]])) {
    .usdt_stop("`task_levels` must use different values for the direct and ",
               "indirect tasks.")
  }
  x
}
