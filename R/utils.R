# utils.R
# This script provides shared tools for the package.
# Author: Ricardo Rey-Sáez
# Last modified: 08-09-2026

# Messages

# This function stops the analysis with an error.
.usdt_stop <- function(...) {
  stop(paste0(...), call. = FALSE)
}

# This function shows a warning.
.usdt_warn <- function(...) {
  warning(paste0(...), call. = FALSE)
}

# This function explains a choice made by the package.
.usdt_msg <- function(...) {
  message(paste0("uSDT: ", ...))
}

# Input checks

# This function checks that a data frame contains rows.
.check_df <- function(x, arg) {
  if (!is.data.frame(x)) {
    .usdt_stop("`", arg, "` must be a data frame, not ", class(x)[1], ".")
  }
  if (nrow(x) == 0L) .usdt_stop("`", arg, "` has no rows.")
  invisible(TRUE)
}

# This function reports missing columns.
.check_cols <- function(data, cols, arg, where) {
  cols <- cols[!is.na(cols)]
  missing <- setdiff(cols, names(data))
  if (length(missing)) {
    .usdt_stop("column", if (length(missing) > 1L) "s" else "", " ",
               paste0("`", missing, "`", collapse = ", "),
               " (from `", arg, "`) not found in ", where, ".\n",
               "  Available: ", paste(utils::head(names(data), 15L), collapse = ", "),
               if (length(names(data)) > 15L) ", ..." else "")
  }
  invisible(TRUE)
}

# This function checks one column name.
.check_string <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    .usdt_stop("`", arg, "` must be a single column name (a character string).")
  }
  invisible(TRUE)
}

# This function checks subject identifiers.
.check_ids <- function(x, col, where) {
  missing <- is.na(x) | !nzchar(trimws(as.character(x)))
  if (any(missing)) {
    .usdt_stop("`", col, "` in ", where, " contains ", sum(missing),
               " missing subject identifier",
               if (sum(missing) == 1L) "." else "s.")
  }
  invisible(TRUE)
}

# This function checks numeric values.
.check_numeric <- function(x, col, where, allow_na = FALSE,
                           nonnegative = FALSE, positive = FALSE,
                           whole = FALSE) {
  if (!is.numeric(x)) {
    .usdt_stop("`", col, "` in ", where, " must be numeric, not ",
               class(x)[1L], ".")
  }
  if (!allow_na && anyNA(x)) {
    .usdt_stop("`", col, "` in ", where, " contains missing values.")
  }
  observed <- x[!is.na(x)]
  if (any(!is.finite(observed))) {
    .usdt_stop("`", col, "` in ", where, " contains infinite values.")
  }
  if (nonnegative && any(observed < 0)) {
    .usdt_stop("`", col, "` in ", where, " contains negative values.")
  }
  if (positive && any(observed <= 0)) {
    .usdt_stop("`", col, "` in ", where, " must contain positive values.")
  }
  if (whole && any(abs(observed - round(observed)) > 1e-8)) {
    .usdt_stop("`", col, "` in ", where, " must contain whole numbers.")
  }
  invisible(TRUE)
}

# This function checks one number against a range.
.check_scalar_number <- function(x, arg, lower = -Inf, upper = Inf,
                                 open_lower = FALSE, open_upper = FALSE,
                                 whole = FALSE) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
  if (valid) {
    valid <- if (open_lower) x > lower else x >= lower
    valid <- valid && if (open_upper) x < upper else x <= upper
    valid <- valid && (!whole || abs(x - round(x)) <= 1e-8)
  }
  if (!valid) .usdt_stop("`", arg, "` has an invalid value.")
  invisible(TRUE)
}

# This function checks a confidence level.
.check_confidence_level <- function(x) {
  .check_scalar_number(x, "level", lower = 0, upper = 1,
                       open_lower = TRUE, open_upper = TRUE)
}

# This function checks two display labels.
.check_labels <- function(x) {
  if (!is.character(x) || length(x) != 2L ||
      !setequal(names(x), c("direct", "indirect")) || anyNA(x) ||
      any(!nzchar(trimws(x))) || x[["direct"]] == x[["indirect"]]) {
    .usdt_stop("`labels` must contain two different names called `direct` ",
               "and `indirect`.")
  }
  x[c("direct", "indirect")]
}

# Argument handling

# This function gives one value to each task.
.per_task <- function(x, arg, allow_null = FALSE) {

  # The function accepts a missing value only when the caller allows it.
  if (is.null(x)) {
    if (allow_null) return(list(direct = NULL, indirect = NULL))
    .usdt_stop("`", arg, "` is required.")
  }

  # Task names provide separate values.
  nms <- names(x)
  if (length(x) == 2L && !is.null(nms) && setequal(nms, c("direct", "indirect"))) {
    return(list(direct = x[["direct"]], indirect = x[["indirect"]]))
  }

  # A partial task name is a typo, not a value shared by both tasks.
  if (any(c("direct", "indirect") %in% nms)) {
    .usdt_stop("`", arg, "` names one task but not the other. Give one value ",
               "for both tasks, or one for each, e.g.\n  ", arg,
               " = list(direct = ..., indirect = ...)\n  Names found: ",
               paste0("`", nms, "`", collapse = ", "), ".")
  }

  # Both tasks use the same value in every other case.
  list(direct = x, indirect = x)
}

# This function gives each task one option from a fixed set.
.per_task_choice <- function(x, arg, choices) {
  v <- .per_task(x, arg)
  lapply(v, function(z) {
    if (!is.character(z) || anyNA(z)) {
      .usdt_stop("`", arg, "` must be one of ",
                 paste0("`", choices, "`", collapse = ", "), ".")
    }
    tryCatch(match.arg(z, choices),
             error = function(e)
               .usdt_stop("`", arg, "` must be one of ",
                          paste0("`", choices, "`", collapse = ", "),
                          ", not `", paste(z, collapse = "`, `"), "`."))
  })
}

# This function decides which tasks need the median split.
.check_dichotomize <- function(x) {

  tasks <- c("direct", "indirect")
  out   <- stats::setNames(logical(2L), tasks)

  # A logical value for each task states the choice directly.
  if (is.logical(x) || (is.list(x) && all(vapply(x, is.logical, TRUE)))) {
    v <- .per_task(x, "dichotomize")
    for (k in tasks) {
      if (length(v[[k]]) != 1L || is.na(v[[k]])) {
        .usdt_stop("`dichotomize` must give one `TRUE` or `FALSE` to each task.")
      }
      out[[k]] <- as.logical(v[[k]])
    }
    return(out)
  }

  # Task names select the tasks to split.
  if (!is.character(x) || !length(x) || anyNA(x)) {
    .usdt_stop("`dichotomize` must be `\"none\"`, `\"direct\"`, `\"indirect\"`, ",
               "`\"both\"`, or one logical value for each task, e.g.\n",
               "  dichotomize = list(direct = FALSE, indirect = TRUE)")
  }
  choices <- c("none", "indirect", "direct", "both")
  if (identical(x, choices)) x <- "none"
  bad <- setdiff(x, choices)
  if (length(bad)) {
    .usdt_stop("`dichotomize` does not accept `", bad[1L], "`. Use `\"none\"`, ",
               "`\"direct\"`, `\"indirect\"`, `\"both\"`, or one logical value ",
               "for each task.")
  }
  if (length(x) > 1L && any(c("none", "both") %in% x)) {
    .usdt_stop("`dichotomize` combines `", paste(x, collapse = "`, `"),
               "`. Name the tasks to split, or use `\"none\"` or `\"both\"`.")
  }
  which_dic <- if ("both" %in% x) tasks else setdiff(x, "none")
  out[which_dic] <- TRUE
  out
}

# This function checks the signal and noise values.
.check_levels <- function(x, arg) {
  if (is.null(x)) return(NULL)
  if (length(x) != 2L || !setequal(names(x), c("signal", "noise"))) {
    .usdt_stop("`", arg, "` must be a length-2 vector named `signal` and ",
               "`noise`, e.g.\n  ", arg, ' = c(signal = "old", noise = "new")')
  }
  if (identical(as.character(x[["signal"]]), as.character(x[["noise"]]))) {
    .usdt_stop("`", arg, "` gives the same value for `signal` and `noise`.")
  }
  x
}

# This function reports that signal and noise roles are missing.
#
# The roles are never inferred from the values themselves. Nothing in a pair
# such as `cued` and `uncued` says which one is the signal, and taking the
# wrong one reverses the sign of every d' that follows.
.require_levels <- function(v, arg, label) {

  # The variable must contain two observed values.
  u <- if (is.factor(v)) levels(droplevels(v)) else sort(unique(v[!is.na(v)]))
  if (length(u) != 2L) {
    .usdt_stop(label, " must have exactly 2 distinct values, but it has ",
               length(u), if (length(u) <= 6L)
                 paste0(" (", paste(u, collapse = ", "), ")") else "", ".\n",
               "  If it is a continuous measure, name that task in ",
               "`dichotomize` so it is median-split first.")
  }

  # The two observed values go into the message so they can be copied.
  .usdt_stop("`", arg, "` is required. ", label, " holds `", u[1L], "` and `",
             u[2L], "`, and which of them is the signal cannot be read from ",
             "the values themselves.\n",
             "  Set `", arg, " = c(signal = <one>, noise = <the other>)`.")
}

# Data aggregation

# This function counts responses and trials in each cell.
.agg_counts <- function(response, by) {

  # One key identifies each cell.
  key <- interaction(by, drop = TRUE, sep = "\r")

  # The code counts responses and trials for each key.
  y <- rowsum(as.numeric(response), key, reorder = TRUE)
  n <- rowsum(rep.int(1, length(key)), key, reorder = TRUE)

  # The labels restore the original grouping columns.
  parts <- do.call(rbind, strsplit(rownames(y), "\r", fixed = TRUE))
  out <- as.data.frame(parts, stringsAsFactors = FALSE)
  names(out) <- names(by)

  # Each grouping column returns to its original type.
  for (j in names(by)) out[[j]] <- .restore_type(out[[j]], by[[j]])

  out$y <- as.integer(round(y[, 1L]))
  out$n <- as.integer(round(n[, 1L]))
  rownames(out) <- NULL
  out
}

# This function restores the original column type.
.restore_type <- function(x, template) {
  if (is.numeric(template)) return(as.numeric(x))
  if (is.logical(template)) return(as.logical(x))
  if (is.factor(template))  return(factor(x, levels = levels(template)))
  x
}

# Console formatting

# This function chooses characters that the console can display.
.usdt_chars <- function() {
  if (isTRUE(l10n_info()[["UTF-8"]])) {
    list(h = "\u2500", dot = "\u00b7", delta = "\u0394")
  } else {
    list(h = "-", dot = "*", delta = "Delta ")
  }
}

# This function draws a simple line around a title.
.rule <- function(title = NULL, right = NULL, width = NULL) {
  ch <- .usdt_chars()$h
  width <- min(if (is.null(width)) getOption("width", 80L) else width, 90L)
  left  <- if (is.null(title)) strrep(ch, 2L) else
    paste0(strrep(ch, 2L), " ", title, " ")
  tail  <- if (is.null(right)) "" else paste0(" ", right, " ", strrep(ch, 2L))
  fill  <- max(width - nchar(left) - nchar(tail), 0L)
  paste0(left, strrep(ch, fill), tail)
}

# This function formats p-values for psychology reports.
.fmt_p <- function(p) {
  ifelse(is.na(p), "     NA",
         ifelse(p < .001, "  <.001",
                formatC(sub("^0\\.", ".", formatC(p, format = "f", digits = 3)),
                        width = 7)))
}

# This function prints a number with fixed decimals.
.fmt_n <- function(x, digits = 4L, width = digits + 3L) {
  ifelse(is.na(x), formatC("NA", width = width),
         formatC(x, format = "f", digits = digits, width = width))
}

# This function prints finite and infinite interval limits.
.fmt_ci <- function(lo, hi, digits = 3L) {
  f <- function(v) {
    if (is.na(v)) return("NA")
    if (is.infinite(v)) return(if (v > 0) "Inf" else "-Inf")
    formatC(v, format = "f", digits = digits)
  }
  br <- if ((!is.na(lo) && is.infinite(lo)) ||
            (!is.na(hi) && is.infinite(hi))) c("(", ")") else c("[", "]")
  paste0(br[1L], formatC(f(lo), width = 7L), ", ",
         formatC(f(hi), width = 6L), br[2L])
}

# This function adds separators to trial counts.
.fmt_int <- function(x) formatC(x, format = "d", big.mark = ",")

