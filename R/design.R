# design.R
# This script builds the model formula and records its structure.
# Author: Ricardo Rey-Sáez
# Last modified: 04-09-2026

# Internal functions

# This function builds the formula and counts its parameters.
.usdt_formula <- function(d, fix_criteria = c("auto", "none"),
                          response = c("cbind", "trials")) {

  # The function checks the requested options.
  fix_criteria <- match.arg(fix_criteria)
  response     <- match.arg(response)

  # The criterion settings select the criterion terms.
  czero <- d$meta$criterion_zero
  keep  <- c(direct = TRUE, indirect = TRUE)
  if (fix_criteria == "auto") keep <- !czero
  crit  <- c("c_D", "c_I")[c(keep[["direct"]], keep[["indirect"]])]
  slope <- c("d_D", "d_I")

  # Each model engine receives its expected response format.
  lhs <- if (response == "cbind") "cbind(y, n - y)" else "y | trials(n)"

  # The fixed part has no overall intercept.
  fixed <- paste(c("0", crit, slope), collapse = " + ")

  # Separate blocks keep the two sensitivities correlated.
  blocks <- c(if (length(crit))
                paste0("(0 + ", paste(crit, collapse = " + "), " | subj)"),
              paste0("(0 + ", paste(slope, collapse = " + "), " | subj)"))

  # The code counts the resulting variance terms.
  nv <- length(crit) * (length(crit) + 1L) / 2L + 3L

  # The result includes the formula and its main features.
  list(
    formula   = stats::as.formula(
      paste(lhs, "~", paste(c(fixed, blocks), collapse = " + ")),
      env = parent.frame()),
    criteria  = crit,
    slopes    = slope,
    dropped   = c("c_D", "c_I")[!c(keep[["direct"]], keep[["indirect"]])],
    n_fixed   = length(crit) + 2L,
    n_var     = nv
  )
}

# This function wraps a long formula for printed output.
.formula_lines <- function(f, width = 74L) {
  txt <- paste(deparse(f, width.cutoff = 500L), collapse = " ")
  txt <- gsub("\\s+", " ", txt)
  if (nchar(txt) <= width) return(txt)

  # A long formula breaks before its first random term.
  at <- regexpr(" + (0 +", txt, fixed = TRUE)
  if (at < 0L) return(strwrap(txt, width))
  head <- substr(txt, 1L, at + 1L)
  tail <- substr(txt, at + 3L, nchar(txt))
  pad  <- strrep(" ", regexpr("~", head, fixed = TRUE) + 1L)
  c(head, paste0(pad, tail))
}
