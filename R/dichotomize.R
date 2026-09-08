# dichotomize.R
# This script turns a continuous measure into a binary response.
# Author: Ricardo Rey-Sáez
# Last modified: 08-09-2026

# Public functions

#' Turn a continuous measure into a binary response
#'
#' Splits a continuous measure, usually response times, at each subject's own
#' median. The median uses all the trials of that subject, without separating
#' the conditions. The binary result can then be analysed on the same
#' sensitivity scale as a direct task that already gives binary responses. The
#' procedure follows Meyen et al. (2022).
#'
#' @param x Numeric vector with the continuous measure, usually response times.
#' @param by Vector identifying the subject of each value in `x`. Every subject
#'   receives their own median.
#' @param signal Which side of the median counts as a signal response. Use
#'   `"faster"` when the signal condition speeds responses up, as in priming
#'   and cueing tasks. Use `"slower"` when the signal condition slows responses
#'   down, as in interference tasks.
#' @param ties What to do with trials that fall exactly on the median.
#'   `"noise"` gives them the noise response. `"random"` assigns them at
#'   random, which keeps the split as close to even as the data allow.
#'
#' @return An integer vector as long as `x`, with `1` for signal responses and
#'   `0` for noise responses. Missing values in `x` stay missing.
#'
#' @details
#' A subject with an odd number of trials cannot be split into two equal
#' halves, because the median is one of the observed values. With
#' `ties = "noise"` that subject gets a proportion of `(n - 1) / (2n)` signal
#' responses instead of 0.5, so the criterion moves slightly away from zero.
#' The difference is `1 / (2n)` and rarely matters. Setting `ties = "random"`
#' removes it, because rounding up or down at random is unbiased across
#' subjects.
#'
#' @references
#' Meyen, S., Zerweck, I. A., Amado, C., von Luxburg, U., & Franz, V. H.
#' (2022). Advancing research on unconscious priming: When can scientists claim
#' an indirect task advantage? *Journal of Experimental Psychology: General*,
#' 151(1), 65-81. \doi{10.1037/xge0001065}
#'
#' @examples
#' rt   <- c(320, 410, 295, 500, 380, 450)
#' subj <- rep(c("s1", "s2"), each = 3)
#' meyen_split(rt, by = subj)
#'
#' @export
meyen_split <- function(x, by,
                        signal = c("faster", "slower"),
                        ties   = c("noise", "random")) {

  # The function checks the requested options and the input data.
  signal <- match.arg(signal)
  ties   <- match.arg(ties)
  if (!is.numeric(x)) {
    .usdt_stop("`x` must be numeric; `meyen_split()` needs the raw continuous ",
               "measure (e.g. response times), not an already-binary response.")
  }
  if (length(by) != length(x)) {
    .usdt_stop("`by` has length ", length(by), " but `x` has length ",
               length(x), "; they must match.")
  }
  if (all(is.na(x))) .usdt_stop("`x` is entirely missing.")

  # Each subject receives one median across both conditions.
  grp <- as.character(by)
  med <- stats::ave(x, grp, FUN = function(v) stats::median(v, na.rm = TRUE))

  # Values on the selected side become signal responses.
  out <- if (signal == "faster") as.integer(x < med) else as.integer(x > med)

  # The random option redistributes values that equal the median.
  if (ties == "random") out <- .break_ties(x, med, grp, out)

  # Missing input values remain missing.
  out[is.na(x)] <- NA_integer_
  out
}

# Internal functions

# This function assigns ties while keeping each split nearly even.
.break_ties <- function(x, med, grp, out) {

  # The function returns early when no value equals the median.
  tied <- which(!is.na(x) & x == med)
  if (!length(tied)) return(out)

  # Each subject's ties are handled together.
  for (g in unique(grp[tied])) {
    idx  <- which(grp == g & !is.na(x))
    tidx <- intersect(tied, idx)
    n    <- length(idx)

    # The target keeps the split as even as possible.
    target <- if (n %% 2L == 0L) n %/% 2L else n %/% 2L + sample.int(2L, 1L) - 1L

    # A random set of tied trials fills the target.
    out[tidx] <- 0L
    need <- max(min(target - sum(out[idx], na.rm = TRUE), length(tidx)), 0L)
    if (need > 0L) out[tidx[sample.int(length(tidx), need)]] <- 1L
  }
  out
}
