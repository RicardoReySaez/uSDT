# dichotomize.R
# This script turns a continuous measure into a binary response.
# Author: Ricardo Rey-Sáez
# Last modified: 04-09-2026

# Public functions

#' Meyen median-split dichotomization
#'
#' Dichotomizes a continuous measure (typically response times) with a median
#' split computed **within subject, pooling conditions**, following Meyen
#' et al. (2022). This places a continuous indirect measure on the same
#' sensitivity scale as a binary direct measure.
#'
#' @param x Numeric vector with the continuous measure, usually response times.
#' @param by Grouping vector identifying the subject. The median is
#'   computed within each group across all of that subject's trials.
#' @param signal Which side of the median counts as a *signal* response.
#'   `"faster"` (the default) suits priming and cueing tasks, where the signal
#'   condition speeds responses up. Use `"slower"` for tasks in which the
#'   signal condition slows responses down, such as interference paradigms.
#' @param ties How to handle trials falling exactly on the median. `"noise"`
#'   (the default) assigns them the noise response, reproducing
#'   `as.integer(median(rt) > rt)`. `"random"` assigns them at random so that
#'   the split is as close to 50/50 as the number of trials allows.
#'
#' @return An integer vector of the same length as `x`, with `1` for signal
#'   responses and `0` for noise responses. `NA` in `x` propagates.
#'
#' @details
#' With an odd number of trials an exact 50/50 split is impossible: the median
#' is itself an observed value. `ties = "noise"` then yields a proportion of
#' `(n - 1) / (2n)` rather than `0.5`, and the criterion is no longer exactly
#' zero. The deviation is negligible in practice (it is `1 / (2n)`), but
#' `ties = "random"` removes it by rounding up or down at random, which is
#' unbiased across subjects.
#'
#' @references
#' Meyen, S., Zerweck, I. A., Amado, C., von Luxburg, U., & Franz, V. H.
#' (2022). Advancing research on unconscious priming: When can scientists claim
#' an indirect task advantage? *Journal of Experimental Psychology: General*.
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
