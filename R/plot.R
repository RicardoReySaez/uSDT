# plot.R
# Visualize results from fitted hierarchical SDT models
# Author: Ricardo Rey-Sáez
# Last modified: 02-11-2026

#' Diagnostic and analytical plots for hierarchical SDT models
#'
#' Generates diagnostic visualizations for a fitted hierarchical SDT model,
#' including observed versus latent regression (H3), shrinkage patterns,
#' participant-level caterpillar intervals, and model-implied ROC curves.
#'
#' @param x An `hsdt` object fitted by [hsdt()], frequentist or Bayesian.
#' @param type Character string indicating the plot type:
#'   * `"regression"`: Compares the observed OLS regression with the latent
#'     regression line (H3).
#'   * `"shrinkage"`: Connects observed \eqn{d'} to model-implied \eqn{d'} for
#'     each participant across bivariate contours.
#'   * `"caterpillar"`: Compares observed \eqn{d'} and model-implied \eqn{d'}
#'     with confidence intervals for each participant alongside zero.
#'   * `"roc"`: Shows model-implied ROC curves with estimated criteria.
#' @param subject_id Identifier for a specific participant when `type = "roc"`.
#'   If omitted, displays population-level curves.
#' @param band Logical. Show confidence bands around regression lines or ROC
#'   curves (default is `TRUE`).
#' @param population_reference Logical. When plotting an individual ROC curve,
#'   add population-average curves as dashed reference lines.
#' @param observed_se Variance formulation for empirical intervals in
#'   `"caterpillar"` plots: `"gg"` (Gourevitch & Galanter, 1967, default) or
#'   `"miller"` (Miller, 1996).
#' @param nrow Number of rows of the two panels when `type = "regression"`:
#'   `1` (default) places them side by side as squares, `2` stacks the
#'   observed panel above the model one, each as wide as the figure.
#' @param caption Logical. Show the caption that explains the plot under it
#'   (default `TRUE`). The shrinkage plot has none.
#' @param ... Additional arguments passed to underlying plotting methods.
#'
#' @return A `ggplot` object. Its underlying data frame is stored in `$data`.
#'
#' @details
#' # Regression plot (`type = "regression"`)
#'
#' Compares observed and latent associations across two panels sharing axes,
#' side by side or, with `nrow = 2`, one above the other.
#' The first panel shows observed \eqn{d'} values and an ordinary least-squares
#' line. When awareness task reliability is low, trial-level sampling noise
#' attenuates this observed slope toward zero. The second panel plots the latent
#' regression line (\eqn{d'_I} on \eqn{d'_A}) from H3, correcting for
#' measurement error. The value of each line at \eqn{d'_A = 0} marks the
#' intercept testing for unconscious processing, and both panels display it the
#' same way: an open circle at the point estimate with a vertical line spanning
#' its confidence interval. The observed marker is the least-squares intercept
#' and the latent marker is its measurement-error-corrected counterpart, so the
#' two panels place the same hypothesis side by side. Confidence bands are
#' computed via the delta method or bootstrap replicates when [usdt_boot()] is
#' present.
#'
#' # Shrinkage plot (`type = "shrinkage"`)
#'
#' Connects each participant's observed \eqn{d'} (from [sdt_moments()] with Hautus
#' correction) to their model-implied \eqn{d'}. The lower the reliability of the
#' measures, the higher the shrinkage of observed estimates toward the
#' group-level mean.
#'
#' # Caterpillar plot (`type = "caterpillar"`)
#'
#' Plots observed \eqn{d'} and model-implied \eqn{d'} with confidence intervals
#' for every participant. Empirical intervals use normal approximations based
#' on `observed_se`. Model-implied intervals incorporate uncertainty from
#' population means, variance components, and participant random effects.
#'
#' # ROC plot (`type = "roc"`)
#'
#' Displays model-implied ROC curves for an average participant or a specific
#' individual, with points marking the estimated response criteria.
#'
#' # Bayesian fits
#'
#' For a fit with `estimation = "bayesian"`, every model-estimated value is a
#' posterior mean and every model band is a central credible interval,
#' computed draw by draw: the latent line and its band, each subject's
#' sensitivity in the shrinkage and caterpillar plots, and the ROC curves and
#' operating points. The latent p-values are the two-sided posterior p-values
#' of the hypothesis table. Under unequal variances the ROC curves bend with
#' the signal standard deviation of each draw.
#'
#' @seealso [hsdt()], [sdt_moments()], [usdt_boot()]
#'
#' @examples
#' # Contextual cuing data from Vadillo et al. (2025)
#' d <- usdt_data_tasks(
#'   awareness = vadillo_awareness,
#'   influence = vadillo_cuing,
#'   subject_col      = "subj",
#'   condition_col    = "condition",
#'   condition_levels = c(signal = "old", noise = "new"),
#'   response_col     = list(awareness = "judged.old", influence = "rt"),
#'   response_levels  = list(awareness = c(signal = 1, noise = 0),
#'                           influence = c(signal = "faster",
#'                                         noise = "slower")),
#'   dichotomize      = list(awareness = FALSE, influence = TRUE)
#' )
#'
#' fit <- hsdt(d)
#'
#' # 1. Observed vs. latent regression (H3), side by side or stacked
#' plot(fit, type = "regression")
#' plot(fit, type = "regression", nrow = 2)
#'
#' # 2. Bivariate shrinkage toward the group mean
#' plot(fit, type = "shrinkage")
#'
#' # 3. Participant-level intervals (observed vs. model-implied)
#' plot(fit, type = "caterpillar")
#'
#' # 4. Model-implied ROC curves
#' plot(fit, type = "roc")
#' plot(fit, type = "roc", subject_id = "2001", population_reference = TRUE)
#'
#' @export
plot.hsdt <- function(x, type = c("regression", "shrinkage",
                                  "caterpillar", "roc"),
                      subject_id = NULL, band = TRUE,
                      population_reference = TRUE, observed_se = NULL,
                      nrow = NULL, caption = TRUE, ...) {

  # The function checks the requested plot.
  type <- match.arg(type)
  .check_caption(caption)
  if (!is.logical(band) || length(band) != 1L || is.na(band)) {
    .usdt_stop("`band` must be `TRUE` or `FALSE`.")
  }
  if (!is.logical(population_reference) ||
      length(population_reference) != 1L || is.na(population_reference)) {
    .usdt_stop("`population_reference` must be `TRUE` or `FALSE`.")
  }
  if (length(list(...))) {
    .usdt_stop("`...` is not used by the selected plot.")
  }
  if (type != "roc" && !is.null(subject_id)) {
    .usdt_stop("`subject_id` is only used by `type = \"roc\"`.")
  }
  if (type != "caterpillar" && !is.null(observed_se)) {
    .usdt_stop("`observed_se` is only used by `type = \"caterpillar\"`.")
  }
  observed_se <- if (is.null(observed_se)) "gg" else
    match.arg(observed_se, c("gg", "miller"))
  if (type != "regression" && !is.null(nrow)) {
    .usdt_stop("`nrow` is only used by `type = \"regression\"`.")
  }
  if (is.null(nrow)) nrow <- 1L
  if (!is.numeric(nrow) || length(nrow) != 1L || !nrow %in% 1:2) {
    .usdt_stop("`nrow` must be 1 (panels side by side) or 2 (stacked).")
  }

  switch(type,
         regression = .plot_regression(x, band, nrow),
         shrinkage = .plot_shrinkage(x),
         caterpillar = .plot_caterpillar(x, observed_se),
         roc = .plot_roc(x, subject_id, band, population_reference)) +
    .hide_caption(caption)
}

# This function checks the `caption` argument of the plot methods.
.check_caption <- function(caption) {
  if (!is.logical(caption) || length(caption) != 1L || is.na(caption)) {
    .usdt_stop("`caption` must be `TRUE` or `FALSE`.")
  }
}

# This function removes the caption of a plot when `caption` is FALSE; added
# to a plot, its NULL otherwise changes nothing.
.hide_caption <- function(caption) {
  if (!caption) ggplot2::labs(caption = NULL)
}

# This function gives the caption and the plot margins of a theme. ggplot2
# never wraps a caption, and one wrapped beforehand at a fixed number of
# characters fills only part of a wide figure, so the caption element wraps
# when it is drawn, to the width of the plot. Its left and right margins,
# those of the plot, are kept from the edges of the figure. `margin` holds the
# top, right, bottom and left plot margins in pt.
.caption_theme <- function(margin) {
  caption <- structure(
    list(family = NULL, face = NULL, colour = "#687076", size = 9.5,
         hjust = 0, vjust = 1, angle = NULL, lineheight = NULL,
         margin = ggplot2::margin(t = 10, r = margin[2L], l = margin[4L]),
         debug = NULL, inherit.blank = FALSE),
    class = c("usdt_caption", "element_text", "element")
  )
  ggplot2::theme(plot.caption = caption, plot.caption.position = "plot",
                 plot.margin = do.call(ggplot2::margin, as.list(margin)))
}

# The caption element draws its text as a grob that keeps its measured words
# and line breaks in an environment, which the copy of the grob that sizes the
# caption row shares with the copy that is drawn. Paragraphs are separated by
# "\n".
#' @exportS3Method ggplot2::element_grob
element_grob.usdt_caption <- function(element, label = "", ...) {
  if (is.null(label) || !any(nzchar(trimws(label)))) return(ggplot2::zeroGrob())
  gp <- list(fontfamily = element$family, fontface = element$face,
             fontsize = element$size, col = element$colour,
             lineheight = element$lineheight)
  grid::gTree(text = paste(label, collapse = "\n"), hjust = element$hjust,
              margin = element$margin, state = new.env(parent = emptyenv()),
              gp = do.call(grid::gpar, gp[!vapply(gp, is.null, TRUE)]),
              cl = "usdt_caption")
}

# grid sizes the caption row in the viewport that holds the whole plot, before
# it knows how wide the plot will be: a square panel makes the plot narrower
# than the figure, and how much narrower depends on the height of the caption
# itself. So the caption reserves the lines it needs across the figure, which
# is never too few.
#' @exportS3Method grid::heightDetails
heightDetails.usdt_caption <- function(x) {
  m <- grid::convertWidth(x$margin, "in", valueOnly = TRUE)
  s <- x$state
  s$words <- .caption_words(x$text)
  s$left <- .device_x(0) + m[4L]
  s$width <- .device_x(1) - m[2L] - s$left
  s$lines <- .caption_lines(s$words, s$width)
  height <- grid::grobHeight(grid::textGrob(paste(s$lines, collapse = "\n")))
  grid::unit(grid::convertHeight(height, "in", valueOnly = TRUE) + m[1L] + m[3L],
             "in")
}

# The caption is drawn as wide as the plot when its text fits there in the
# lines it reserved, which it always does when the plot fills the figure;
# otherwise it takes the narrowest width at which it does fit, centred under
# the plot.
#' @exportS3Method grid::makeContent
makeContent.usdt_caption <- function(x) {
  s <- x$state
  if (is.null(s$lines)) heightDetails.usdt_caption(x)
  plot <- c(.device_x(0), .device_x(1))
  fits <- function(w) length(.caption_lines(s$words, w)) <= length(s$lines)
  low <- min(diff(plot), s$width)
  width <- low
  if (!fits(low)) {
    high <- s$width
    while (high - low > 0.01) {
      mid <- (low + high) / 2
      if (fits(mid)) high <- mid else low <- mid
    }
    width <- high
  }
  left <- if (width <= diff(plot)) plot[1L] else mean(plot) - width / 2
  left <- min(max(left, s$left), s$left + s$width - width)
  top <- grid::convertHeight(x$margin[1L], "in", valueOnly = TRUE)
  grid::setChildren(x, grid::gList(grid::textGrob(
    paste(.caption_lines(s$words, width), collapse = "\n"),
    x = grid::unit(left - plot[1L] + x$hjust * width, "in"),
    y = grid::unit(1, "npc") - grid::unit(top, "in"),
    hjust = x$hjust, vjust = 1)))
}

# This function gives the position on the device, in inches, of a point at a
# fraction `npc` of the width of the current viewport.
.device_x <- function(npc) {
  grid::deviceLoc(grid::unit(npc, "npc"), grid::unit(0, "npc"),
                  valueOnly = TRUE)$x
}

# This function measures, in the current font, the words of each paragraph of
# a caption and the space between them, in inches.
.caption_words <- function(text) {
  paragraphs <- trimws(strsplit(text, "\n", fixed = TRUE)[[1L]])
  inches <- function(s) grid::convertWidth(grid::stringWidth(s), "in",
                                           valueOnly = TRUE)
  words <- lapply(paragraphs[nzchar(paragraphs)], function(p) {
    words <- strsplit(p, " +")[[1L]]
    list(words = words, sizes = inches(words))
  })
  list(paragraphs = words, space = inches(" "))
}

# This function breaks the measured paragraphs into lines no wider than
# `width` inches.
.caption_lines <- function(words, width) {
  space <- words$space
  unlist(lapply(words$paragraphs, function(p) {
    end <- cumsum(p$sizes + space) - space
    rest <- p$words
    out <- character()
    while (length(rest)) {
      n <- max(1L, sum(end <= width))
      out <- c(out, paste(rest[seq_len(n)], collapse = " "))
      rest <- rest[-seq_len(n)]
      end <- end[-seq_len(n)] - end[n] - space
    }
    out
  }))
}

# This function evaluates the latent line and its band over a grid of x.
.latent_line <- function(object, x, band = TRUE) {

  # The band follows whichever method the hypothesis table used, and names the
  # reason whenever it cannot be drawn.
  limits <- matrix(NA_real_, length(x), 2L)
  method <- "none"
  reason <- if (band) NA_character_ else "it was turned off"

  if (.is_bayes(object)) {
    # A posterior gives the line its mean and its band draw by draw.
    Q <- .draw_quantities(.usdt_bivariate(.population_draws(object)))
    lines <- Q[, c("intercept", "slope")] %*% rbind(1, x)
    fit <- colMeans(lines)
    slope <- mean(Q[, "slope"])
    intercept <- mean(Q[, "intercept"])
    if (band) {
      limits <- .boot_ci(lines, fit, object$level, "perc")
      method <- "posterior"
    }
  } else {
    e <- object$pars$est
    q <- .usdt_quantities(e)
    slope <- q[, "slope"]
    intercept <- q[, "intercept"]
    fit <- e[["gamma_I"]] + slope * (x - e[["gamma_A"]])
    if (!band) {
      # The reason above already says why there is no band.
    } else if (.has_boot_summary(object)) {
      replicates <- object$boot$t[object$boot$ok, c("intercept", "slope"),
                                  drop = FALSE]
      limits <- .boot_ci(replicates %*% rbind(1, x), fit, object$level,
                         object$boot$type)
      method <- object$boot$type
    } else if (isTRUE(object$pars$joint_ok)) {
      critical <- stats::qnorm(1 - (1 - object$level) / 2)
      se <- .delta_se(.reg_grad(e, x), object$pars$vcov)
      limits <- cbind(fit - critical * se, fit + critical * se)
      method <- "delta"
    } else {
      reason <- object$pars$inference_reason
    }
  }

  list(line = data.frame(x = x, fit = fit, conf.low = limits[, 1L],
                         conf.high = limits[, 2L], stringsAsFactors = FALSE),
       slope = slope, intercept = intercept, method = method, reason = reason)
}

# This function formats one p-value as a plotmath fragment.
.plot_p <- function(p) {
  if (!is.finite(p)) return('italic(p) == "NA"')
  if (p < .001) 'italic(p) < ".001"' else
    paste0('italic(p) == "', sub("^0", "", sprintf("%.3f", p)), '"')
}

# This function writes only the intercept and slope as a two-line expression.
.regression_label <- function(intercept, slope, p_intercept, p_slope) {
  value <- function(x) if (is.finite(x)) sprintf("%.3f", x) else '"NA"'
  paste0("atop(italic(b)[0] == ", value(intercept),
         "~(", .plot_p(p_intercept), "),",
         "italic(b)[1] == ", value(slope),
         "~(", .plot_p(p_slope), "))")
}

# This function collects the observed and model-estimated clouds and their lines.
.regression_data <- function(object, band = TRUE) {

  points <- .shrinkage_data(object)
  panels <- c("Observed estimates", "Model-estimated sensitivities")
  panel_factor <- function(x) factor(x, levels = panels)

  # Least squares supplies the observed line, its confidence band and its tests.
  spread <- stats::var(points$observed_awareness)
  naive_ok <- nrow(points) >= 2L && isTRUE(spread > 0)
  naive_fit <- if (naive_ok) {
    stats::lm(observed_influence ~ observed_awareness, data = points)
  } else NULL
  # Zero joins the observed grid for the same reason it joins the latent one:
  # the value of the line there is the quantity both panels put on trial.
  naive_span <- if (naive_ok) range(c(0, points$observed_awareness)) else
    numeric(0)
  naive_x <- if (naive_ok) {
    sort(unique(c(seq(naive_span[1L], naive_span[2L], length.out = 127L), 0)))
  } else numeric(0)
  naive_line <- data.frame(
    x = naive_x, fit = rep(NA_real_, length(naive_x)),
    conf.low = rep(NA_real_, length(naive_x)),
    conf.high = rep(NA_real_, length(naive_x)),
    panel = rep(panels[1L], length(naive_x))
  )
  naive_line$panel <- panel_factor(naive_line$panel)
  naive_coef <- c(intercept = NA_real_, slope = NA_real_)
  naive_p <- c(intercept = NA_real_, slope = NA_real_)
  if (naive_ok) {
    naive_coef[] <- stats::coef(naive_fit)
    coef_table <- summary(naive_fit)$coefficients
    naive_p[] <- coef_table[, "Pr(>|t|)"]
    prediction <- suppressWarnings(stats::predict(
      naive_fit, newdata = data.frame(observed_awareness = naive_x),
      interval = if (band) "confidence" else "none", level = object$level
    ))
    if (band && is.matrix(prediction)) {
      naive_line[c("fit", "conf.low", "conf.high")] <-
        prediction[, c("fit", "lwr", "upr"), drop = FALSE]
    } else {
      naive_line$fit <- as.numeric(prediction)
    }
  }

  # Read at zero, the observed band returns the least-squares intercept and its
  # interval, so marker and ribbon agree exactly as they do in the latent panel.
  naive_origin <- if (naive_ok) {
    naive_line[match(0, naive_x), , drop = FALSE]
  } else naive_line[0L, , drop = FALSE]

  # The latent line is only drawn over the central fitted population range,
  # plus zero because its value there is the hypothesis of interest. A
  # posterior places that range with its mean values.
  e <- if (.is_bayes(object)) {
    colMeans(.draw_primitives(.usdt_bivariate(.population_draws(object))))
  } else {
    object$pars$est
  }
  latent_span <- if (is.finite(e[["s2_A"]]) && e[["s2_A"]] > 0) {
    e[["gamma_A"]] + stats::qnorm(c(.025, .975)) * sqrt(e[["s2_A"]])
  } else {
    range(points$model_awareness)
  }
  latent_span <- range(c(0, latent_span))
  latent_x <- sort(unique(c(seq(latent_span[1L], latent_span[2L],
                                length.out = 127L), 0)))
  latent <- .latent_line(object, latent_x, band)
  latent$line$panel <- panel_factor(rep(panels[2L], nrow(latent$line)))
  origin <- latent$line[match(0, latent_x), , drop = FALSE]

  # The H3 table is the single source for latent p-values, including bootstrap.
  latent_tests <- object$tests[object$tests$hypothesis == "H3", ]
  latent_p <- stats::setNames(latent_tests$p.value, latent_tests$term)
  annotations <- data.frame(
    panel = panel_factor(panels), x = Inf, y = Inf,
    label = c(
      .regression_label(naive_coef[["intercept"]], naive_coef[["slope"]],
                        naive_p[["intercept"]], naive_p[["slope"]]),
      .regression_label(latent$intercept, latent$slope,
                        latent_p[["intercept"]], latent_p[["slope"]])
    ),
    stringsAsFactors = FALSE
  )

  # Repeating the subjects in both facets preserves a useful data component.
  observed_panel <- transform(points, panel = panels[1L])
  model_panel <- transform(points, panel = panels[2L])
  observed_panel$panel <- panel_factor(observed_panel$panel)
  model_panel$panel <- panel_factor(model_panel$panel)
  panel_points <- rbind(observed_panel, model_panel)
  panel_points$panel <- panel_factor(panel_points$panel)

  list(points = points, panel_points = panel_points,
       observed_points = observed_panel, model_points = model_panel,
       naive_line = naive_line, naive_intercept = naive_coef[["intercept"]],
       naive_slope = naive_coef[["slope"]], naive_p = naive_p,
       naive_origin = naive_origin,
       line = latent$line, slope = latent$slope, intercept = latent$intercept,
       method = latent$method, reason = latent$reason, origin = origin,
       annotations = annotations, panels = panels)
}

# This function separates the observed regression from the model representation.
.plot_regression <- function(object, band, nrow = 1L) {

  values <- .regression_data(object, band)
  labels <- object$data$meta$labels
  origin <- values$origin
  naive_origin <- values$naive_origin
  banded <- values$method != "none"
  naive <- is.finite(values$naive_slope)
  naive_banded <- band && naive && nrow(naive_origin) == 1L &&
    all(is.finite(c(naive_origin$conf.low, naive_origin$conf.high)))

  background <- "#FFFFFF"
  observed <- .observed_style()
  observed_colour <- observed$colour
  latent_colour <- "#116B60"
  naive_colour <- "#B6543A"
  # Common axes make the raw and partially pooled subject positions comparable.
  plot <- ggplot2::ggplot(values$panel_points) +
    ggplot2::geom_hline(
      yintercept = 0, colour = "#555B60",
      linetype = "dotted", linewidth = 0.5
    ) +
    ggplot2::geom_vline(
      xintercept = 0, colour = "#555B60",
      linetype = "dotted", linewidth = 0.5
    )

  # The model panel shows where partial pooling moves every observed subject.
  # The segments and the dotted lines pass under the points, whose discs are
  # cleared first; the bands come after, so they tint the points as they tint
  # the background.
  plot <- plot +
    ggplot2::geom_segment(
      data = values$model_points,
      ggplot2::aes(x = .data[["observed_awareness"]],
                   y = .data[["observed_influence"]],
                   xend = .data[["model_awareness"]],
                   yend = .data[["model_influence"]]),
      inherit.aes = FALSE, colour = observed_colour, alpha = 0.34,
      linewidth = 0.42
    ) +
    .knockout(values$panel_points, "observed_awareness", "observed_influence",
              size = 2.4, stroke = 0.7) +
    .knockout(values$model_points, "model_awareness", "model_influence",
              size = 2.65, stroke = 0.75)

  # Each panel receives the interval appropriate to the line it contains.
  if (band && naive && all(is.finite(values$naive_line$conf.low))) {
    plot <- plot + ggplot2::geom_ribbon(
      data = values$naive_line,
      ggplot2::aes(x = .data[["x"]], ymin = .data[["conf.low"]],
                   ymax = .data[["conf.high"]]),
      inherit.aes = FALSE, fill = naive_colour, alpha = 0.11
    )
  }
  if (banded) {
    plot <- plot + ggplot2::geom_ribbon(
      data = values$line,
      ggplot2::aes(x = .data[["x"]], ymin = .data[["conf.low"]],
                   ymax = .data[["conf.high"]]),
      inherit.aes = FALSE, fill = latent_colour, alpha = 0.13
    )
  }

  plot <- plot +
    ggplot2::geom_point(
      data = values$observed_points,
      ggplot2::aes(x = .data[["observed_awareness"]],
                   y = .data[["observed_influence"]]),
      inherit.aes = FALSE, shape = 21, size = 2.4, stroke = 0.7,
      colour = observed_colour, fill = observed$fill, alpha = observed$alpha
    ) +
    ggplot2::geom_point(
      data = values$model_points,
      ggplot2::aes(x = .data[["observed_awareness"]],
                   y = .data[["observed_influence"]]),
      inherit.aes = FALSE, shape = 21, size = 2.4, stroke = 0.7,
      colour = observed_colour, fill = observed$fill,
      alpha = observed$alpha / 2
    ) +
    ggplot2::geom_point(
      data = values$model_points,
      ggplot2::aes(x = .data[["model_awareness"]],
                   y = .data[["model_influence"]]),
      inherit.aes = FALSE,
      shape = 21, size = 2.65, stroke = 0.75,
      colour = latent_colour, fill = "#3FA88E", alpha = 0.9
    )
  if (naive) {
    plot <- plot + ggplot2::geom_line(
      data = values$naive_line,
      ggplot2::aes(x = .data[["x"]], y = .data[["fit"]]),
      inherit.aes = FALSE, colour = naive_colour,
      linewidth = 0.9, linetype = "dashed"
    )
  }
  plot <- plot + ggplot2::geom_line(
    data = values$line,
    ggplot2::aes(x = .data[["x"]], y = .data[["fit"]]),
    inherit.aes = FALSE, colour = latent_colour, linewidth = 1.05
  )

  # The intercept is the quantity H3 tests, so it gets its own marker in both
  # panels: the attenuated observed one and the corrected latent one answer the
  # same question, and drawing them alike makes the contrast the reader's to
  # make. Their limits are each band read at zero, which keeps every marker
  # identical to the ribbon it sits in.
  if (naive_banded) {
    plot <- plot + ggplot2::geom_linerange(
      data = naive_origin,
      ggplot2::aes(x = .data[["x"]], ymin = .data[["conf.low"]],
                   ymax = .data[["conf.high"]]),
      inherit.aes = FALSE, colour = naive_colour, linewidth = 1
    )
  }
  if (banded && all(is.finite(c(origin$conf.low, origin$conf.high)))) {
    plot <- plot + ggplot2::geom_linerange(
      data = origin,
      ggplot2::aes(x = .data[["x"]], ymin = .data[["conf.low"]],
                   ymax = .data[["conf.high"]]),
      inherit.aes = FALSE, colour = latent_colour, linewidth = 1
    )
  }
  if (naive && nrow(naive_origin) == 1L && is.finite(naive_origin$fit)) {
    plot <- plot + ggplot2::geom_point(
      data = naive_origin,
      ggplot2::aes(x = .data[["x"]], y = .data[["fit"]]),
      inherit.aes = FALSE, shape = 21, size = 3.2, stroke = 1.05,
      colour = naive_colour, fill = background
    )
  }

  # The caption explains interpretation rather than repeating coefficients.
  band_note <- if (!band) {
    "Bands turned off."
  } else switch(values$method,
    none = paste0("The latent band is unavailable: ", values$reason, "."),
    delta = sprintf("Bands are %.0f%% OLS and delta-method confidence intervals.",
                    100 * object$level),
    posterior = sprintf(paste0("Bands are a %.0f%% OLS confidence interval and ",
                               "a %.0f%% credible interval."),
                        100 * object$level, 100 * object$level),
    sprintf("Bands are %.0f%% OLS and bootstrap %s confidence intervals.",
            100 * object$level,
            switch(values$method, perc = "percentile", norm = "normal",
                   basic = "basic")))
  bayes <- .is_bayes(object)
  caption <- paste0(
    "Segments show shrinkage from observed to ",
    if (bayes) "posterior mean" else "conditional model", " estimates. ",
    "The latent line is implied by the fitted random-effects distribution, ",
    "not fitted to the green points. ",
    "An open circle marks the intercept of each panel at zero. ", band_note,
    " P-values test the intercept and slope against zero",
    if (bayes) " (two-sided posterior p-values in the latent panel)" else "",
    "; p < .05 indicates significance."
  )

  plot +
    ggplot2::geom_point(
      data = origin,
      ggplot2::aes(x = .data[["x"]], y = .data[["fit"]]),
      inherit.aes = FALSE, shape = 21, size = 3.2, stroke = 1.05,
      colour = latent_colour, fill = background
    ) +
    ggplot2::geom_text(
      data = values$annotations,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                   label = .data[["label"]]),
      inherit.aes = FALSE, hjust = 1.08, vjust = 1.2,
      lineheight = 1.05, size = 3.35, colour = "#3E4347", parse = TRUE
    ) +
    ggplot2::facet_wrap(~panel, nrow = nrow) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.06)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.06)) +
    ggplot2::labs(
      x = bquote(.(labels[["awareness"]]) ~ sensitivity ~ (italic(d) * minute)),
      y = bquote(.(labels[["influence"]]) ~ sensitivity ~ (italic(d) * minute)),
      caption = caption
    ) +
    .usdt_theme() +
    # Side by side, the panels are squares; stacked, they take the width of
    # the figure, so a wide figure stretches the lines.
    ggplot2::theme(aspect.ratio = if (nrow == 1L) 1,
                   panel.spacing = grid::unit(24, "pt"))
}

# This function pairs the observed sensitivities of every complete subject.
.observed_pairs <- function(object) {

  # The observed values come from the same prepared data as the model.
  observed <- sdt_moments(object$data, correction = "hautus")
  labels <- object$data$meta$labels
  awareness <- observed[observed$task == labels[["awareness"]],
                        c("subj", "dprime")]
  influence <- observed[observed$task == labels[["influence"]],
                        c("subj", "dprime")]
  awareness$subj <- as.character(awareness$subj)
  influence$subj <- as.character(influence$subj)

  # Both observed values place a subject on the plot.
  influence_row <- match(awareness$subj, influence$subj)
  complete <- !is.na(influence_row)
  omitted <- length(unique(c(awareness$subj, influence$subj))) - sum(complete)
  if (omitted > 0L) {
    .usdt_warn(omitted, " subject", if (omitted == 1L) " was" else "s were",
               " omitted because both observed sensitivities are required.")
  }

  data.frame(
    subject = awareness$subj[complete],
    observed_awareness = awareness$dprime[complete],
    observed_influence = influence$dprime[influence_row[complete]],
    stringsAsFactors = FALSE
  )
}

# This function adds the conditional model value of each subject.
.shrinkage_data <- function(object) {

  values <- .observed_pairs(object)

  # Each model value combines the task average with the subject difference;
  # a posterior gives its mean.
  model <- if (.is_bayes(object)) {
    cbind(colMeans(.subject_draws(object, "d", 1L, values$subject)),
          colMeans(.subject_draws(object, "d", 2L, values$subject)))
  } else {
    .task_effects(object$fit, values$subject)$subjects[, c("d_A", "d_I")]
  }
  values$model_awareness <- unname(model[, 1L])
  values$model_influence <- unname(model[, 2L])
  values
}

# This function gives the look of the observed subjects in the regression and
# shrinkage plots: grey circles with a light, translucent fill, so that
# overlapping subjects darken where they gather.
.observed_style <- function() {
  list(colour = "#8F969B", fill = "#C4CACE", alpha = 0.55)
}

# This function draws white discs where points will go, so that the lines
# drawn before them stop at their edge while translucent points still blend
# with one another.
.knockout <- function(data, x, y, size, stroke) {
  ggplot2::geom_point(
    data = data, ggplot2::aes(x = .data[[x]], y = .data[[y]]),
    inherit.aes = FALSE, shape = 21, size = size, stroke = stroke,
    colour = "white", fill = "white"
  )
}

# This function draws the shrinkage comparison.
.plot_shrinkage <- function(object) {

  values <- .shrinkage_data(object)
  labels <- object$data$meta$labels

  # The contours show the shape of both groups of points.
  observed_ellipse <- .ellipse_points(
    values$observed_awareness, values$observed_influence,
    levels = c(0.95, 0.80, 0.50)
  )
  model_ellipse <- .ellipse_points(
    values$model_awareness, values$model_influence,
    levels = 0.90
  )

  # Warm grey and muted teal separate raw values from model estimates.
  observed <- .observed_style()
  observed_colour <- observed$colour
  model_colour <- "#116B60"
  model_fill <- "#3FA88E"

  # Every line passes under the points, whose discs are cleared first; the
  # ellipse fills come after, so they tint the points as they tint the
  # background.
  ellipse <- function(data, ...) {
    if (nrow(data)) ggplot2::geom_polygon(
      data = data,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                   group = .data[["coverage"]]),
      inherit.aes = FALSE, ...
    )
  }
  ggplot2::ggplot(values) +
    ggplot2::geom_hline(
      yintercept = 0, colour = "#555B60",
      linetype = "dotted", linewidth = 0.5
    ) +
    ggplot2::geom_vline(
      xintercept = 0, colour = "#555B60",
      linetype = "dotted", linewidth = 0.5
    ) +
    ellipse(model_ellipse, fill = NA, colour = model_colour, alpha = 0.17,
            linewidth = 0.4) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data[["observed_awareness"]],
                   y = .data[["observed_influence"]],
                   xend = .data[["model_awareness"]],
                   yend = .data[["model_influence"]]),
      colour = observed_colour, alpha = 0.48, linewidth = 0.45
    ) +
    .knockout(values, "observed_awareness", "observed_influence",
              size = 2.55, stroke = 0.75) +
    .knockout(values, "model_awareness", "model_influence",
              size = 2.8, stroke = 0.75) +
    ellipse(observed_ellipse, fill = observed_colour, alpha = 0.06,
            colour = NA) +
    ellipse(model_ellipse, fill = model_fill, alpha = 0.17, colour = NA) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data[["observed_awareness"]],
                   y = .data[["observed_influence"]]),
      shape = 21, size = 2.55, stroke = 0.75,
      colour = observed_colour, fill = observed$fill, alpha = observed$alpha
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data[["model_awareness"]],
                   y = .data[["model_influence"]]),
      shape = 21, size = 2.8, stroke = 0.75,
      colour = model_colour, fill = model_fill, alpha = 0.95
    ) +
    ggplot2::geom_blank(
      data = data.frame(x = 0, y = 0),
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
      inherit.aes = FALSE
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.08)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.08)) +
    ggplot2::labs(
      x = bquote(.(labels[["awareness"]]) ~ sensitivity ~ (italic(d) * minute)),
      y = bquote(.(labels[["influence"]]) ~ sensitivity ~ (italic(d) * minute))
    ) +
    .usdt_theme() +
    ggplot2::theme(aspect.ratio = 1, legend.position = "none")
}

# This function collects both intervals for each subject and task.
.caterpillar_data <- function(object, observed_se = "gg") {

  # The selected variance gives the trial error around each observed value.
  se_col <- paste0("se_", observed_se)
  observed <- sdt_moments(object$data, correction = "hautus", variances = TRUE)
  labels <- object$data$meta$labels
  critical <- stats::qnorm(1 - (1 - object$level) / 2)

  # The model supplies each subject estimate and its interval.
  model <- .subject_intervals(object, critical)

  # Each task keeps the order given by its observed values.
  tasks <- list(
    awareness = list(label = labels[["awareness"]], term = "d_A"),
    influence = list(label = labels[["influence"]], term = "d_I")
  )
  result <- lapply(tasks, function(task) {
    raw <- observed[observed$task == task$label, ]
    raw$subj <- as.character(raw$subj)
    rows <- model[[task$term]]
    model_row <- match(raw$subj, rows$subject)
    if (anyNA(model_row) || any(!is.finite(rows$se[model_row]))) {
      .usdt_stop("the fitted model does not provide intervals for every subject.")
    }

    order_row <- order(raw$dprime, raw$subj)
    rank <- integer(nrow(raw))
    rank[order_row] <- seq_len(nrow(raw))
    observed_rows <- data.frame(
      subject = raw$subj,
      task = task$label,
      method = "Observed",
      estimate = raw$dprime,
      se = raw[[se_col]],
      rank = rank,
      stringsAsFactors = FALSE
    )
    observed_rows$conf.low <- observed_rows$estimate - critical * observed_rows$se
    observed_rows$conf.high <- observed_rows$estimate + critical * observed_rows$se
    model_rows <- data.frame(
      subject = raw$subj,
      task = task$label,
      method = "Model-estimated",
      estimate = rows$estimate[model_row],
      se = rows$se[model_row],
      rank = rank,
      conf.low = rows$conf.low[model_row],
      conf.high = rows$conf.high[model_row],
      stringsAsFactors = FALSE
    )
    rbind(observed_rows, model_rows)
  })

  # The interval colour shows whether zero remains plausible.
  values <- do.call(rbind, result)
  includes_zero <- values$conf.low <= 0 & values$conf.high >= 0
  values$status <- factor(
    ifelse(includes_zero, "Includes zero", "Excludes zero"),
    levels = c("Includes zero", "Excludes zero")
  )
  values$task <- factor(values$task,
                        levels = c(labels[["awareness"]],
                                   labels[["influence"]]))
  values$method <- factor(values$method,
                          levels = c("Observed", "Model-estimated"))
  # Panel titles are drawn on a graphics device, so they stay ASCII: a device
  # in a single-byte locale cannot convert every character, and a failed
  # conversion is an error at render time rather than a missing glyph.
  panel_levels <- c(
    paste0(labels[["awareness"]], " task: Observed"),
    paste0(labels[["awareness"]], " task: Model-estimated"),
    paste0(labels[["influence"]], " task: Observed"),
    paste0(labels[["influence"]], " task: Model-estimated")
  )
  values$panel <- factor(
    paste0(values$task, " task: ", values$method),
    levels = panel_levels
  )
  rownames(values) <- NULL
  values
}

# This function gives each subject's model-estimated sensitivity, with its
# uncertainty and interval, one table per task: the conditional estimate with
# its joint standard error and Wald interval, or the posterior mean, SD and
# central credible interval.
.subject_intervals <- function(object, critical) {
  if (.is_bayes(object)) {
    a <- (1 - object$level) / 2
    out <- lapply(1:2, function(j) {
      draws <- .subject_draws(object, "d", j)
      limits <- apply(draws, 2L, stats::quantile, probs = c(a, 1 - a),
                      names = FALSE)
      data.frame(subject = object$subjects, estimate = colMeans(draws),
                 se = apply(draws, 2L, stats::sd), conf.low = limits[1L, ],
                 conf.high = limits[2L, ], stringsAsFactors = FALSE)
    })
    return(stats::setNames(out, c("d_A", "d_I")))
  }

  # The whole interval needs reliable joint uncertainty.
  if (!isTRUE(object$pars$joint_ok)) {
    .usdt_stop("subject intervals are unavailable because ",
               object$pars$inference_reason, ".")
  }
  conditional <- .conditional_se(object$fit, object$pars$V_full, object$devfun)
  fixed <- lme4::fixef(object$fit)
  out <- lapply(c("d_A", "d_I"), function(term) {
    rows <- conditional[conditional$grpvar == "subj" &
                          conditional$term == term, ]
    estimate <- fixed[[term]] + rows$condval
    data.frame(subject = as.character(rows$grp), estimate = estimate,
               se = rows$se, conf.low = estimate - critical * rows$se,
               conf.high = estimate + critical * rows$se,
               stringsAsFactors = FALSE)
  })
  stats::setNames(out, c("d_A", "d_I"))
}

# This function describes how many intervals include zero in each panel.
.caterpillar_zero_summary <- function(values) {
  panels <- levels(values$panel)
  result <- lapply(panels, function(panel) {
    rows <- values$panel == panel
    percentage <- 100 * mean(values$status[rows] == "Includes zero")
    data.frame(
      panel = panel, percentage = percentage,
      label = sprintf("Intervals including 0: %.0f%%", percentage),
      x = Inf, y = -Inf, stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, result)
  result$panel <- factor(result$panel, levels = panels)
  rownames(result) <- NULL
  result
}

# This function draws the interval comparison.
.plot_caterpillar <- function(object, observed_se = "gg") {

  values <- .caterpillar_data(object, observed_se)
  zero_summary <- .caterpillar_zero_summary(values)
  level <- object$level
  zero_colour <- "#873B45"
  away_colour <- "#176B60"

  ggplot2::ggplot(values) +
    ggplot2::geom_vline(
      xintercept = 0, colour = "#697076",
      linetype = "dotted", linewidth = 0.5
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data[["conf.low"]], xend = .data[["conf.high"]],
                   y = .data[["rank"]], yend = .data[["rank"]],
                   colour = .data[["status"]]),
      linewidth = 0.6, alpha = 0.78, lineend = "round"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data[["estimate"]], y = .data[["rank"]],
                   colour = .data[["status"]]),
      size = 1.8, alpha = 0.95
    ) +
    ggplot2::geom_text(
      data = zero_summary,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                   label = .data[["label"]]),
      inherit.aes = FALSE, hjust = 1.08, vjust = -0.8,
      colour = "#687076", size = 3.05
    ) +
    ggplot2::geom_blank(
      data = unique(data.frame(task = values$task, method = values$method,
                               x = 0, y = 1)),
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
      inherit.aes = FALSE
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data[["panel"]]),
      ncol = 2, scales = "free_y"
    ) +
    ggplot2::scale_colour_manual(
      values = c("Includes zero" = zero_colour,
                 "Excludes zero" = away_colour),
      name = sprintf("%.0f%% interval", 100 * level)
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.05)) +
    ggplot2::scale_y_continuous(
      breaks = NULL, expand = ggplot2::expansion(mult = c(0.045, 0.015))
    ) +
    ggplot2::labs(
      x = expression(Sensitivity ~ (italic(d) * minute)),
      y = NULL,
      caption = paste0(
        "Subjects are ordered by observed d' within each task. ",
        if (.is_bayes(object))
          "Model estimates are posterior means with central credible intervals. ",
        "Percentages descriptively summarise the displayed intervals."
      )
    ) +
    .panel_theme()
}

# This function gives the theme every uSDT plot starts from: classic axes in
# dark grey, panel titles in plain text on the left over a thin grey line, the
# legend on top and the caption under the plot.
.usdt_theme <- function() {
  ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.spacing.x = grid::unit(18, "pt"),
      panel.spacing.y = grid::unit(14, "pt"),
      axis.line = ggplot2::element_line(colour = "#292D30", linewidth = 0.65),
      axis.ticks = ggplot2::element_line(colour = "#292D30", linewidth = 0.5),
      axis.text = ggplot2::element_text(colour = "#3E4347", size = 10),
      axis.title = ggplot2::element_text(colour = "#292D30", size = 13.5),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
      strip.background = .strip_rule(),
      strip.text = ggplot2::element_text(
        colour = "#252A2D", size = 12.5, hjust = 0,
        margin = ggplot2::margin(2, 0, 6, 0)
      ),
      plot.title = ggplot2::element_text(colour = "#252A2D", size = 13.5),
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 10.5),
      legend.text = ggplot2::element_text(size = 10),
      legend.key.width = grid::unit(18, "pt")
    ) +
    .caption_theme(c(14, 20, 14, 16))
}

# This function gives the background of the panel titles: instead of a box, a
# thin grey line along its bottom edge, as wide as the panel, that separates
# the title from the panel under it.
.strip_rule <- function() {
  structure(
    list(fill = NA, colour = "#B9C0C5", linewidth = 0.45, linetype = 1,
         linejoin = "mitre", inherit.blank = FALSE),
    class = c("usdt_rule", "element_rect", "element")
  )
}

#' @exportS3Method ggplot2::element_grob
element_grob.usdt_rule <- function(element, ...) {
  grid::segmentsGrob(0, 0, 1, 0, gp = grid::gpar(
    col = element$colour, lwd = element$linewidth * ggplot2::.pt,
    lineend = "butt"))
}

# This function gives the theme of the faceted plots whose y axis carries no
# scale of its own: subject ranks in the caterpillar plot and densities in the
# prior predictive check and the Bayes factor plot.
.panel_theme <- function() {
  .usdt_theme() +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank()
    )
}

# This function collects the parameters that define both ROC curves.
.roc_parameters <- function(object, subject_id = NULL, conditional = FALSE) {

  # Four common names keep estimated and fixed criteria interchangeable.
  parameters <- c("c_A", "c_I", "d_A", "d_I")
  se <- stats::setNames(rep.int(NA_real_, length(parameters)), parameters)

  # The population curve uses fixed effects and their analytic sensitivity SEs.
  if (is.null(subject_id)) {
    fixed_se <- sqrt(diag(as.matrix(stats::vcov(object$fit))))
    estimated <- intersect(parameters, names(fixed_se))
    se[estimated] <- fixed_se[estimated]
    population <- .task_effects(object$fit, parameters = parameters)$population
    return(list(subject = NULL, estimate = population, se = se))
  }

  # One checked identifier selects the subject's conditional modes.
  if (length(subject_id) != 1L || is.na(subject_id)) {
    .usdt_stop("`subject_id` must identify exactly one subject.")
  }
  subject <- as.character(subject_id)
  estimate <- .task_effects(object$fit, subject, parameters)$subjects[1L, ]
  if (anyNA(estimate)) {
    .usdt_stop("subject `", subject, "` was not found in the fitted model.")
  }

  # Conditional standard deviations describe uncertainty around this mode.
  if (conditional) {
    values <- as.data.frame(lme4::ranef(object$fit, condVar = TRUE))
    values <- values[values$grpvar == "subj" &
                       as.character(values$grp) == subject, ]
    value_row <- match(parameters, values$term)
    available <- !is.na(value_row)
    se[available] <- values$condsd[value_row[available]]
  }
  list(subject = subject, estimate = estimate, se = se)
}

# This function describes which task curves are observed for one plot.
.roc_tasks <- function(object, subject = NULL) {

  labels <- object$data$meta$labels
  tasks <- data.frame(
    code = c("A", "I"),
    task = unname(labels[c("awareness", "influence")]),
    criterion = c("c_A", "c_I"),
    sensitivity = c("d_A", "d_I"),
    stringsAsFactors = FALSE
  )
  if (is.null(subject)) return(tasks)

  # A subject curve is shown only for tasks that supplied observations.
  observed <- unique(as.character(object$data$agg$task[
    as.character(object$data$agg$subj) == subject]))
  tasks <- tasks[tasks$code %in% observed, , drop = FALSE]
  if (!nrow(tasks)) {
    .usdt_stop("subject `", subject, "` has no observations in the fitted data.")
  }
  if (nrow(tasks) < 2L) {
    .usdt_warn("subject `", subject, "` appears in only one task; only that ROC curve is shown.")
  }
  tasks
}

# This function transforms a sensitivity interval into an ROC band.
.roc_analytic_band <- function(far, estimate, se, level) {

  critical <- stats::qnorm(1 - (1 - level) / 2)
  limits <- estimate + c(-1, 1) * critical * se
  zfar <- stats::qnorm(far)
  cbind(low = stats::pnorm(zfar + limits[1L]),
        high = stats::pnorm(zfar + limits[2L]))
}

# This function transforms bootstrap sensitivities into an ROC interval.
.roc_bootstrap_band <- function(far, estimate, samples, level, type) {

  fitted <- stats::pnorm(stats::qnorm(far) + estimate)
  simulated <- stats::pnorm(outer(samples, stats::qnorm(far), "+"))
  interval <- pmin(pmax(.boot_ci(simulated, fitted, level, type), 0), 1)
  colnames(interval) <- c("low", "high")
  interval
}

# This function creates the curves, intervals and criterion points.
.roc_data <- function(object, subject_id = NULL, band = TRUE) {

  if (.is_bayes(object)) return(.roc_data_bayes(object, subject_id, band))
  subject_fit <- .roc_parameters(object, subject_id, conditional = band)
  tasks <- .roc_tasks(object, subject_fit$subject)
  level <- object$level
  far <- seq(0, 1, length.out = 201L)
  coding <- object$data$meta$coding
  boot_band <- is.null(subject_fit$subject) && band &&
    .has_boot_summary(object)

  # Every task uses the same equal-variance ROC transformation.
  curves <- vector("list", nrow(tasks))
  points <- vector("list", nrow(tasks))
  for (i in seq_len(nrow(tasks))) {
    task <- tasks[i, ]
    dprime <- unname(subject_fit$estimate[[task$sensitivity]])
    criterion <- unname(subject_fit$estimate[[task$criterion]])
    auc <- stats::pnorm(dprime / sqrt(2))
    curve_label <- sprintf("%s: d' = %.2f | AUC = %.2f",
                           task$task, dprime, auc)
    hit <- stats::pnorm(stats::qnorm(far) + dprime)
    interval <- matrix(NA_real_, nrow = length(far), ncol = 2L,
                       dimnames = list(NULL, c("low", "high")))
    method <- "none"

    if (boot_band) {
      samples <- object$boot$population$t[object$boot$ok,
                                           task$sensitivity]
      interval <- .roc_bootstrap_band(far, dprime, samples, level,
                                      object$boot$type)
      interval_name <- switch(object$boot$type,
                              perc = "percentile", norm = "normal",
                              basic = "basic")
      method <- paste("bootstrap", interval_name)
    } else if (band) {
      se <- unname(subject_fit$se[[task$sensitivity]])
      if (!is.finite(se)) {
        .usdt_stop("the fitted model does not provide an uncertainty estimate for `",
                   task$sensitivity, "`.")
      }
      interval <- .roc_analytic_band(far, dprime, se, level)
      method <- if (is.null(subject_fit$subject)) "Wald" else "conditional"
    }

    # The criterion chooses one operating point on the model-implied curve.
    linear_noise <- if (coding == "deviation") criterion - dprime / 2 else criterion
    linear_signal <- if (coding == "deviation") criterion + dprime / 2 else
      criterion + dprime
    point <- data.frame(
      task = task$task, curve = curve_label,
      false_alarm = stats::pnorm(linear_noise),
      hit = stats::pnorm(linear_signal), criterion = -criterion,
      stringsAsFactors = FALSE
    )
    curves[[i]] <- data.frame(
      task = task$task, curve = curve_label, false_alarm = far,
      hit = hit, conf.low = interval[, "low"], conf.high = interval[, "high"],
      dprime = dprime, auc = auc, criterion = point$criterion,
      operating_false_alarm = point$false_alarm, operating_hit = point$hit,
      interval = method, subject = subject_fit$subject %||% "Population",
      stringsAsFactors = FALSE
    )
    points[[i]] <- point
  }
  .roc_bundle(curves, points, tasks, subject_fit$subject)
}

# This function joins the curves and operating points of every task.
.roc_bundle <- function(curves, points, tasks, subject) {
  curves <- do.call(rbind, curves)
  points <- do.call(rbind, points)
  curves$curve <- factor(curves$curve, levels = unique(curves$curve))
  points$curve <- factor(points$curve, levels = levels(curves$curve))
  rownames(curves) <- NULL
  rownames(points) <- NULL
  list(curves = curves, points = points, tasks = tasks, subject = subject)
}

# This function builds the ROC curves of a Bayesian fit draw by draw. Each draw
# gives one curve and one operating point; the plot shows their posterior means
# with central credible bands. Under unequal variances each curve bends with
# the signal SD of its draw.
.roc_data_bayes <- function(object, subject_id = NULL, band = TRUE) {

  subject <- NULL
  if (!is.null(subject_id)) {
    if (length(subject_id) != 1L || is.na(subject_id)) {
      .usdt_stop("`subject_id` must identify exactly one subject.")
    }
    subject <- as.character(subject_id)
    if (!subject %in% object$subjects) {
      .usdt_stop("subject `", subject, "` was not found in the fitted model.")
    }
  }
  tasks <- .roc_tasks(object, subject)
  far <- seq(0, 1, length.out = 201L)
  a <- (1 - object$level) / 2
  free <- object$design$free_c
  population <- .population_draws(object)
  draw <- function(name) {
    as.numeric(posterior::extract_variable(population, name))
  }

  curves <- vector("list", nrow(tasks))
  points <- vector("list", nrow(tasks))
  for (i in seq_len(nrow(tasks))) {
    task <- tasks[i, ]
    j <- match(task$code, c("A", "I"))

    # The population curve uses the means and a subject curve its own values.
    if (is.null(subject)) {
      d <- draw(sprintf("mu_d[%d]", j))
      criterion <- if (free[j] == 1L) draw(sprintf("mu_c[%d]", sum(free[seq_len(j)])))
    } else {
      d <- as.vector(.subject_draws(object, "d", j, subject))
      criterion <- as.vector(.subject_draws(object, "c", j, subject))
    }
    s <- if (object$design$unequal_variances) draw(sprintf("sigma_s[%d]", j)) else
      rep(1, length(d))

    # A criterion fixed by the Meyen split makes HR + FAR = 1.
    if (is.null(criterion)) criterion <- d * (1 - s) / (2 * (1 + s))

    hits <- stats::pnorm(outer(d, stats::qnorm(far), "+") / s)
    interval <- matrix(NA_real_, nrow = length(far), ncol = 2L)
    if (band) {
      interval <- t(apply(hits, 2L, stats::quantile, probs = c(a, 1 - a),
                          names = FALSE))
    }
    dprime <- mean(d)
    auc <- mean(stats::pnorm(d / sqrt(1 + s^2)))
    curve_label <- sprintf("%s: d' = %.2f | AUC = %.2f",
                           task$task, dprime, auc)
    point <- data.frame(
      task = task$task, curve = curve_label,
      false_alarm = mean(stats::pnorm(-d / 2 - criterion)),
      hit = mean(stats::pnorm((d / 2 - criterion) / s)),
      criterion = mean(criterion), stringsAsFactors = FALSE
    )
    curves[[i]] <- data.frame(
      task = task$task, curve = curve_label, false_alarm = far,
      hit = colMeans(hits), conf.low = interval[, 1L],
      conf.high = interval[, 2L], dprime = dprime, auc = auc,
      criterion = point$criterion,
      operating_false_alarm = point$false_alarm, operating_hit = point$hit,
      interval = if (band) "credible" else "none",
      subject = subject %||% "Population", stringsAsFactors = FALSE
    )
    points[[i]] <- point
  }
  .roc_bundle(curves, points, tasks, subject)
}

# This function draws the model-implied ROC comparison.
.plot_roc <- function(object, subject_id, band, population_reference) {

  values <- .roc_data(object, subject_id, band)
  curves <- values$curves
  points <- values$points
  reference <- NULL
  if (!is.null(values$subject) && population_reference) {
    reference <- .roc_data(object, NULL, band = FALSE)$curves
    reference <- reference[reference$task %in% values$tasks$task, , drop = FALSE]
  }

  task_colours <- stats::setNames(
    c("#176B60", "#B6543A"),
    unname(object$data$meta$labels[c("awareness", "influence")])
  )
  curve_tasks <- unique(curves[c("curve", "task")])
  colours <- task_colours[curve_tasks$task]
  names(colours) <- as.character(curve_tasks$curve)
  plot <- ggplot2::ggplot(curves)
  if (band) {
    plot <- plot + ggplot2::geom_ribbon(
      ggplot2::aes(x = .data[["false_alarm"]], ymin = .data[["conf.low"]],
                   ymax = .data[["conf.high"]], fill = .data[["curve"]]),
      alpha = 0.15, colour = NA, show.legend = FALSE
    )
  }
  if (!is.null(reference)) {
    reference$curve <- factor(reference$task, levels = values$tasks$task,
                              labels = levels(curves$curve))
    plot <- plot + ggplot2::geom_line(
      data = reference,
      ggplot2::aes(x = .data[["false_alarm"]], y = .data[["hit"]],
                   colour = .data[["curve"]]),
      inherit.aes = FALSE, linewidth = 0.55, linetype = "dashed",
      alpha = 0.65, show.legend = FALSE
    )
  }

  # Only a subject's plot has a title; the caption names the population.
  title <- if (!is.null(values$subject)) paste("Subject", values$subject)
  variance <- if (isTRUE(object$design$unequal_variances)) "unequal" else "equal"
  caption <- paste0(
    if (is.null(values$subject)) "Population ROC curves" else "ROC curves",
    " under ", variance, "-variance SDT",
    if (band) sprintf(", with %.0f%% %s bands", 100 * object$level,
                      curves$interval[1L]),
    "; points mark the fitted criteria.",
    if (!is.null(reference)) " Dashed lines show the population curves."
  )

  # The legend sits inside the square, in the corner the curves keep away
  # from: under the diagonal unless the sensitivities lean negative.
  dprime <- range(curves$dprime)
  corner <- if (sum(dprime) >= 0) c(1, 0) else c(0, 1)

  plot +
    ggplot2::geom_abline(
      intercept = 0, slope = 1, colour = "#8A9196",
      linewidth = 0.55, linetype = "dotted"
    ) +
    ggplot2::geom_line(
      ggplot2::aes(x = .data[["false_alarm"]], y = .data[["hit"]],
                   colour = .data[["curve"]]),
      linewidth = 1.05, lineend = "round"
    ) +
    ggplot2::geom_point(
      data = points,
      ggplot2::aes(x = .data[["false_alarm"]], y = .data[["hit"]],
                   colour = .data[["curve"]]),
      inherit.aes = FALSE, shape = 21, fill = "white",
      size = 3.2, stroke = 1.05
    ) +
    ggplot2::scale_colour_manual(values = colours, name = NULL) +
    (if (band) ggplot2::scale_fill_manual(values = colours, name = NULL)) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
      expand = ggplot2::expansion(mult = 0.01)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
      expand = ggplot2::expansion(mult = 0.01)
    ) +
    ggplot2::coord_equal() +
    ggplot2::annotation_custom(.roc_legend(curves, colours, corner)) +
    ggplot2::labs(
      title = title, x = "False-alarm rate", y = "Hit rate", caption = caption
    ) +
    .usdt_theme() +
    ggplot2::theme(aspect.ratio = 1, legend.position = "none")
}

# This function draws the legend of the ROC plot as a small table, so that the
# d' and the AUC of the tasks line up in columns under their headers: each row
# has the key of a task (its line and open circle), its name and its values.
# The table sits in `corner` of the panel: c(1, 0) is the bottom right.
.roc_legend <- function(curves, colours, corner) {
  rows <- unique(curves[c("curve", "task", "dprime", "auc")])
  text <- function(label, right = FALSE, colour = "#252A2D") {
    grid::textGrob(label, x = as.numeric(right), hjust = as.numeric(right),
                   gp = grid::gpar(fontsize = 10, col = colour))
  }
  key <- function(colour) {
    stroke <- 1.05 * ggplot2::.stroke / 2
    grid::gTree(children = grid::gList(
      grid::segmentsGrob(0, 0.5, 1, 0.5, gp = grid::gpar(
        col = colour, lwd = 1.05 * ggplot2::.pt, lineend = "round")),
      grid::pointsGrob(0.5, 0.5, pch = 21, gp = grid::gpar(
        col = colour, fill = "white", lwd = stroke,
        fontsize = 3.2 * ggplot2::.pt + stroke))))
  }
  space <- function(left = 0, right = 0) grid::unit(c(2, left, 2, right), "pt")

  # The values align on the right under headers in the grey of the captions.
  # A frame grows one column at a time, so the header row starts with empty
  # cells over the keys and the names.
  frame <- grid::frameGrob()
  frame <- grid::packGrob(frame, grid::nullGrob(), row = 1L, col = 1L)
  frame <- grid::packGrob(frame, grid::nullGrob(), row = 1L, col = 2L)
  frame <- grid::packGrob(frame, text(quote(italic(d) * minute), TRUE,
                                      "#687076"),
                          row = 1L, col = 3L, border = space(12))
  frame <- grid::packGrob(frame, text("AUC", TRUE, "#687076"),
                          row = 1L, col = 4L, border = space(12))
  for (i in seq_len(nrow(rows))) {
    r <- i + 1L
    frame <- grid::packGrob(frame, key(colours[[as.character(rows$curve[i])]]),
                            row = r, col = 1L, width = grid::unit(22, "pt"),
                            height = grid::unit(14, "pt"), border = space(0, 6))
    frame <- grid::packGrob(frame, text(rows$task[i]), row = r, col = 2L,
                            border = space())
    frame <- grid::packGrob(frame, text(sprintf("%.2f", rows$dprime[i]), TRUE),
                            row = r, col = 3L, border = space(12))
    frame <- grid::packGrob(frame, text(sprintf("%.2f", rows$auc[i]), TRUE),
                            row = r, col = 4L, border = space(12))
  }

  # annotation_custom() gives the outer tree the viewport of the panel, so the
  # table takes its corner in a viewport of its own inside it.
  inset <- grid::unit((1 - 2 * corner) * 6, "pt")
  table <- grid::gTree(children = grid::gList(frame), vp = grid::viewport(
    x = grid::unit(corner[1L], "npc") + inset[1L],
    y = grid::unit(corner[2L], "npc") + inset[2L],
    width = grid::grobWidth(frame), height = grid::grobHeight(frame),
    just = corner))
  grid::gTree(children = grid::gList(table))
}

# This function creates the points that form each contour.
.ellipse_points <- function(x, y, levels, n = 121L) {

  values <- stats::na.omit(cbind(x = x, y = y))
  if (nrow(values) < 3L) return(data.frame())
  covariance <- stats::cov(values)
  decomposition <- eigen(covariance, symmetric = TRUE)
  spread <- pmax(decomposition$values, 0)
  if (!all(is.finite(spread)) || max(spread) <= 0) return(data.frame())

  # A circle takes the shape and direction of the values.
  angle <- seq(0, 2 * pi, length.out = n)
  circle <- cbind(cos(angle), sin(angle))
  transform <- diag(sqrt(spread), nrow = 2L) %*% t(decomposition$vectors)
  centre <- colMeans(values)

  result <- lapply(levels, function(level) {
    points <- sqrt(stats::qchisq(level, df = 2L)) * circle %*% transform
    points <- sweep(points, 2L, centre, "+")
    data.frame(x = points[, 1L], y = points[, 2L],
               coverage = factor(level, levels = levels))
  })
  do.call(rbind, result)
}
