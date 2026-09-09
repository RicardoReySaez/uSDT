# plot.R
# Visualize results from fitted hierarchical SDT models
# Author: Ricardo Rey-Sáez
# Last modified: 08-09-2026

#' Diagnostic and analytical plots for hierarchical SDT models
#'
#' Generates diagnostic visualizations for a fitted hierarchical SDT model,
#' including observed versus latent regression (H3), shrinkage patterns,
#' participant-level caterpillar intervals, and model-implied ROC curves.
#'
#' @param x An `hsdt` object fitted by [hsdt()].
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
#' @param ... Additional arguments passed to underlying plotting methods.
#'
#' @return A `ggplot` object. Its underlying data frame is stored in `$data`.
#'
#' @details
#' # Regression plot (`type = "regression"`)
#'
#' Compares observed and latent associations across two panels sharing axes.
#' The left panel shows observed \eqn{d'} values and an ordinary least-squares
#' line. When direct task reliability is low, trial-level sampling noise
#' attenuates this observed slope toward zero. The right panel plots the latent
#' regression line (\eqn{d'_I} on \eqn{d'_D}) from H3, correcting for
#' measurement error. The value of each line at \eqn{d'_D = 0} marks the
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
#' @seealso [hsdt()], [sdt_moments()], [usdt_boot()]
#'
#' @examples
#' \donttest{
#' # Contextual cuing data from Vadillo et al. (2025)
#' d <- usdt_data_tasks(
#'   direct   = vadillo_awareness,
#'   indirect = vadillo_cuing,
#'   subject_col      = "subj",
#'   condition_col    = "condition",
#'   condition_levels = c(signal = "old", noise = "new"),
#'   response_col     = list(direct = "judged.old", indirect = "rt"),
#'   response_levels  = list(direct   = c(signal = 1, noise = 0),
#'                           indirect = c(signal = "faster", noise = "slower")),
#'   dichotomize      = list(direct = FALSE, indirect = TRUE)
#' )
#'
#' fit <- hsdt(d)
#'
#' # 1. Observed vs. latent regression (H3)
#' plot(fit, type = "regression")
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
#' }
#'
#' @export
plot.hsdt <- function(x, type = c("regression", "shrinkage",
                                  "caterpillar", "roc"),
                      subject_id = NULL, band = TRUE,
                      population_reference = TRUE, observed_se = NULL, ...) {

  # The function checks the requested plot.
  type <- match.arg(type)
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

  switch(type,
         regression = .plot_regression(x, band),
         shrinkage = .plot_shrinkage(x),
         caterpillar = .plot_caterpillar(x, observed_se),
         roc = .plot_roc(x, subject_id, band, population_reference))
}

# This function breaks a caption into lines that fit inside the plot.
.wrap_caption <- function(..., width = 88L) {

  # ggplot2 lays a caption out at its natural width and never wraps it, so a
  # long one simply runs past the edge of the device and its tail is lost. The
  # width is counted in characters rather than inches because the caption is
  # drawn at a fixed point size: at the 9.5pt these plots use, this many
  # characters stay inside the panels at the figure sizes they are drawn at.
  # A caption that is still too wide for a very narrow figure will overflow,
  # which is why the value is conservative.
  paragraphs <- unlist(strsplit(paste0(...), "\n", fixed = TRUE))
  paragraphs <- paragraphs[nzchar(trimws(paragraphs))]

  # Each paragraph wraps on its own so the author's own breaks survive.
  wrapped <- vapply(
    paragraphs,
    function(p) paste(strwrap(p, width = width), collapse = "\n"),
    character(1L)
  )
  paste(wrapped, collapse = "\n")
}

# This function evaluates the latent line and its band over a grid of x.
.latent_line <- function(object, x, band = TRUE) {

  e <- object$pars$est
  slope <- e[["s_DI"]] / e[["s2_D"]]
  fit <- e[["gamma_I"]] + slope * (x - e[["gamma_D"]])
  limits <- matrix(NA_real_, length(x), 2L)

  # The band follows whichever method the hypothesis table used, and names the
  # reason whenever it cannot be drawn.
  method <- "none"
  reason <- NA_character_
  if (!band) {
    reason <- "it was turned off"
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

  list(line = data.frame(x = x, fit = fit, conf.low = limits[, 1L],
                         conf.high = limits[, 2L], stringsAsFactors = FALSE),
       slope = slope, intercept = e[["gamma_I"]] - slope * e[["gamma_D"]],
       method = method, reason = reason)
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
  spread <- stats::var(points$observed_direct)
  naive_ok <- nrow(points) >= 2L && isTRUE(spread > 0)
  naive_fit <- if (naive_ok) {
    stats::lm(observed_indirect ~ observed_direct, data = points)
  } else NULL
  # Zero joins the observed grid for the same reason it joins the latent one:
  # the value of the line there is the quantity both panels put on trial.
  naive_span <- if (naive_ok) range(c(0, points$observed_direct)) else numeric(0)
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
      naive_fit, newdata = data.frame(observed_direct = naive_x),
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
  # plus zero because its value there is the hypothesis of interest.
  e <- object$pars$est
  latent_span <- if (is.finite(e[["s2_D"]]) && e[["s2_D"]] > 0) {
    e[["gamma_D"]] + stats::qnorm(c(.025, .975)) * sqrt(e[["s2_D"]])
  } else {
    range(points$model_direct)
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
.plot_regression <- function(object, band) {

  values <- .regression_data(object, band)
  labels <- object$data$meta$labels
  origin <- values$origin
  naive_origin <- values$naive_origin
  banded <- values$method != "none"
  naive <- is.finite(values$naive_slope)
  naive_banded <- band && naive && nrow(naive_origin) == 1L &&
    all(is.finite(c(naive_origin$conf.low, naive_origin$conf.high)))

  background <- "#FFFFFF"
  observed_colour <- "#A8AEB3"
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

  # The right panel shows where partial pooling moves every observed subject.
  plot <- plot +
    ggplot2::geom_segment(
      data = values$model_points,
      ggplot2::aes(x = .data[["observed_direct"]],
                   y = .data[["observed_indirect"]],
                   xend = .data[["model_direct"]],
                   yend = .data[["model_indirect"]]),
      inherit.aes = FALSE, colour = observed_colour, alpha = 0.34,
      linewidth = 0.42
    ) +
    ggplot2::geom_point(
      data = values$observed_points,
      ggplot2::aes(x = .data[["observed_direct"]],
                   y = .data[["observed_indirect"]]),
      inherit.aes = FALSE,
      shape = 21, size = 2.4, stroke = 0.7,
      colour = observed_colour, fill = background, alpha = 0.8
    ) +
    ggplot2::geom_point(
      data = values$model_points,
      ggplot2::aes(x = .data[["observed_direct"]],
                   y = .data[["observed_indirect"]]),
      inherit.aes = FALSE,
      shape = 21, size = 2.25, stroke = 0.65,
      colour = observed_colour, fill = background, alpha = 0.42
    ) +
    ggplot2::geom_point(
      data = values$model_points,
      ggplot2::aes(x = .data[["model_direct"]],
                   y = .data[["model_indirect"]]),
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
    sprintf("Bands are %.0f%% OLS and bootstrap %s confidence intervals.",
            100 * object$level,
            switch(values$method, perc = "percentile", norm = "normal",
                   basic = "basic")))
  caption <- .wrap_caption(
    "Segments show shrinkage from observed to conditional model estimates. ",
    "The latent line is implied by the fitted random-effects distribution, ",
    "not fitted to the green points.\n",
    "An open circle marks the intercept of each panel at zero. ", band_note,
    " P-values test the intercept and slope against zero; p < .05 indicates ",
    "significance."
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
    ggplot2::facet_wrap(~panel, nrow = 1L) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.06)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.06)) +
    ggplot2::labs(
      x = bquote(.(labels[["direct"]]) ~ sensitivity ~ (italic(d) * minute)),
      y = bquote(.(labels[["indirect"]]) ~ sensitivity ~ (italic(d) * minute)),
      caption = caption
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.background = ggplot2::element_rect(fill = background, colour = NA),
      plot.background = ggplot2::element_rect(fill = background, colour = NA),
      axis.line = ggplot2::element_line(colour = "#292D30", linewidth = 0.72),
      axis.ticks = ggplot2::element_line(colour = "#292D30", linewidth = 0.55),
      axis.text = ggplot2::element_text(colour = "#3E4347", size = 10.5),
      axis.title = ggplot2::element_text(colour = "#292D30", size = 14.5),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 12)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12)),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(colour = "#292D30", size = 14.5),
      panel.spacing = grid::unit(24, "pt"),
      plot.caption = ggplot2::element_text(
        colour = "#687076", size = 9.5, hjust = 0,
        margin = ggplot2::margin(t = 10)
      ),
      plot.margin = ggplot2::margin(18, 22, 18, 18)
    )
}

# This function pairs the observed sensitivities of every complete subject.
.observed_pairs <- function(object) {

  # The observed values come from the same prepared data as the model.
  observed <- sdt_moments(object$data, correction = "hautus")
  labels <- object$data$meta$labels
  direct <- observed[observed$task == labels[["direct"]], c("subj", "dprime")]
  indirect <- observed[observed$task == labels[["indirect"]], c("subj", "dprime")]
  direct$subj <- as.character(direct$subj)
  indirect$subj <- as.character(indirect$subj)

  # Both observed values place a subject on the plot.
  indirect_row <- match(direct$subj, indirect$subj)
  complete <- !is.na(indirect_row)
  omitted <- length(unique(c(direct$subj, indirect$subj))) - sum(complete)
  if (omitted > 0L) {
    .usdt_warn(omitted, " subject", if (omitted == 1L) " was" else "s were",
               " omitted because both observed sensitivities are required.")
  }

  data.frame(
    subject = direct$subj[complete],
    observed_direct = direct$dprime[complete],
    observed_indirect = indirect$dprime[indirect_row[complete]],
    stringsAsFactors = FALSE
  )
}

# This function adds the conditional model value of each subject.
.shrinkage_data <- function(object) {

  values <- .observed_pairs(object)

  # Each model value combines the task average with the subject difference.
  random <- lme4::ranef(object$fit, condVar = FALSE)[["subj"]]
  fixed <- lme4::fixef(object$fit)
  random_row <- match(values$subject, rownames(random))
  if (anyNA(random_row) || !all(c("d_D", "d_I") %in% names(random))) {
    .usdt_stop("the fitted model does not contain the expected subject sensitivities.")
  }
  values$model_direct <- fixed[["d_D"]] + random$d_D[random_row]
  values$model_indirect <- fixed[["d_I"]] + random$d_I[random_row]
  values
}

# This function draws the shrinkage comparison.
.plot_shrinkage <- function(object) {

  values <- .shrinkage_data(object)
  labels <- object$data$meta$labels

  # The contours show the shape of both groups of points.
  observed_ellipse <- .ellipse_points(
    values$observed_direct, values$observed_indirect,
    levels = c(0.95, 0.80, 0.50)
  )
  model_ellipse <- .ellipse_points(
    values$model_direct, values$model_indirect,
    levels = 0.90
  )

  # Warm grey and muted teal separate raw values from model estimates.
  background <- "#FFFFFF"
  observed_colour <- "#A8AEB3"
  model_colour <- "#116B60"
  model_fill <- "#3FA88E"

  plot <- ggplot2::ggplot(values)
  if (nrow(observed_ellipse)) {
    plot <- plot + ggplot2::geom_polygon(
      data = observed_ellipse,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                   group = .data[["coverage"]]),
      inherit.aes = FALSE,
      fill = observed_colour, alpha = 0.06, colour = NA
    )
  }
  if (nrow(model_ellipse)) {
    plot <- plot + ggplot2::geom_polygon(
      data = model_ellipse,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                   group = .data[["coverage"]]),
      inherit.aes = FALSE,
      fill = model_fill, colour = model_colour,
      alpha = 0.17, linewidth = 0.4
    )
  }

  plot +
    ggplot2::geom_hline(
      yintercept = 0, colour = "#555B60",
      linetype = "dotted", linewidth = 0.5
    ) +
    ggplot2::geom_vline(
      xintercept = 0, colour = "#555B60",
      linetype = "dotted", linewidth = 0.5
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data[["observed_direct"]],
                   y = .data[["observed_indirect"]],
                   xend = .data[["model_direct"]],
                   yend = .data[["model_indirect"]]),
      colour = observed_colour, alpha = 0.48, linewidth = 0.45
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data[["observed_direct"]],
                   y = .data[["observed_indirect"]]),
      shape = 21, size = 2.55, stroke = 0.75,
      colour = observed_colour, fill = background, alpha = 0.75
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data[["model_direct"]],
                   y = .data[["model_indirect"]]),
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
      x = bquote(.(labels[["direct"]]) ~ sensitivity ~ (italic(d) * minute)),
      y = bquote(.(labels[["indirect"]]) ~ sensitivity ~ (italic(d) * minute))
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.background = ggplot2::element_rect(fill = background, colour = NA),
      plot.background = ggplot2::element_rect(fill = background, colour = NA),
      axis.line = ggplot2::element_line(colour = "#292D30", linewidth = 0.72),
      axis.ticks = ggplot2::element_line(colour = "#292D30", linewidth = 0.55),
      axis.text = ggplot2::element_text(colour = "#3E4347", size = 10.5),
      axis.title = ggplot2::element_text(colour = "#292D30", size = 14.5),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 12)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12)),
      legend.position = "none",
      plot.margin = ggplot2::margin(18, 22, 18, 18)
    )
}

# This function collects both intervals for each subject and task.
.caterpillar_data <- function(object, observed_se = "gg") {

  # The whole interval needs reliable joint uncertainty.
  if (!isTRUE(object$pars$joint_ok)) {
    .usdt_stop("subject intervals are unavailable because ",
               object$pars$inference_reason, ".")
  }

  # The selected variance gives the trial error around each observed value.
  se_col <- paste0("se_", observed_se)
  observed <- sdt_moments(object$data, correction = "hautus", variances = TRUE)
  labels <- object$data$meta$labels
  level <- object$level
  critical <- stats::qnorm(1 - (1 - level) / 2)

  # The model supplies each subject estimate and its uncertainty.
  conditional <- .conditional_se(object$fit, object$pars$V_full, object$devfun)
  conditional <- conditional[
    conditional$grpvar == "subj" & conditional$term %in% c("d_D", "d_I"), ]
  fixed <- lme4::fixef(object$fit)

  # Each task keeps the order given by its observed values.
  tasks <- list(
    direct = list(label = labels[["direct"]], term = "d_D"),
    indirect = list(label = labels[["indirect"]], term = "d_I")
  )
  result <- lapply(tasks, function(task) {
    raw <- observed[observed$task == task$label, ]
    raw$subj <- as.character(raw$subj)
    model <- conditional[conditional$term == task$term, ]
    model_row <- match(raw$subj, as.character(model$grp))
    if (anyNA(model_row) || any(!is.finite(model$se[model_row]))) {
      .usdt_stop("the fitted model does not provide conditional intervals for every subject.")
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
    model_rows <- data.frame(
      subject = raw$subj,
      task = task$label,
      method = "Model-estimated",
      estimate = fixed[[task$term]] + model$condval[model_row],
      se = model$se[model_row],
      rank = rank,
      stringsAsFactors = FALSE
    )
    rbind(observed_rows, model_rows)
  })

  # The interval colour shows whether zero remains plausible.
  values <- do.call(rbind, result)
  values$conf.low <- values$estimate - critical * values$se
  values$conf.high <- values$estimate + critical * values$se
  includes_zero <- values$conf.low <= 0 & values$conf.high >= 0
  values$status <- factor(
    ifelse(includes_zero, "Includes zero", "Excludes zero"),
    levels = c("Includes zero", "Excludes zero")
  )
  values$task <- factor(values$task,
                        levels = c(labels[["direct"]], labels[["indirect"]]))
  values$method <- factor(values$method,
                          levels = c("Observed", "Model-estimated"))
  # Panel titles are drawn on a graphics device, so they stay ASCII: a device
  # in a single-byte locale cannot convert every character, and a failed
  # conversion is an error at render time rather than a missing glyph.
  panel_levels <- c(
    paste0(labels[["direct"]], " task: Observed"),
    paste0(labels[["direct"]], " task: Model-estimated"),
    paste0(labels[["indirect"]], " task: Observed"),
    paste0(labels[["indirect"]], " task: Model-estimated")
  )
  values$panel <- factor(
    paste0(values$task, " task: ", values$method),
    levels = panel_levels
  )
  rownames(values) <- NULL
  values
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
      caption = .wrap_caption(
        "Subjects are ordered by observed d' within each task. ",
        "Percentages descriptively summarise the displayed intervals."
      )
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.spacing.x = grid::unit(18, "pt"),
      panel.spacing.y = grid::unit(14, "pt"),
      axis.line = ggplot2::element_line(colour = "#292D30", linewidth = 0.65),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(colour = "#292D30", linewidth = 0.5),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = "#3E4347", size = 10),
      axis.text.y = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(colour = "#292D30", size = 13.5),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
      strip.placement = "inside",
      strip.background = ggplot2::element_rect(fill = "#F5F7F6", colour = NA),
      strip.text = ggplot2::element_text(
        colour = "#252A2D", face = "bold", size = 11,
        margin = ggplot2::margin(7, 8, 7, 8)
      ),
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold", size = 10.5),
      legend.text = ggplot2::element_text(size = 10),
      legend.key.width = grid::unit(18, "pt"),
      plot.caption = ggplot2::element_text(
        colour = "#687076", size = 9.5, hjust = 0,
        margin = ggplot2::margin(t = 10)
      ),
      text = ggplot2::element_text(family = "sans"),
      plot.margin = ggplot2::margin(14, 20, 14, 16)
    )
}

# This function collects the parameters that define both ROC curves.
.roc_parameters <- function(object, subject_id = NULL, conditional = FALSE) {

  # Four common names keep estimated and fixed criteria interchangeable.
  parameters <- c("c_D", "c_I", "d_D", "d_I")
  fixed <- stats::setNames(rep.int(0, length(parameters)), parameters)
  beta <- lme4::fixef(object$fit)
  in_model <- intersect(parameters, names(beta))
  fixed[in_model] <- beta[in_model]

  # The population curve uses fixed effects and their analytic sensitivity SEs.
  if (is.null(subject_id)) {
    fixed_se <- sqrt(diag(as.matrix(stats::vcov(object$fit))))
    se <- stats::setNames(rep.int(NA_real_, length(parameters)), parameters)
    estimated <- intersect(parameters, names(fixed_se))
    se[estimated] <- fixed_se[estimated]
    return(list(subject = NULL, estimate = fixed, se = se))
  }

  # One checked identifier selects the subject's conditional modes.
  if (length(subject_id) != 1L || is.na(subject_id)) {
    .usdt_stop("`subject_id` must identify exactly one subject.")
  }
  subject <- as.character(subject_id)
  random <- lme4::ranef(object$fit, condVar = FALSE)[["subj"]]
  subject_row <- match(subject, rownames(random))
  if (is.na(subject_row)) {
    .usdt_stop("subject `", subject, "` was not found in the fitted model.")
  }
  random_parameters <- intersect(parameters, colnames(random))
  estimate <- fixed
  estimate[random_parameters] <- estimate[random_parameters] +
    unlist(random[subject_row, random_parameters, drop = FALSE], use.names = FALSE)

  # Conditional standard deviations describe uncertainty around this mode.
  se <- stats::setNames(rep.int(NA_real_, length(parameters)), parameters)
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
    code = c("D", "I"),
    task = unname(labels[c("direct", "indirect")]),
    criterion = c("c_D", "c_I"),
    sensitivity = c("d_D", "d_I"),
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

  curves <- do.call(rbind, curves)
  points <- do.call(rbind, points)
  curves$curve <- factor(curves$curve, levels = unique(curves$curve))
  points$curve <- factor(points$curve, levels = levels(curves$curve))
  rownames(curves) <- NULL
  rownames(points) <- NULL
  list(curves = curves, points = points, tasks = tasks,
       subject = subject_fit$subject)
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
    unname(object$data$meta$labels[c("direct", "indirect")])
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

  title <- if (is.null(values$subject)) "Population ROC curves" else
    paste0("ROC curves for subject ", values$subject)
  interval <- unique(curves$interval)
  caption <- if (!band) {
    "Curves assume equal-variance SDT; points mark the fitted criteria."
  } else if (is.null(values$subject)) {
    paste0(sprintf("%.0f%% %s bands; points mark the fitted criteria.",
                   100 * object$level, interval[1L]),
           "\nCurves assume equal-variance SDT.")
  } else {
    paste0(sprintf("%.0f%% conditional bands; points mark the fitted criteria.",
                   100 * object$level),
           if (population_reference) "\nDashed lines show the population curves." else "",
           " Equal-variance SDT.")
  }

  # The square panel makes this the narrowest plot in the package, so its
  # caption wraps sooner than the others.
  caption <- .wrap_caption(caption, width = 58L)

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
    # A fixed aspect ratio makes the panel square and therefore narrower than
    # the device. The legend is laid out at its own natural width, so two
    # entries side by side overflow that panel and spill past the background
    # of the plot. Stacking them keeps the legend as wide as its widest single
    # label, which fits whatever width the square panel ends up with.
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1L)) +
    ggplot2::labs(
      title = title, x = "False-alarm rate", y = "Hit rate", caption = caption
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      axis.line = ggplot2::element_line(colour = "#292D30", linewidth = 0.65),
      axis.ticks = ggplot2::element_line(colour = "#292D30", linewidth = 0.5),
      axis.text = ggplot2::element_text(colour = "#3E4347", size = 10),
      axis.title = ggplot2::element_text(colour = "#292D30", size = 13),
      plot.title = ggplot2::element_text(colour = "#252A2D", face = "bold", size = 14),
      legend.position = "top",
      legend.text = ggplot2::element_text(size = 10),
      legend.key.width = grid::unit(22, "pt"),
      plot.caption = ggplot2::element_text(
        colour = "#687076", size = 9.5, hjust = 0,
        margin = ggplot2::margin(t = 10)
      ),
      text = ggplot2::element_text(family = "sans"),
      plot.margin = ggplot2::margin(14, 20, 14, 16)
    )
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
