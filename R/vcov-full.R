# vcov-full.R
# This script estimates uncertainty for fixed and random model parameters.
# Author: Ricardo Rey-Sáez
# Last modified: 04-09-2026

# Numerical derivatives

# The deviance and the conditional modes need numerical derivatives.

# This function estimates the curvature of one result.
.num_hess <- function(f, x, h = 1e-4) {
  k <- length(x)
  H <- matrix(0, k, k)
  for (i in seq_len(k)) {
    for (j in i:k) {
      e1 <- e2 <- e3 <- e4 <- x
      e1[i] <- e1[i] + h; e1[j] <- e1[j] + h
      e2[i] <- e2[i] + h; e2[j] <- e2[j] - h
      e3[i] <- e3[i] - h; e3[j] <- e3[j] + h
      e4[i] <- e4[i] - h; e4[j] <- e4[j] - h
      H[i, j] <- H[j, i] <- (f(e1) - f(e2) - f(e3) + f(e4)) / (4 * h * h)
    }
  }
  H
}

# This function estimates how several results change together.
.num_jac <- function(f, x, h = 1e-5) {
  f0 <- f(x)
  J  <- matrix(0, length(f0), length(x))
  for (j in seq_along(x)) {
    xp <- xm <- x
    xp[j] <- xp[j] + h
    xm[j] <- xm[j] - h
    J[, j] <- (f(xp) - f(xm)) / (2 * h)
  }
  dimnames(J) <- list(names(f0), names(x))
  J
}

# Joint covariance

# This function turns a deviance Hessian into a covariance matrix.
.hessian_vcov <- function(H, par) {
  if (is.null(H) || !identical(dim(H), c(length(par), length(par))) ||
      any(!is.finite(H))) {
    return(NULL)
  }

  info <- (H + t(H)) / 4
  values <- tryCatch(eigen(info, symmetric = TRUE, only.values = TRUE)$values,
                     error = function(e) numeric(0))
  scale <- max(abs(values), 1)
  if (!length(values) || min(values) <= sqrt(.Machine$double.eps) * scale ||
      !is.finite(rcond(info)) || rcond(info) < 1e-10) {
    return(NULL)
  }

  V <- tryCatch(solve(info), error = function(e) NULL)
  if (is.null(V) || any(!is.finite(V))) return(NULL)
  dimnames(V) <- list(names(par), names(par))
  V
}

# This function builds the full covariance matrix of a fitted model.
.full_vcov <- function(fit, devfun = NULL) {

  # The parameters follow the order expected by the model calculation.
  theta <- lme4::getME(fit, "theta")
  beta  <- lme4::fixef(fit)
  par   <- c(theta, beta)

  # Wald uncertainty is unsafe near a variance boundary.
  lower <- lme4::getME(fit, "lower")
  at_bound <- is.finite(lower) & theta - lower < 1e-3
  if (any(at_bound)) {
    .usdt_stop("the joint covariance is unavailable because a variance ",
               "component is on or near its boundary.")
  }

  # The fitted Hessian is the most accurate route.
  H <- tryCatch(fit@optinfo$derivs$Hessian, error = function(e) NULL)
  V <- .hessian_vcov(H, par)
  source <- "lme4"

  # A numerical Hessian covers unusual interior fits.
  if (is.null(V)) {
    if (is.null(devfun)) devfun <- .rebuild_devfun(fit)
    H <- tryCatch(.num_hess(devfun, par), error = function(e) NULL)
    V <- .hessian_vcov(H, par)
    source <- "numerical"
  }
  if (is.null(V)) {
    .usdt_stop("the observed information matrix is unavailable or unstable.")
  }

  attr(V, "par") <- par
  attr(V, "source") <- source
  V
}

# This function rebuilds a calculation that an external model did not store.
.rebuild_devfun <- function(fit) {

  # The function first asks the fitted model to rebuild the calculation.
  out <- tryCatch(stats::update(fit, devFunOnly = TRUE), error = function(e) NULL)
  if (is.function(out)) return(out)

  # The function refits the model when the first route fails.
  out <- tryCatch({
    cl  <- stats::getCall(fit)
    dat <- eval(cl$data, environment(stats::formula(fit)))
    lme4::glmer(stats::formula(fit), data = dat, family = stats::family(fit),
                control = lme4::glmerControl(optimizer = "bobyqa"),
                nAGQ = 1L, devFunOnly = TRUE)
  }, error = function(e) NULL)
  if (is.function(out)) return(out)

  .usdt_stop("could not rebuild the deviance function of this model, so the ",
             "covariance of the variance components is unavailable.\n",
             "  This happens when the data the model was fitted to can no ",
             "longer be found.\n",
             "  Either refit the model where its data are visible, or pass ",
             "the object returned by hsdt() instead of hsdt()$fit, ",
             "which carries the deviance function with it.")
}

# Random-effect variances and covariances

# This function keeps a correlation inside its range.
.cor_clamp <- function(rho) max(-1, min(1, rho))

# This function finds the random term that contains both sensitivities.
.locate_block <- function(fit, direct, indirect) {

  # The model lists the columns of each random term in parameter order.
  cnms  <- lme4::getME(fit, "cnms")
  sizes <- vapply(cnms, length, integer(1L))
  ends  <- cumsum(sizes * (sizes + 1L) / 2L)
  start <- c(0L, ends[-length(ends)])

  # The selected term must contain both sensitivities.
  hit <- which(vapply(cnms, function(z) all(c(direct, indirect) %in% z), TRUE))
  if (!length(hit)) {
    .usdt_stop("no random-effects term contains both `", direct, "` and `",
               indirect, "`.\n  Terms found: ",
               paste(vapply(cnms, function(z) paste0("(", paste(z, collapse = " + "),
                                                     " | ...)"), ""),
                     collapse = ", "), "\n",
               "  The latent correlation requires the two sensitivities to ",
               "vary together in the same term.")
  }
  if (length(hit) > 1L) hit <- hit[1L]

  # The result identifies the parameters for the selected term.
  q <- sizes[[hit]]
  list(idx  = seq.int(start[hit] + 1L, ends[hit]),
       q    = q,
       iD   = match(direct,   cnms[[hit]]),
       iI   = match(indirect, cnms[[hit]]),
       name = names(cnms)[hit])
}

# This function rebuilds one random-effects covariance matrix.
.theta_to_sigma <- function(theta_block, q) {
  L <- matrix(0, q, q)
  # The values fill the lower half of the matrix by columns.
  L[lower.tri(L, diag = TRUE)] <- theta_block
  tcrossprod(L)
}

# This function derives a covariance block from its fitted parameters.
#
#   d Sigma_ij / d L_kl = L_jl when i == k, plus L_il when j == k,
.sigma_gradient <- function(theta_block, q) {
  L  <- matrix(0, q, q)
  L[lower.tri(L, diag = TRUE)] <- theta_block
  at <- which(lower.tri(L, diag = TRUE), arr.ind = TRUE)
  lapply(seq_len(nrow(at)), function(m) {
    k <- at[m, 1L]
    column <- L[, at[m, 2L]]
    E <- matrix(0, q, q)
    E[k, ] <- column
    E[, k] <- E[, k] + column
    E
  })
}

# This function locates every random-effects block.
.theta_blocks <- function(fit) {
  cnms  <- lme4::getME(fit, "cnms")
  sizes <- vapply(cnms, length, integer(1L))
  ends  <- cumsum(sizes * (sizes + 1L) / 2L)
  Map(function(start, end, q) list(idx = seq.int(start + 1L, end), q = q),
      c(0L, ends[-length(ends)]), ends, sizes)
}

# This function estimates uncertainty for each random-effects deviation.
.re_summary <- function(fit, V = NULL, par) {

  blocks <- .theta_blocks(fit)
  labs   <- unlist(lme4::getME(fit, "cnms"), use.names = FALSE)

  est <- numeric(0)
  J   <- matrix(0, length(labs), length(par))
  row <- 0L
  for (b in blocks) {
    S <- .theta_to_sigma(par[b$idx], b$q)
    G <- .sigma_gradient(par[b$idx], b$q)
    for (i in seq_len(b$q)) {
      row <- row + 1L
      sd  <- sqrt(S[i, i])
      est <- c(est, sd)
      if (sd > 0) {
        J[row, b$idx] <- vapply(G, function(g) g[i, i], 0) / (2 * sd)
      }
    }
  }

  # The covariance gives the uncertainty of those deviations.
  if (is.null(V)) {
    se <- rep(NA_real_, length(est))
  } else {
    v  <- diag(J %*% V %*% t(J))
    se <- ifelse(is.finite(v) & v > 0, sqrt(v), NA_real_)
  }

  data.frame(term = labs, sd = est, se = unname(se),
             stringsAsFactors = FALSE)
}

# This function describes one correlation from a random-effects block.
.random_correlation <- function(fit, V, par, first, second) {
  blk <- .locate_block(fit, first, second)
  S   <- .theta_to_sigma(par[blk$idx], blk$q)
  a   <- S[blk$iD, blk$iD]; b <- S[blk$iI, blk$iI]; cc <- S[blk$iD, blk$iI]
  estimate <- cc / sqrt(a * b)

  se <- NA_real_
  if (!is.null(V) && is.finite(estimate) && a > 0 && b > 0) {
    G <- .sigma_gradient(par[blk$idx], blk$q)
    J <- matrix(0, 1L, length(par))
    J[1L, blk$idx] <- vapply(G, function(g) {
      g[blk$iD, blk$iI] / sqrt(a * b) -
        estimate * (g[blk$iD, blk$iD] / a + g[blk$iI, blk$iI] / b) / 2
    }, 0)
    variance <- drop(J %*% V %*% t(J))
    if (is.finite(variance) && variance >= 0) se <- sqrt(variance)
  }
  list(estimate = .cor_clamp(estimate), se = se)
}

# Conditional effects

# This function estimates uncertainty for each group's total effect.
#
#   Var(beta_k + b_j) = condVar_jj + t(g) V g,   g = d b_j / d phi + e_k
.conditional_se <- function(fit, V = NULL, devfun = NULL) {

  # lme4 supplies the uncertainty conditional on fitted parameters.
  re <- as.data.frame(lme4::ranef(fit, condVar = TRUE))

  # The joint covariance carries the remaining uncertainty.
  theta <- lme4::getME(fit, "theta")
  fixed <- lme4::fixef(fit)
  phi   <- c(theta, fixed)
  if (is.null(V)) {
    .usdt_stop("the joint covariance is required for subject intervals.")
  }

  # The model calculation shows how each mode changes.
  if (is.null(devfun)) devfun <- .rebuild_devfun(fit)
  modes <- function(p) {
    devfun(p)
    as.vector(environment(devfun)$pp$b(1))
  }
  D <- .num_jac(modes, phi)

  # Each reported effect is matched to the model vector.
  cnms   <- lme4::getME(fit, "cnms")
  flist  <- lme4::getME(fit, "flist")
  sizes  <- vapply(cnms, length, 1L)
  counts <- vapply(names(cnms), function(g) nlevels(flist[[g]]), 1L)
  starts <- c(0L, cumsum(sizes * counts))

  index <- integer(nrow(re))
  for (k in seq_along(cnms)) {
    levels_k <- levels(flist[[names(cnms)[k]]])
    for (j in seq_along(cnms[[k]])) {
      rows <- which(re$grpvar == names(cnms)[k] & re$term == cnms[[k]][j])
      position <- match(as.character(re$grp[rows]), levels_k)
      index[rows] <- starts[k] + (position - 1L) * sizes[k] + j
    }
  }

  # The fixed effect completes each subject estimate.
  at <- length(theta) + match(as.character(re$term), names(fixed))
  if (any(index < 1L) || anyNA(at)) {
    .usdt_stop("the fitted fixed and random terms could not be matched.")
  }
  re$se <- vapply(seq_len(nrow(re)), function(i) {
    g <- D[index[i], ]
    g[at[i]] <- g[at[i]] + 1
    variance <- re$condsd[i]^2 + drop(crossprod(g, V %*% g))
    if (!is.finite(variance) || variance < 0) {
      .usdt_stop("the subject variance is unavailable or unstable.")
    }
    sqrt(variance)
  }, 0)
  re
}

# This function calculates the five values used by the hypotheses.
.quantity_parts <- function(fit, direct, indirect, par) {
  blk <- .locate_block(fit, direct, indirect)
  nt  <- length(lme4::getME(fit, "theta"))
  fx  <- names(lme4::fixef(fit))
  S   <- .theta_to_sigma(par[blk$idx], blk$q)
  G   <- .sigma_gradient(par[blk$idx], blk$q)

  est <- c(gamma_D = unname(par[nt + match(direct, fx)]),
           gamma_I = unname(par[nt + match(indirect, fx)]),
           s2_D    = S[blk$iD, blk$iD],
           s2_I    = S[blk$iI, blk$iI],
           s_DI    = S[blk$iD, blk$iI])

  J <- matrix(0, 5L, length(par), dimnames = list(names(est), names(par)))
  J[1L, nt + match(direct, fx)]   <- 1
  J[2L, nt + match(indirect, fx)] <- 1
  J[3L, blk$idx] <- vapply(G, function(g) g[blk$iD, blk$iD], 0)
  J[4L, blk$idx] <- vapply(G, function(g) g[blk$iI, blk$iI], 0)
  J[5L, blk$idx] <- vapply(G, function(g) g[blk$iD, blk$iI], 0)

  list(est = est, jacobian = J)
}

# This function collects the estimates and their joint covariance.
.usdt_pars <- function(fit, direct, indirect, devfun = NULL) {

  # The model must contain both sensitivities as fixed effects.
  fx <- lme4::fixef(fit)
  missing <- setdiff(c(direct, indirect), names(fx))
  if (length(missing)) {
    .usdt_stop("fixed effect", if (length(missing) > 1L) "s" else "", " ",
               paste0("`", missing, "`", collapse = " and "),
               " not found in the model.\n  Available: ",
               paste(names(fx), collapse = ", "))
  }

  # The fitted model always provides the five point estimates.
  par   <- c(lme4::getME(fit, "theta"), fx)
  parts <- .quantity_parts(fit, direct, indirect, par)
  est   <- parts$est

  # The fixed-effects covariance remains available without the full matrix.
  V_fixed <- tryCatch(as.matrix(stats::vcov(fit)), error = function(e) NULL)

  # The full covariance can fail without losing the fitted model.
  full_error <- NULL
  V <- tryCatch(.full_vcov(fit, devfun = devfun),
                error = function(e) {
                  full_error <<- conditionMessage(e)
                  NULL
                })
  Vq <- matrix(NA_real_, length(est), length(est),
               dimnames = list(names(est), names(est)))
  if (!is.null(V)) {
    Vq <- parts$jacobian %*% V %*% t(parts$jacobian)
    dimnames(Vq) <- list(names(est), names(est))
  }

  # The sensitivity block must be inside its valid parameter range.
  rho <- .cor_clamp(est[["s_DI"]] / sqrt(est[["s2_D"]] * est[["s2_I"]]))
  S <- matrix(c(est[["s2_D"]], est[["s_DI"]],
                est[["s_DI"]], est[["s2_I"]]), 2L)
  s_values <- tryCatch(eigen(S, symmetric = TRUE, only.values = TRUE)$values,
                       error = function(e) NA_real_)
  sensitivity_boundary <- any(!is.finite(c(est[c("s2_D", "s2_I", "s_DI")],
                                               rho, s_values))) ||
    min(est[["s2_D"]], est[["s2_I"]]) < 1e-8 ||
    min(s_values) < 1e-8 || abs(rho) > 1 - 1e-6

  lower <- lme4::getME(fit, "lower")
  theta <- lme4::getME(fit, "theta")
  at_bound <- names(theta)[theta - lower < 1e-3]
  joint_ok <- !sensitivity_boundary && !is.null(V)
  reason <- if (sensitivity_boundary) {
    "the sensitivity covariance is on the boundary"
  } else if (is.null(V)) {
    full_error %||% "the full covariance matrix is unavailable"
  } else {
    NA_character_
  }

  criterion_cor <- if (all(c("c_D", "c_I") %in% names(fx))) {
    .random_correlation(fit, V, par, "c_D", "c_I")
  } else {
    NULL
  }

  list(est = est, vcov = Vq, fixed_vcov = V_fixed,
       fixed_names = c(direct, indirect), rho = rho,
       V_full = V, par = par, re = .re_summary(fit, V, par),
       criterion_cor = criterion_cor,
       at_bound = at_bound,
       singular = isTRUE(lme4::isSingular(fit, tol = 1e-4)),
       sensitivity_boundary = sensitivity_boundary,
       joint_ok = joint_ok, inference_reason = reason,
       boundary = !joint_ok)
}
