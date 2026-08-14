# Generation: draw outcomes for a realized schedule from a DGM.
#
# Complete data are always generated first; missingness is applied as a
# separate layer afterwards. This is what permits the information loss
# from dropout to be isolated by comparison against the complete-data
# run, and it is what makes response-dependent (MAR, MNAR) mechanisms
# expressible at all.
#
# Subjects may be observed on different grids, because staggered
# accrual with a common close-out gives early enrollees more visits.
# Each subject's covariance is therefore evaluated at that subject's
# own times, which is why `cov_at()` is a function of time rather than
# a fixed matrix.

#' Generate outcomes for a trial
#'
#' @param design A data frame of one row per subject-visit, as returned
#'   by [runin_design()] or by [realize_schedule()] passed through
#'   [runin_design()]. Must contain `id` and `time`, plus any columns
#'   named in `beta`.
#' @param dgm A `dgm` object supplying the within-subject covariance.
#' @param beta Named numeric vector of fixed effects. Names must match
#'   columns of `design`; an `"(Intercept)"` entry is added implicitly
#'   if absent from `design`.
#' @param intercept Numeric. Overall mean at time zero.
#' @return `design` with `mu`, the fixed-effect mean, and `y`, the
#'   generated outcome, added. Every path additionally attaches the
#'   attributes `beta`, `intercept`, and `family`, which
#'   [reference_based()] requires.
#'
#'   Two paths add more. A non-Gaussian `dgm_conditional` (the GLMM
#'   construction) also adds the columns `eta`, the linear predictor,
#'   and `cmean`, the conditional mean on the response scale, plus a
#'   `ranef` attribute holding the drawn random effects. A
#'   `dgm_mixture` also adds the column `component`, the latent class
#'   each subject was drawn from, and an attribute of the same name.
#' @details
#' Draws are made per subject from a multivariate normal with mean
#' `mu` and covariance `cov_at(dgm, times_i)`. Subjects with a single
#' observation are handled as a univariate draw.
#' @examples
#' s <- trial_schedule(run_in = 1, treatment = 4, interval = 3)
#' arm <- factor(rep(c("placebo", "active"), each = 5),
#'               levels = c("placebo", "active"))
#' d <- runin_design(s, arm, reference = "placebo")
#' g <- dgm_conditional(G = diag(c(4, 0.01)), sigma2 = 2)
#' out <- generate_outcomes(d, g,
#'   beta = c(x_slope = 0.5, x_hinge = -0.1, x_trt_active = -0.2))
#' head(out)
#' @export
generate_outcomes <- function(design, dgm, beta, intercept = 0) {
  stopifnot(is.data.frame(design), inherits(dgm, "dgm"))
  if (is.null(names(beta)) || any(names(beta) == "")) {
    stop("`beta` must be fully named.")
  }
  miss <- setdiff(names(beta), names(design))
  if (length(miss)) {
    stop("`beta` names not found in `design`: ",
         paste(miss, collapse = ", "))
  }

  X <- as.matrix(design[, names(beta), drop = FALSE])
  design$mu <- intercept + as.vector(X %*% beta)

  # Non-Gaussian families cannot be drawn from an induced marginal,
  # because no multivariate analogue of the binomial or Poisson exists.
  # The random effects are therefore drawn explicitly and the response
  # conditionally, which is the GLMM construction.
  out <- if (inherits(dgm, "dgm_conditional") && !isTRUE(dgm$gaussian)) {
    .generate_glmm(design, dgm)
  } else if (inherits(dgm, "dgm_mixture")) {
    .generate_mixture(design, dgm)
  } else {
    design$y <- .draw_gaussian(design, dgm)
    attr(design, "family") <- "gaussian"
    design
  }

  # `beta` and `intercept` are attached here rather than inside each
  # branch, so that every generated data set carries them regardless of
  # which generator produced it. `reference_based()` requires both, and
  # attaching them on the Gaussian path alone silently restricted that
  # function to Gaussian data. `family` is set by each branch, since
  # only the branch knows it.
  attr(out, "beta") <- beta
  attr(out, "intercept") <- intercept
  out
}


# Draw a Gaussian response for a design that already carries `mu`.
#
# Subjects sharing a time vector share a covariance and are drawn in
# one call, so a balanced design collapses to a single `mvrnorm()`;
# only genuinely ragged grids, as staggered accrual produces, need more
# than one. Shared by the plain Gaussian path and by each component of
# a mixture.
.draw_gaussian <- function(design, dgm) {
  ord <- order(design$id, design$time)
  sig <- vapply(split(design$time[ord], design$id[ord]),
                function(t) paste0(format(t, digits = 12),
                                   collapse = ","),
                character(1))
  ids <- names(sig)
  rows_by_id <- split(ord, design$id[ord])

  y <- numeric(nrow(design))
  for (s in unique(sig)) {
    grp <- ids[sig == s]
    k_list <- rows_by_id[grp]
    t_i <- design$time[k_list[[1]]]
    V_i <- cov_at(dgm, t_i)
    m <- length(grp)
    p <- length(t_i)
    draws <- if (p == 1L) {
      matrix(stats::rnorm(m, 0, sqrt(V_i[1, 1])), nrow = m)
    } else {
      matrix(MASS::mvrnorm(m, mu = rep(0, p), Sigma = V_i),
             nrow = m, ncol = p)
    }
    idx <- unlist(k_list, use.names = FALSE)
    y[idx] <- design$mu[idx] + as.vector(t(draws))
  }
  y
}

# GLMM path: draw b_i ~ N(0, G), form the linear predictor, apply the
# inverse link, then draw the response conditionally. Adds `eta` (the
# full linear predictor) and `cmean` (the conditional mean) alongside
# `mu` (the fixed-effect part), so all three scales are inspectable.
.generate_glmm <- function(design, dgm) {
  ids <- unique(design$id)
  n <- length(ids)
  q <- dgm$q

  b <- if (q == 1L) {
    matrix(stats::rnorm(n, 0, sqrt(dgm$G[1, 1])), ncol = 1L)
  } else {
    matrix(MASS::mvrnorm(n, mu = rep(0, q), Sigma = dgm$G),
           nrow = n, ncol = q)
  }
  rownames(b) <- as.character(ids)

  eta <- numeric(nrow(design))
  for (i in ids) {
    k <- which(design$id == i)
    Z <- matrix(dgm$z(design$time[k]), nrow = length(k))
    if (ncol(Z) != q) {
      stop("`z()` returned ", ncol(Z), " columns but `G` is ",
           q, " x ", q, ".")
    }
    eta[k] <- design$mu[k] + as.vector(Z %*% b[as.character(i), ])
  }

  cmean <- dgm$family$linkinv(eta)
  design$eta <- eta
  design$cmean <- cmean
  design$y <- .rfamily(dgm, cmean)
  attr(design, "ranef") <- b
  attr(design, "family") <- dgm$family$family
  design
}

#' Name the missingness mechanism implied by a hazard specification
#'
#' The three mechanisms are nested restrictions of one model, not three
#' separate models.
#'
#' @param psi1 Coefficient on the previously observed response.
#' @param psi2 Coefficient on the current, unobserved response.
#' @return `"MCAR"`, `"MAR"`, or `"MNAR"`.
#' @export
dropout_mechanism <- function(psi1 = 0, psi2 = 0) {
  if (psi2 != 0) return("MNAR")
  if (psi1 != 0) return("MAR")
  "MCAR"
}

# Per-visit dropout hazards for every subject-visit eligible for
# dropout, given a hazard specification. Returns a list of numeric
# vectors, one per subject, aligned to that subject's eligible visits.
.hazards <- function(dat, psi0, psi1, psi2, center, from) {
  ids <- unique(dat$id)
  lapply(ids, function(i) {
    k <- which(dat$id == i)
    k <- k[order(dat$time[k])]
    yv <- dat$y[k]
    tv <- dat$time[k]
    elig <- which(tv > from)
    if (!length(elig)) return(numeric(0))
    # No predecessor means no previously observed response to condition
    # on, so the history term contributes nothing. Using the current
    # response here instead would be MNAR while still reporting MAR.
    # See the matching note in `.haz_rows()` in R/missing.R.
    prev <- ifelse(elig > 1L, yv[pmax(elig - 1L, 1L)], center)
    stats::plogis(psi0 + psi1 * (prev - center) +
                    psi2 * (yv[elig] - center))
  })
}

# Expected proportion of subjects ever dropping out, given complete
# data. Exact, because the hazards depend only on complete-data
# responses: P(ever drop) = 1 - prod(1 - p_j). Smooth and strictly
# increasing in psi0, so uniroot is well behaved.
.expected_dropout <- function(dat, psi0, psi1, psi2, center, from) {
  h <- .hazards(dat, psi0, psi1, psi2, center, from)
  mean(vapply(h, function(p) {
    if (!length(p)) 0 else 1 - prod(1 - p)
  }, numeric(1)))
}

#' Apply monotone dropout to generated data
#'
#' Missingness is applied after complete data are generated, following
#' the selection-model factorization `f(Y, R) = f(Y) f(R | Y)`, so the
#' dropout hazard may depend on realized responses.
#'
#' @param dat Output of [generate_outcomes()].
#' @param target Numeric in `(0, 1)`. Target proportion of subjects who
#'   ever drop out. `psi0` is calibrated to achieve it. Ignored if
#'   `psi0` is supplied.
#' @param psi0 Numeric. Intercept of the dropout logit. Supply directly
#'   to bypass calibration.
#' @param psi1 Numeric. Coefficient on the previously observed
#'   response, centred. Non-zero gives MAR.
#' @param psi2 Numeric. Coefficient on the current, unobserved
#'   response, centred. Non-zero gives MNAR.
#' @param center Numeric. Value at which the response is centred in the
#'   hazard. Defaults to the mean of the complete-data response, which
#'   keeps `psi0` interpretable as a baseline log-odds.
#' @param from Numeric. Dropout is only permitted at visits with
#'   `time > from`, so run-in and baseline visits are retained.
#' @return `dat` with `y` set to `NA` at and after each subject's
#'   dropout visit and a `dropped` logical column added. Attributes
#'   `dropout_spec` (the realized `psi0`, `psi1`, `psi2`, `center`,
#'   mechanism, target and expected proportion) are attached.
#' @details
#' The hazard follows Diggle and Kenward (1994) \doi{10.2307/2986113}:
#' for a subject still under observation at visit `j`,
#' \deqn{\mathrm{logit}\, P(\mathrm{drop\ at\ } j) =
#'   \psi_0 + \psi_1 (y_{j-1} - c) + \psi_2 (y_j - c).}
#' MCAR, MAR and MNAR are the nested restrictions
#' `psi1 = psi2 = 0`, `psi2 = 0`, and `psi2` free. MNAR therefore
#' *extends* MAR rather than replacing its history dependence.
#'
#' Centring at `c` matters: on an untransformed ADAS-Cog or CDR-SB
#' scale an uncentred response saturates the logit, so that any
#' non-trivial `psi1` drives the dropout proportion to one regardless
#' of the intended level.
#'
#' Calibration solves for `psi0` such that the expected proportion of
#' subjects ever dropping equals `target`. Because the hazards depend
#' only on complete-data responses, that expectation is available in
#' closed form as `1 - prod(1 - p_j)` averaged over subjects, and is
#' strictly increasing in `psi0`.
#' @examples
#' s <- trial_schedule(run_in = 1, treatment = 5, interval = 3)
#' arm <- factor(rep(c("placebo", "active"), each = 50),
#'               levels = c("placebo", "active"))
#' d <- runin_design(s, arm, reference = "placebo")
#' g <- dgm_conditional(G = diag(c(9, 0.04)), sigma2 = 4)
#' out <- generate_outcomes(d, g,
#'   beta = c(x_slope = 0.5, x_hinge = -0.1, x_trt_active = -0.25),
#'   intercept = 20)
#' mar <- apply_dropout(out, target = 0.25, psi1 = 0.15)
#' attr(mar, "dropout_spec")$mechanism
#' @export
apply_dropout <- function(dat,
                          target = 0.25,
                          psi0 = NULL,
                          psi1 = 0,
                          psi2 = 0,
                          center = NULL,
                          from = 0) {
  stopifnot(is.data.frame(dat), "y" %in% names(dat))
  .check_family_for_hazard(dat, psi1, psi2, center)
  if (is.null(center)) center <- mean(dat$y, na.rm = TRUE)

  # Recorded before `psi0` is overwritten by calibration below, so the
  # spec can say whether `target` was actually used. Supplying `psi0`
  # bypasses calibration, and reporting the default `target` in that
  # case would misdescribe the mechanism that was applied.
  calibrated <- is.null(psi0)

  if (is.null(psi0)) {
    stopifnot(target > 0, target < 1)
    f <- function(p) {
      .expected_dropout(dat, p, psi1, psi2, center, from) - target
    }
    lo <- -50; hi <- 50
    if (f(lo) > 0 || f(hi) < 0) {
      stop("Cannot calibrate psi0 to a target of ", target,
           " with psi1 = ", psi1, ", psi2 = ", psi2,
           ". The attainable range is [",
           format(f(lo) + target, digits = 3), ", ",
           format(f(hi) + target, digits = 3), "].")
    }
    psi0 <- stats::uniroot(f, lower = lo, upper = hi,
                           tol = 1e-8)$root
  }

  expected <- .expected_dropout(dat, psi0, psi1, psi2, center, from)
  spec <- list(psi0 = psi0, psi1 = psi1, psi2 = psi2,
               center = center,
               mechanism = dropout_mechanism(psi1, psi2),
               target = if (calibrated) target else NA_real_,
               expected = expected)

  dat$dropped <- FALSE
  haz <- .hazards(dat, psi0, psi1, psi2, center, from)
  ids <- unique(dat$id)

  for (a in seq_along(ids)) {
    i <- ids[a]
    k <- which(dat$id == i)
    k <- k[order(dat$time[k])]
    tv <- dat$time[k]
    elig <- which(tv > from)
    p <- haz[[a]]
    if (!length(p)) next
    u <- stats::runif(length(p))
    hit <- which(u < p)
    if (!length(hit)) next
    m <- elig[hit[1]]
    drop_idx <- k[m:length(k)]
    dat$y[drop_idx] <- NA_real_
    dat$dropped[drop_idx] <- TRUE
  }

  attr(dat, "dropout_spec") <- spec
  dat
}
