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
#'   by [runin_design()] -- either from a balanced [trial_schedule()],
#'   or from a ragged [realize_schedule()] under staggered accrual,
#'   which [runin_design()] also accepts. Must contain `id` and
#'   `time`, plus any columns named in `beta`.
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
