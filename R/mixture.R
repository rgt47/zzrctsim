# Finite mixtures of multivariate normals.
#
# A component label is drawn once per subject and held for all that
# subject's visits, so a mixture describes *latent subpopulations* --
# fast and slow progressors, say -- rather than visit-level noise. That
# is the structure repo 29 needs, and it is the one that makes the
# mixture visible in the between-subject covariance rather than
# averaging away within it.
#
# Because every component has finite moments, the mixture has finite
# moments too, and the covariance machinery of the Gaussian case
# continues to apply: `cov_at()` returns the marginal covariance of the
# mixture, which is the usual within-plus-between decomposition
#
#     V = sum_k w_k V_k  +  sum_k w_k (mu_k - mubar)(mu_k - mubar)'
#
# and is exact, not an approximation. This is what distinguishes a
# normal mixture from the stable mixtures of the requirements document,
# which have no variance to decompose.

#' Specify a finite mixture of multivariate normals
#'
#' Each subject is assigned to a component with probability `weights`,
#' and their whole response vector is drawn from that component. The
#' assignment is per subject, not per visit.
#'
#' @param components A list of `dgm` objects, one per mixture
#'   component. Each supplies that component's covariance through
#'   [cov_at()], so components may differ in covariance structure, not
#'   only in mean.
#' @param weights Numeric vector of mixing probabilities, summing to
#'   one. Defaults to equal weights.
#' @param mean_shift Numeric vector of component-specific additive
#'   shifts to the mean, one per component, or a list of functions of
#'   the time vector for time-varying shifts. Defaults to zero, in
#'   which case the components differ only in covariance.
#' @return An object of class `dgm_mixture`.
#' @details
#' The `mean_shift` is added to the fixed-effect mean supplied to
#' [generate_outcomes()]. A shift that is constant in time moves the
#' whole trajectory; a shift given as a function of time can express
#' components that differ in *slope*, which is the usual way to encode
#' fast and slow progressors.
#'
#' Identifiability is the user's problem, not the package's. Two
#' components with equal weights, equal covariances and mean shifts
#' differing by less than a residual standard deviation are a mixture
#' in name only, and no analysis will recover them. Simulating such a
#' case is legitimate -- demonstrating that it cannot be recovered is a
#' finding -- but it should be deliberate.
#' @examples
#' slow <- dgm_conditional(G = diag(c(4, 0.01)), sigma2 = 2)
#' fast <- dgm_conditional(G = diag(c(9, 0.04)), sigma2 = 2)
#' mix <- dgm_mixture(list(slow, fast), weights = c(0.7, 0.3),
#'                    mean_shift = list(function(t) 0 * t,
#'                                      function(t) 0.4 * t))
#' mix
#' @export
dgm_mixture <- function(components, weights = NULL,
                        mean_shift = NULL) {
  stopifnot(is.list(components), length(components) >= 2)
  if (!all(vapply(components, inherits, logical(1), "dgm"))) {
    stop("every element of `components` must be a `dgm` object.")
  }
  if (any(vapply(components, inherits, logical(1), "dgm_tte"))) {
    stop("`dgm_mixture()` mixes response distributions; a `dgm_tte` ",
         "component has no `cov_at()` and cannot be mixed this way.")
  }
  K <- length(components)

  if (is.null(weights)) weights <- rep(1 / K, K)
  stopifnot(length(weights) == K, all(weights > 0))
  if (abs(sum(weights) - 1) > 1e-8) {
    stop("`weights` must sum to 1; they sum to ",
         format(sum(weights), digits = 8), ".")
  }

  if (is.null(mean_shift)) {
    mean_shift <- rep(list(function(t) rep(0, length(t))), K)
  } else if (is.numeric(mean_shift)) {
    stopifnot(length(mean_shift) == K)
    consts <- as.list(mean_shift)
    mean_shift <- lapply(consts, function(cst) {
      force(cst)
      function(t) rep(cst, length(t))
    })
  } else {
    stopifnot(is.list(mean_shift), length(mean_shift) == K,
              all(vapply(mean_shift, is.function, logical(1))))
  }

  structure(list(components = components, weights = weights,
                 mean_shift = mean_shift, K = K),
            class = c("dgm_mixture", "dgm"))
}

#' @export
print.dgm_mixture <- function(x, ...) {
  cat("<dgm_mixture>", x$K, "components\n")
  cat("  weights:", paste(format(x$weights, digits = 3),
                          collapse = "  "), "\n")
  shifts <- vapply(x$mean_shift, function(f) f(1)[1], numeric(1))
  cat("  mean shift at t = 1:",
      paste(format(shifts, digits = 3), collapse = "  "), "\n")
  invisible(x)
}

#' @export
cov_at.dgm_mixture <- function(dgm, times) {
  # Marginal covariance of the mixture: within-component covariance
  # averaged over components, plus the between-component spread of the
  # mean shifts. Both terms are exact.
  w <- dgm$weights
  Vs <- lapply(dgm$components, cov_at, times = times)
  shifts <- lapply(dgm$mean_shift, function(f) f(times))
  mbar <- Reduce(`+`, Map(function(wi, si) wi * si, w, shifts))

  within <- Reduce(`+`, Map(function(wi, Vi) wi * Vi, w, Vs))
  between <- Reduce(`+`, Map(function(wi, si) {
    dv <- si - mbar
    wi * outer(dv, dv)
  }, w, shifts))
  within + between
}

# Mixture generation: assign a component per subject, then draw that
# subject's whole vector from it. Subjects sharing a component and a
# time vector are drawn together.
.generate_mixture <- function(design, dgm) {
  ids <- unique(design$id)
  n <- length(ids)
  comp <- sample.int(dgm$K, size = n, replace = TRUE,
                     prob = dgm$weights)
  names(comp) <- as.character(ids)

  design$component <- comp[as.character(design$id)]
  shift <- numeric(nrow(design))
  for (k in seq_len(dgm$K)) {
    sel <- design$component == k
    if (any(sel)) {
      shift[sel] <- dgm$mean_shift[[k]](design$time[sel])
    }
  }
  design$mu <- design$mu + shift

  y <- numeric(nrow(design))
  for (k in seq_len(dgm$K)) {
    rows <- which(design$component == k)
    if (!length(rows)) next
    y[rows] <- .draw_gaussian(design[rows, , drop = FALSE],
                              dgm$components[[k]])
  }
  design$y <- y
  attr(design, "component") <- comp
  attr(design, "family") <- "gaussian"
  design
}
