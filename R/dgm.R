# Data-generating mechanisms in either parameterization.
#
# The same longitudinal DGM can be written two ways:
#
#   conditional  (Z, G, sigma^2)  random effects plus i.i.d. residuals
#   marginal     V                the response covariance directly
#
# Forward conversion is exact: V = Z G Z' + sigma^2 I. Backward is not
# always possible, and `certify()` settles when it is.
#
# Because staggered accrual makes the observed grid subject specific,
# a marginal specification must be a covariance *function* of the
# observation times, not a fixed matrix, unless the design is balanced.

#' Continuous-time covariance structures
#'
#' Return a function of a time vector giving the covariance at those
#' times. Because they are defined on continuous time rather than on a
#' visit index, they extend to unequally spaced and subject-specific
#' grids, which is what staggered accrual requires.
#'
#' @param sd Numeric scalar, numeric vector matching the reference
#'   grid, or a function of time returning standard deviations.
#' @param rho Correlation parameter.
#' @return A function of a numeric time vector returning a covariance
#'   matrix.
#' @name cov_structures
NULL

.sd_at <- function(sd, t) {
  if (is.function(sd)) return(sd(t))
  # A scalar is broadcast, but a vector must match the grid exactly.
  # Partial recycling is refused: these covariances are evaluated on
  # subject-specific grids of varying length under staggered accrual,
  # where a recycled `sd` silently assigns the wrong variance to the
  # wrong visit.
  if (length(sd) != 1L && length(sd) != length(t)) {
    stop("`sd` has length ", length(sd), " but the time grid has ",
         length(t), " point(s); supply a scalar, a vector of matching ",
         "length, or a function of time.")
  }
  rep(sd, length.out = length(t))
}

#' @rdname cov_structures
#' @details `cov_ar1()` uses `rho^|t_i - t_j|`, the continuous-time
#'   AR(1) or exponential correlation, which reduces to the usual
#'   lag-based form on an equally spaced grid.
#' @export
cov_ar1 <- function(sd = 1, rho = 0.6) {
  stopifnot(rho > 0, rho < 1)
  function(t) {
    s <- .sd_at(sd, t)
    diag(s, nrow = length(s)) %*%
      rho^abs(outer(t, t, "-")) %*%
      diag(s, nrow = length(s))
  }
}

#' @rdname cov_structures
#' @export
cov_cs <- function(sd = 1, rho = 0.5) {
  stopifnot(rho > -1, rho < 1)
  function(t) {
    s <- .sd_at(sd, t)
    R <- matrix(rho, length(t), length(t))
    diag(R) <- 1
    diag(s, nrow = length(s)) %*% R %*% diag(s, nrow = length(s))
  }
}

#' Specify a data-generating mechanism
#'
#' `dgm_conditional()` declares random effects; `dgm_marginal()`
#' declares the response covariance directly. Both produce a `dgm`
#' object that answers [cov_at()].
#'
#' @param G Random-effects covariance, a `q x q` matrix (or scalar for
#'   a single random effect).
#' @param sigma2 Residual variance. Required for the Gaussian family
#'   and ignored otherwise, since the conditional variance of the other
#'   families is determined by their mean.
#' @param z Function of a time vector returning the `n x q`
#'   random-effects design. Defaults to a random intercept and slope,
#'   `cbind(1, t)`.
#' @param family A [stats::family()] object giving the conditional
#'   distribution of the response given the random effects, and the
#'   link. Supported: `gaussian`, `binomial`, `poisson`, and
#'   `MASS::negative.binomial`.
#' @param trials Integer. Number of Bernoulli trials per observation
#'   for the binomial family. `1` gives a binary outcome.
#' @param theta Numeric. Dispersion for the negative binomial, in the
#'   `size` parameterization of [stats::rnbinom()].
#' @return An object of class `dgm`. The two constructors return
#'   different element sets, since they hold different
#'   parameterizations.
#'
#'   `dgm_conditional()` returns a list of class
#'   `c("dgm_conditional", "dgm")` with
#'   \describe{
#'     \item{G}{the `q x q` random-effects covariance, coerced to a
#'       matrix}
#'     \item{sigma2}{the residual variance; `NA_real_` for a
#'       non-Gaussian family, whose conditional variance is set by its
#'       mean}
#'     \item{z}{the random-effects design function of time, as supplied}
#'     \item{q}{integer, `nrow(G)`, the number of random effects}
#'     \item{family}{the [stats::family()] object, after coercion from
#'       a name or a generator}
#'     \item{trials}{number of Bernoulli trials per observation, used
#'       only by the binomial family}
#'     \item{theta}{negative-binomial dispersion, or `NULL`}
#'     \item{gaussian}{logical, whether the family is Gaussian; the
#'       flag [cov_at()] and the generator branch on}
#'   }
#'
#'   `dgm_marginal()` returns a list of class
#'   `c("dgm_marginal", "dgm")` with either
#'   \describe{
#'     \item{Vfun}{the covariance function of time, when `V` was given
#'       as a function, with `times` set to `NULL` so the object is
#'       evaluable on any grid}
#'     \item{V}{the covariance matrix, when `V` was given as a matrix}
#'     \item{times}{the reference grid `V` refers to, in the matrix
#'       case; [cov_at()] then refuses times off this grid}
#'   }
#' @details
#' With a non-Gaussian family the response is drawn conditionally on
#' the random effects: `b_i ~ N(0, G)`, then
#' `y_ij ~ family(g^-1(eta_ij))` with
#' `eta_ij = x_ij'beta + z_ij'b_i`. Two consequences follow and are
#' the user's to manage rather than the package's to hide.
#'
#' First, dependence between visits arises solely through the shared
#' `b_i`, so the reachable dependence is exactly what the
#' random-effects structure can induce. There is no marginal escape
#' hatch: [dgm_marginal()] has no meaning outside the Gaussian case.
#'
#' Second, under a non-identity link the marginal mean is not
#' `g^-1(x'beta)`. The fixed effects are conditional on the random
#' effects, so a conditional and a marginal treatment effect differ in
#' magnitude. An [estimand()] used with a non-Gaussian family must
#' state which of the two it refers to.
#' @export
dgm_conditional <- function(G, sigma2 = NULL,
                            z = function(t) cbind(1, t),
                            family = stats::gaussian(),
                            trials = 1L, theta = NULL) {
  G <- as.matrix(G)
  stopifnot(nrow(G) == ncol(G), is.function(z))
  if (is.character(family)) family <- get(family, mode = "function")()
  if (is.function(family)) family <- family()
  if (!inherits(family, "family")) {
    stop("`family` must be a family object, e.g. binomial('logit').")
  }
  if (!isSymmetric(unname(G), tol = sqrt(.Machine$double.eps))) {
    stop("`G` must be symmetric.")
  }
  if (min(eigen(G, symmetric = TRUE, only.values = TRUE)$values) <
      -1e-9) {
    stop("`G` must be positive semi-definite.")
  }

  fam <- family$family
  gaussian_fam <- identical(fam, "gaussian")
  if (gaussian_fam) {
    if (is.null(sigma2)) {
      stop("`sigma2` is required for the Gaussian family.")
    }
    stopifnot(sigma2 >= 0)
  } else if (is.null(sigma2)) {
    sigma2 <- NA_real_
  }
  if (grepl("^Negative Binomial", fam) && is.null(theta)) {
    stop("`theta` is required for the negative binomial family.")
  }
  known <- gaussian_fam || fam %in% c("binomial", "poisson") ||
    grepl("^Negative Binomial", fam)
  if (!known) {
    stop("Unsupported family '", fam, "'. Supported: gaussian, ",
         "binomial, poisson, negative binomial.")
  }

  structure(list(G = G, sigma2 = sigma2, z = z, q = nrow(G),
                 family = family, trials = trials, theta = theta,
                 gaussian = gaussian_fam),
            class = c("dgm_conditional", "dgm"))
}

# Draw from the conditional distribution given a mean vector.
.rfamily <- function(dgm, mu) {
  n <- length(mu)
  fam <- dgm$family$family
  if (dgm$gaussian) {
    return(stats::rnorm(n, mean = mu, sd = sqrt(dgm$sigma2)))
  }
  if (fam == "binomial") {
    return(stats::rbinom(n, size = dgm$trials, prob = mu) /
             if (dgm$trials == 1L) 1L else dgm$trials)
  }
  if (fam == "poisson") {
    return(stats::rpois(n, lambda = mu))
  }
  stats::rnbinom(n, mu = mu, size = dgm$theta)
}

#' @rdname dgm_conditional
#' @param V Either a covariance matrix, in which case `times` must be
#'   supplied to say which grid it refers to, or a function of a time
#'   vector returning a covariance matrix.
#' @param times Numeric vector. The reference grid, required when `V`
#'   is a matrix.
#' @export
dgm_marginal <- function(V, times = NULL) {
  if (is.function(V)) {
    return(structure(list(Vfun = V, times = NULL),
                     class = c("dgm_marginal", "dgm")))
  }
  V <- as.matrix(V)
  if (is.null(times)) {
    stop("`times` is required when `V` is a matrix; it identifies ",
         "the grid the matrix refers to.")
  }
  stopifnot(nrow(V) == length(times))
  if (!isSymmetric(unname(V), tol = sqrt(.Machine$double.eps))) {
    stop("`V` must be symmetric.")
  }
  # Checked here rather than left to the draw: a non-PSD `V` otherwise
  # fails inside MASS::mvrnorm once per replicate, deep in the call
  # stack and with no mention of `V` or of this constructor.
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  tol <- sqrt(.Machine$double.eps) * max(1, abs(ev[1]))
  if (any(ev < -tol)) {
    stop("`V` must be positive semi-definite; its smallest ",
         "eigenvalue is ", format(min(ev), digits = 4), ".")
  }
  structure(list(V = V, times = times),
            class = c("dgm_marginal", "dgm"))
}

#' Covariance of a data-generating mechanism at given times
#'
#' The single polymorphic operation both parameterizations answer.
#'
#' @param dgm A `dgm` object.
#' @param times Numeric vector of observation times.
#' @return A covariance matrix of dimension `length(times)`.
#' @export
cov_at <- function(dgm, times) UseMethod("cov_at")

#' @export
cov_at.dgm_conditional <- function(dgm, times) {
  if (!isTRUE(dgm$gaussian)) {
    stop("`cov_at()` is the covariance of the response and is defined ",
         "only for the Gaussian family. For '", dgm$family$family,
         "' the response variance is determined by its mean, so it ",
         "cannot be specified separately. Use `linpred_cov()` for the ",
         "covariance of the linear predictor.")
  }
  .zgz(dgm, times) + dgm$sigma2 * diag(length(times))
}

.zgz <- function(dgm, times) {
  Z <- matrix(dgm$z(times), nrow = length(times))
  if (ncol(Z) != dgm$q) {
    stop("`z()` returned ", ncol(Z), " columns but `G` is ",
         dgm$q, " x ", dgm$q, ".")
  }
  Z %*% dgm$G %*% t(Z)
}

#' Covariance of the linear predictor
#'
#' `Z G Z'`, the covariance induced on the linear predictor scale by
#' the random effects. Defined for every family, unlike [cov_at()],
#' which is the covariance of the response and therefore Gaussian-only.
#'
#' @param dgm A `dgm_conditional`.
#' @param times Numeric vector of observation times.
#' @return A covariance matrix on the linear-predictor scale.
#' @examples
#' # A Poisson GLMM with a random intercept and slope. `cov_at()` is
#' # undefined here, but the linear predictor still has a covariance.
#' g <- dgm_conditional(G = diag(c(4, 0.05)), sigma2 = NULL,
#'                      family = poisson())
#' linpred_cov(g, times = c(0, 3, 6))
#' @export
linpred_cov <- function(dgm, times) {
  stopifnot(inherits(dgm, "dgm_conditional"))
  .zgz(dgm, times)
}

#' @export
cov_at.dgm_marginal <- function(dgm, times) {
  if (!is.null(dgm$times)) {
    idx <- match(times, dgm$times)
    if (anyNA(idx)) {
      stop("A marginal DGM given as a matrix can only be evaluated ",
           "on its own grid. Times not on the grid: ",
           paste(utils::head(times[is.na(idx)], 5), collapse = ", "),
           ". Supply `V` as a function of time for irregular or ",
           "extended grids.")
    }
    return(dgm$V[idx, idx, drop = FALSE])
  }
  dgm$Vfun(times)
}

#' Convert a DGM to the marginal parameterization
#'
#' Always exact: `V = Z G Z' + sigma^2 I`.
#'
#' @param dgm A `dgm`.
#' @param times Grid on which to evaluate.
#' @return A `dgm_marginal`.
#' @export
as_marginal <- function(dgm, times) {
  dgm_marginal(cov_at(dgm, times), times = times)
}

#' Convert a DGM to the conditional parameterization
#'
#' Possible only when the target covariance is reachable from `q`
#' random effects. Uses [certify()] and refuses rather than silently
#' returning an approximation.
#'
#' @param dgm A `dgm`.
#' @param times Grid on which to evaluate.
#' @param q Integer. Random-effects budget. Defaults to the minimum
#'   required.
#' @return A `dgm_conditional` whose `z()` returns the certified
#'   loadings on `times`.
#' @export
as_conditional <- function(dgm, times, q = NULL) {
  V <- cov_at(dgm, times)
  cert <- certify(V)
  if (is.null(q)) q <- cert$q_min
  if (q < cert$q_min) {
    err <- reach_error(V, q = q)
    stop("Target covariance is not reachable from ", q,
         " random effect(s): it requires q_min = ", cert$q_min,
         ".\n  Best attainable relative Frobenius error at q = ", q,
         " is ", format(err$rel_frobenius, digits = 3),
         ".\n  Use `reach_error()` to inspect, or raise `q`.")
  }
  Zc <- cert$Z
  dgm_conditional(
    G = cert$G, sigma2 = cert$sigma2,
    z = function(t) {
      idx <- match(t, times)
      if (anyNA(idx)) {
        stop("Certified loadings are defined only on the grid used ",
             "by `as_conditional()`.")
      }
      Zc[idx, , drop = FALSE]
    }
  )
}
