# Reachability of a marginal covariance from a random-effects model.
#
# A random-effects model with q random effects and i.i.d. residuals
# induces V = Z G Z' + sigma^2 I, where Z G Z' is positive
# semi-definite of rank at most q. A target V is therefore reachable
# only if some sigma^2 >= 0 makes V - sigma^2 I positive semi-definite
# with rank at most q.
#
# Subtracting sigma^2 I shifts every eigenvalue by the same amount, so
# the question is settled by the eigenvalue spectrum. Taking sigma^2 as
# the smallest eigenvalue zeroes as many eigenvalues as its
# multiplicity and leaves the rest non-negative; the resulting rank is
# the minimum q.
#
# The construction is classical: V = Lambda Lambda' + sigma^2 I solved
# by eigendecomposition is probabilistic PCA (Tipping and Bishop,
# 1999), equivalently factor analysis with isotropic uniquenesses.
# Nothing here is new mathematics. Its value is as a routine, automatic
# check inside a simulation workflow.

# Shared input guard for the certificate functions. Both `certify()` and
# `reach_error()` interpret the eigenvalue spectrum of `V` and are
# meaningless on an input that is not a symmetric PSD covariance: the
# spectrum of a non-symmetric matrix need not even be real, and a
# negative eigenvalue is reported as a negative residual variance, which
# is not a result but a silent contradiction.
.check_square_psd <- function(V, tol = 1e-9, arg = "V") {
  if (nrow(V) != ncol(V)) {
    stop("`", arg, "` must be a square covariance matrix; it is ",
         nrow(V), " x ", ncol(V),
         ". Supply a symmetric matrix with equal row and column counts.")
  }
  if (!isSymmetric(unname(V), tol = sqrt(.Machine$double.eps))) {
    stop("`", arg, "` must be symmetric, but it is not. Symmetrize it ",
         "with `(V + t(V)) / 2` if the asymmetry is numerical noise.")
  }
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < -tol) {
    stop("`", arg, "` must be positive semi-definite; its smallest ",
         "eigenvalue is ", format(min(ev), digits = 4),
         ". Supply a valid covariance matrix, or raise `tol` if the ",
         "negative eigenvalue is numerical noise.")
  }
  invisible(TRUE)
}

#' Certify whether a marginal covariance is reachable from random effects
#'
#' Computes the minimum number of random effects `q` required to
#' represent a target covariance as `Z G Z' + sigma^2 I`, together with
#' the construction that attains it.
#'
#' @param V A symmetric positive semi-definite covariance matrix. A
#'   singular `V` is accepted: it simply certifies a `sigma2` of zero,
#'   which is the noise-free case rather than an error.
#' @param tol Numeric. Eigenvalues below this, after shifting, are
#'   treated as zero. Also the tolerance by which a slightly negative
#'   eigenvalue is forgiven before `V` is rejected as not positive
#'   semi-definite.
#' @return An object of class `reach_certificate`, a list with
#'   \describe{
#'     \item{q_min}{integer. The minimum number of random effects that
#'       can represent `V` as `Z G Z' + sigma2 I`, namely the number of
#'       eigenvalues exceeding `sigma2` by more than `tol`.}
#'     \item{sigma2}{the residual variance attaining it: the smallest
#'       eigenvalue of `V`, floored at zero.}
#'     \item{Z}{a `p x q_min` matrix of loadings, the leading `q_min`
#'       eigenvectors of `V`. Orthonormal, and in general *not* design
#'       columns such as an intercept and a slope; see Details.}
#'     \item{G}{a `q_min x q_min` diagonal matrix whose entries are the
#'       leading eigenvalues shifted down by `sigma2`.}
#'     \item{eigenvalues}{the full spectrum of `V`, all `p` values in
#'       decreasing order.}
#'     \item{p}{integer, the dimension of `V`.}
#'   }
#'   [reconstruct()] returns `Z G Z' + sigma2 I` from these, equal to
#'   `V` to numerical precision whenever the certificate is exact.
#' @details
#' `q_min` is **necessary but not sufficient for a particular model**.
#' The certificate lets `Z` be arbitrary, namely the eigenvectors of
#' `V`. A real random-effects specification fixes `Z` to design columns
#' such as an intercept and a slope in time. Compound symmetry has a
#' flat leading eigenvector, which is why a random intercept reproduces
#' it exactly; AR(1) has a bowed one, so even the optimal
#' one-dimensional direction is unavailable to a random intercept.
#'
#' The certificate therefore proves impossibility cleanly -- no model
#' with `q < q_min` can represent the target -- but does not prove
#' possibility for a given `Z`. Use [reach_error()] to quantify how
#' close a given budget can get.
#' @examples
#' p <- 6
#' cs <- matrix(0.5, p, p); diag(cs) <- 1
#' certify(cs)$q_min          # 1: a random intercept suffices
#' ar1 <- 0.6 ^ abs(outer(1:p, 1:p, "-"))
#' certify(ar1)$q_min         # 5: not reachable parsimoniously
#' @export
certify <- function(V, tol = 1e-9) {
  V <- as.matrix(V)
  .check_square_psd(V, tol = tol)
  e <- eigen(V, symmetric = TRUE)
  s2 <- max(min(e$values), 0)
  lam <- e$values - s2
  q <- sum(lam > tol)

  out <- list(
    q_min = q,
    sigma2 = s2,
    Z = e$vectors[, seq_len(q), drop = FALSE],
    G = diag(lam[seq_len(q)], nrow = q),
    eigenvalues = e$values,
    p = nrow(V)
  )
  class(out) <- "reach_certificate"
  out
}

#' @export
print.reach_certificate <- function(x, ...) {
  cat("<reach_certificate>\n")
  cat("  p        =", x$p, "\n")
  cat("  q_min    =", x$q_min,
      if (x$q_min >= x$p - 1L && x$p > 2L) " (not parsimonious)" else "",
      "\n")
  cat("  sigma^2  =", format(x$sigma2, digits = 6), "\n")
  cat("  spectrum =", paste(format(x$eigenvalues, digits = 4),
                            collapse = " "), "\n")
  invisible(x)
}

#' Reconstruct the marginal covariance from a certificate
#'
#' @param cert A `reach_certificate`.
#' @return The covariance `Z G Z' + sigma^2 I`, equal to the original
#'   target to numerical precision.
#' @export
reconstruct <- function(cert) {
  if (!inherits(cert, "reach_certificate")) {
    stop("`cert` must be a `reach_certificate` object as returned by ",
         "`certify()`, but it has class ",
         paste(class(cert), collapse = "/"),
         ". Pass `certify(V)` rather than `V` itself.")
  }
  cert$Z %*% cert$G %*% t(cert$Z) +
    cert$sigma2 * diag(cert$p)
}

#' Error incurred by a restricted random-effects budget
#'
#' Best attainable rank-`q` plus `sigma^2 I` approximation to a target
#' covariance, and the resulting discrepancy. Converts "unreachable"
#' into a statement of how far wrong a given budget would be.
#'
#' @param V Target covariance.
#' @param q Integer vector of random-effect budgets to evaluate. Every
#'   element must lie in `1:(p-1)`: a budget of `p` or more reaches `V`
#'   exactly and has no approximation error to report. Defaults to
#'   `1:(p-1)`.
#' @param tol Numeric. Tolerance by which a slightly negative eigenvalue
#'   is forgiven before `V` is rejected as not positive semi-definite,
#'   as in [certify()].
#' @return A data frame with columns `q`, `sigma2`, `max_abs_error`,
#'   and `rel_frobenius`.
#' @details
#' The bound is optimistic for a *random intercept* specifically: it
#' uses the leading eigenvector of `V`, whereas a random intercept is
#' forced to the all-ones direction and can only do worse or equal.
#' @examples
#' ar1 <- 0.6 ^ abs(outer(1:6, 1:6, "-"))
#' reach_error(ar1)
#' @export
reach_error <- function(V, q = NULL, tol = 1e-9) {
  V <- as.matrix(V)
  .check_square_psd(V, tol = tol)
  p <- nrow(V)
  if (is.null(q)) q <- seq_len(p - 1L)
  # A budget of p or more is not an approximation problem: `V` is
  # reachable exactly and there is no residual to average over. Left
  # unchecked, `(qq + 1L):p` counts backwards and every reported
  # quantity is silently NA.
  if (!length(q) || !is.numeric(q) || anyNA(q)) {
    stop("`q` must be a non-empty numeric vector of random-effect ",
         "budgets, with no missing values.")
  }
  if (any(q < 1)) {
    stop("`q` must be at least 1; the smallest value supplied is ",
         min(q), ". A budget of zero random effects is not a rank-`q` ",
         "approximation.")
  }
  if (any(q >= p)) {
    stop("`q` must be less than `nrow(V)` = ", p,
         "; the largest value supplied is ", max(q),
         " and the maximum permitted is ", p - 1L,
         ". A budget of ", p, " or more reaches `V` exactly, so there ",
         "is no approximation error to report. Drop those values ",
         "from `q`, or use `certify()` for the exact construction.")
  }
  e <- eigen(V, symmetric = TRUE)
  lam <- e$values
  U <- e$vectors

  rows <- lapply(q, function(qq) {
    s2 <- mean(lam[(qq + 1L):p])
    lam_hat <- c(pmax(lam[seq_len(qq)], s2), rep(s2, p - qq))
    Vq <- U %*% diag(lam_hat) %*% t(U)
    data.frame(
      q = qq,
      sigma2 = s2,
      max_abs_error = max(abs(Vq - V)),
      rel_frobenius = norm(Vq - V, "F") / norm(V, "F")
    )
  })
  do.call(rbind, rows)
}
