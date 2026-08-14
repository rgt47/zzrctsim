# Missingness as two composable layers.
#
# A *mask* says which observations are missing. A *filter* applies a
# mask to data. Keeping them separate means a mask may be generated
# stochastically, imported from a completed trial, or constructed by
# hand, and any of those can be applied to any simulated data set with
# the same design.
#
# An earlier single-call wrapper fused generation and application. That
# made an externally supplied pattern impossible to use, and made it
# awkward to hold the missingness pattern fixed across replicates under
# common random numbers, so it was removed in favor of this pair.

# The response family a data set was generated from, or "gaussian" when
# unrecorded (data not produced by generate_outcomes()).
.response_family <- function(dat) {
  f <- attr(dat, "family")
  if (is.null(f)) "gaussian" else f
}

# Response-dependent missingness reads `y` on a continuous scale and
# centers it. That is well defined for a Gaussian response and needs a
# deliberate choice otherwise, so a non-Gaussian family is allowed only
# when `center` is supplied explicitly. MCAR is always allowed, because
# its hazard never touches `y`.
.check_family_for_hazard <- function(dat, psi1, psi2, center) {
  fam <- .response_family(dat)
  if (identical(fam, "gaussian")) return(invisible(TRUE))
  if (psi1 == 0 && psi2 == 0) return(invisible(TRUE))
  if (!is.null(center)) return(invisible(TRUE))
  stop("Response-dependent missingness (psi1 or psi2 non-zero) on a '",
       fam, "' response requires `center` to be supplied explicitly.\n",
       "  The default centering is the mean of `y`, and `psi1` is a ",
       "change in log-odds per unit of `y`; neither is interpretable ",
       "on a non-Gaussian scale without a deliberate choice.\n",
       "  MCAR (psi1 = psi2 = 0) is available for every family.",
       call. = FALSE)
}

#' Generate a missingness mask
#'
#' Draws a missingness pattern from a Diggle-Kenward selection model.
#' The mask is a first-class object: save it, inspect it, reuse it
#' across replicates, or apply it to a different data set with the same
#' design.
#'
#' @param dat Complete data, as returned by [generate_outcomes()].
#' @param target Numeric. Target proportion of subjects with any
#'   missing observation. Either a single value, or a vector named by
#'   the levels of `by` for group-specific targets.
#' @param psi0 Numeric. Dropout logit intercept, supplied directly
#'   instead of calibrating from `target`. May be named by the levels
#'   of `by`.
#' @param psi1,psi2 Numeric. Coefficients on the previously observed
#'   and the current, unobserved response, each centered at `center`.
#'   Non-zero `psi1` gives MAR; non-zero `psi2` gives MNAR. See
#'   Details.
#' @param psi_cov Named numeric vector of coefficients on additional
#'   columns of `dat`, giving covariate-dependent missingness. Columns
#'   are used as supplied, without centering.
#' @param by Character. Column of `dat` defining groups with separate
#'   `target` or `psi0`, typically `"arm"`.
#' @param monotone Logical. `TRUE` gives absorbing dropout: once
#'   missing, missing thereafter. `FALSE` gives intermittent
#'   missingness, each eligible visit independent.
#' @param center Numeric. Centering value for the response.
#' @param from Numeric. Missingness only permitted where `time > from`.
#' @return An object of class `missing_mask`: a data frame of `id`,
#'   `time`, and logical `missing`, with a `spec` attribute.
#' @details
#' The hazard follows Diggle and Kenward (1994) \doi{10.2307/2986113}:
#' for a subject still under observation at visit `j`,
#' \deqn{\mathrm{logit}\, P(\mathrm{drop\ at\ } j) =
#'   \psi_0 + \psi_1 (y_{j-1} - c) + \psi_2 (y_j - c).}
#' MCAR, MAR and MNAR are the nested restrictions
#' `psi1 = psi2 = 0`, `psi2 = 0`, and `psi2` free. MNAR therefore
#' *extends* MAR rather than replacing its history dependence.
#'
#' At a visit with no predecessor there is no previously observed
#' response to condition on, so the `psi1` term contributes nothing
#' there rather than falling back on the current response, which would
#' make a mechanism reported as MAR behave as MNAR at that visit.
#'
#' Centering at `c` matters: on an untransformed ADAS-Cog or CDR-SB
#' scale an uncentered response saturates the logit, so that any
#' non-trivial `psi1` drives the dropout proportion to one regardless
#' of the intended level.
#'
#' Calibration solves for `psi0` such that the expected proportion of
#' subjects ever missing equals `target`. Because the hazards depend
#' only on complete-data responses, that expectation is available in
#' closed form as `1 - prod(1 - p_j)` averaged over subjects, and is
#' strictly increasing in `psi0`, so no inner Monte Carlo loop is
#' needed. With `by` supplied, calibration is solved separately within
#' each group.
#' @seealso [apply_mask()] to apply the mask to data,
#'   [mask_from_data()] to import an empirical pattern, and
#'   [reference_based()] for post-discontinuation trajectories.
#' @export
dropout_mask <- function(dat,
                         target = 0.25,
                         psi0 = NULL,
                         psi1 = 0,
                         psi2 = 0,
                         psi_cov = NULL,
                         by = NULL,
                         monotone = TRUE,
                         center = NULL,
                         from = 0) {
  stopifnot(is.data.frame(dat), all(c("id", "time", "y") %in% names(dat)))
  .check_family_for_hazard(dat, psi1, psi2, center)
  if (is.null(center)) center <- mean(dat$y, na.rm = TRUE)

  if (!is.null(psi_cov)) {
    miss <- setdiff(names(psi_cov), names(dat))
    if (length(miss)) {
      stop("`psi_cov` names not found in `dat`: ",
           paste(miss, collapse = ", "))
    }
  }

  grp <- if (is.null(by)) {
    factor(rep("all", nrow(dat)))
  } else {
    if (!by %in% names(dat)) stop("`by` column not found: ", by)
    factor(dat[[by]])
  }
  levs <- levels(grp)

  .expand <- function(v, what) {
    if (length(v) == 1L && is.null(names(v))) {
      return(stats::setNames(rep(v, length(levs)), levs))
    }
    if (is.null(names(v))) names(v) <- levs
    if (!all(levs %in% names(v))) {
      stop("`", what, "` must be length 1 or named for every level of `by`.")
    }
    v[levs]
  }

  # Calibrate or accept psi0, group by group.
  psi0_g <- if (!is.null(psi0)) {
    .expand(psi0, "psi0")
  } else {
    tg <- .expand(target, "target")
    vapply(levs, function(g) {
      sub <- dat[grp == g, , drop = FALSE]
      f <- function(p) {
        .expected_missing(sub, p, psi1, psi2, psi_cov, center, from,
                          monotone) - tg[[g]]
      }
      if (f(-50) > 0 || f(50) < 0) {
        stop("Cannot calibrate psi0 for group '", g,
             "' to a target of ", tg[[g]], ".")
      }
      stats::uniroot(f, lower = -50, upper = 50, tol = 1e-8)$root
    }, numeric(1))
  }

  expected <- vapply(levs, function(g) {
    .expected_missing(dat[grp == g, , drop = FALSE], psi0_g[[g]],
                      psi1, psi2, psi_cov, center, from, monotone)
  }, numeric(1))

  # Draw the pattern.
  missing <- logical(nrow(dat))
  for (i in unique(dat$id)) {
    k <- which(dat$id == i)
    k <- k[order(dat$time[k])]
    g <- as.character(grp[k[1]])
    p <- .haz_rows(dat, k, psi0_g[[g]], psi1, psi2, psi_cov, center,
                   from)
    elig <- which(dat$time[k] > from)
    if (!length(elig)) next
    u <- stats::runif(length(elig))
    hit <- u < p
    if (!any(hit)) next
    if (monotone) {
      m <- elig[which(hit)[1]]
      missing[k[m:length(k)]] <- TRUE
    } else {
      missing[k[elig[hit]]] <- TRUE
    }
  }

  out <- data.frame(id = dat$id, time = dat$time, missing = missing)
  attr(out, "spec") <- list(
    psi0 = psi0_g, psi1 = psi1, psi2 = psi2, psi_cov = psi_cov,
    by = by, center = center, from = from, monotone = monotone,
    mechanism = dropout_mechanism(psi1, psi2),
    target = if (is.null(psi0)) target else NA_real_,
    expected = expected
  )
  class(out) <- c("missing_mask", "data.frame")
  out
}

# Hazards for one subject's eligible visits.
.haz_rows <- function(dat, k, psi0, psi1, psi2, psi_cov, center, from) {
  yv <- dat$y[k]
  tv <- dat$time[k]
  elig <- which(tv > from)
  if (!length(elig)) return(numeric(0))
  # At a visit with no predecessor there is no previously observed
  # response to condition on, so the history term contributes nothing:
  # `prev` is set to `center`, making `psi1 * (prev - center)` zero.
  # Substituting the *current* response here instead would make the
  # hazard depend on an unobserved value, which is MNAR, while the
  # mechanism would still be reported as MAR.
  prev <- ifelse(elig > 1L, yv[pmax(elig - 1L, 1L)], center)
  eta <- psi0 + psi1 * (prev - center) + psi2 * (yv[elig] - center)
  if (!is.null(psi_cov)) {
    for (nm in names(psi_cov)) {
      eta <- eta + psi_cov[[nm]] * as.numeric(dat[[nm]][k][elig])
    }
  }
  stats::plogis(eta)
}

# Expected proportion of subjects with any missing observation.
# Exact given complete data: monotone and intermittent coincide here,
# since both are "at least one hazard fires".
.expected_missing <- function(dat, psi0, psi1, psi2, psi_cov, center,
                              from, monotone) {
  ids <- unique(dat$id)
  mean(vapply(ids, function(i) {
    k <- which(dat$id == i)
    k <- k[order(dat$time[k])]
    p <- .haz_rows(dat, k, psi0, psi1, psi2, psi_cov, center, from)
    if (!length(p)) 0 else 1 - prod(1 - p)
  }, numeric(1)))
}

#' @export
print.missing_mask <- function(x, ...) {
  sp <- attr(x, "spec")
  cat("<missing_mask>", sp$mechanism,
      if (sp$monotone) "monotone" else "intermittent", "\n")
  cat("  observations missing:", sum(x$missing), "/", nrow(x), "\n")
  cat("  subjects affected:",
      length(unique(x$id[x$missing])), "/",
      length(unique(x$id)), "\n")
  cat("  psi0:", paste(names(sp$psi0), format(sp$psi0, digits = 4),
                       sep = "=", collapse = "  "), "\n")
  invisible(x)
}

#' Build a missingness mask from an observed data set
#'
#' Imports an empirical missingness pattern so that simulated data can
#' be masked exactly as a completed trial was. The subject identifiers
#' and times must line up with the data the mask will be applied to.
#'
#' @param observed A data frame with `id`, `time`, and a response
#'   column.
#' @param response Character. Name of the response column whose `NA`
#'   pattern defines the mask.
#' @return A `missing_mask`.
#' @export
mask_from_data <- function(observed, response = "y") {
  stopifnot(all(c("id", "time", response) %in% names(observed)))
  out <- data.frame(id = observed$id, time = observed$time,
                    missing = is.na(observed[[response]]))
  mono <- all(vapply(split(out, out$id), function(z) {
    m <- z$missing[order(z$time)]
    !any(diff(m) < 0)
  }, logical(1)))
  attr(out, "spec") <- list(
    psi0 = NA_real_, psi1 = NA_real_, psi2 = NA_real_,
    psi_cov = NULL, by = NULL, center = NA_real_, from = NA_real_,
    monotone = mono, mechanism = "empirical",
    target = NA_real_,
    expected = mean(tapply(out$missing, out$id, any))
  )
  class(out) <- c("missing_mask", "data.frame")
  out
}

#' Apply a missingness mask to data
#'
#' @param dat Complete data.
#' @param mask A `missing_mask`, or a logical vector of `nrow(dat)`.
#' @param response Character. Column to set to `NA`.
#' @return `dat` with masked entries set to `NA` and a `missing`
#'   logical column added.
#' @export
apply_mask <- function(dat, mask, response = "y") {
  if (inherits(mask, "missing_mask")) {
    key_d <- paste(dat$id, dat$time)
    key_m <- paste(mask$id, mask$time)
    idx <- match(key_d, key_m)
    if (anyNA(idx)) {
      stop("Mask does not cover every row of `dat`: ",
           sum(is.na(idx)), " unmatched id/time combinations.")
    }
    m <- mask$missing[idx]
  } else {
    stopifnot(is.logical(mask), length(mask) == nrow(dat))
    m <- mask
  }
  dat[[response]][m] <- NA_real_
  dat$missing <- m
  dat
}

#' Reference-based post-discontinuation trajectories
#'
#' Modifies the mean trajectory after treatment discontinuation
#' according to a pattern-mixture assumption. This is a *generator*:
#' it produces data whose true post-withdrawal behavior follows the
#' named assumption, so that analyses which assume something else can
#' be evaluated for bias.
#'
#' @param dat Complete data from [generate_outcomes()], carrying `beta`
#'   and `intercept` attributes.
#' @param mask A `missing_mask` identifying discontinuation. The first
#'   masked visit for each subject is taken as the discontinuation
#'   point.
#' @param method Character. `"j2r"` jump to reference, `"cir"` copy
#'   increments in reference, `"cr"` copy reference, `"lmcf"` last mean
#'   carried forward.
#' @param trt_cols Character vector of treatment-effect columns.
#'   Defaults to columns beginning `x_trt_`.
#' @return `dat` with `y` and `mu` modified after discontinuation, and
#'   a `discontinued` logical column added.
#' @details
#' The reference mean is the fitted mean with every treatment column
#' set to zero, so it is the trajectory the subject would have followed
#' in the control arm. Modifications are applied to the mean and the
#' subject's realized residual is preserved, which keeps the
#' within-subject correlation intact.
#'
#' Note that discontinuation and missingness are distinct: under a
#' treatment-policy estimand post-discontinuation values are observed,
#' under a hypothetical estimand they are not. This function alters the
#' trajectory only. Use [apply_mask()] separately if the values are
#' also to be hidden.
#' @export
reference_based <- function(dat, mask,
                            method = c("j2r", "cir", "cr", "lmcf"),
                            trt_cols = NULL) {
  method <- match.arg(method)
  fam <- .response_family(dat)
  if (!identical(fam, "gaussian")) {
    stop("`reference_based()` is defined only for a Gaussian response; ",
         "this data set was generated from '", fam, "'.\n",
         "  The construction preserves each subject's additive residual ",
         "`y - mu`, which is not meaningful for a count or binary ",
         "outcome and would produce values outside the support.",
         call. = FALSE)
  }
  beta <- attr(dat, "beta")
  intercept <- attr(dat, "intercept")
  if (is.null(beta)) {
    stop("`dat` must carry a `beta` attribute; use generate_outcomes().")
  }
  if (is.null(trt_cols)) {
    trt_cols <- grep("^x_trt_", names(beta), value = TRUE)
  }
  if (!length(trt_cols)) {
    stop("No treatment columns found; supply `trt_cols`.")
  }

  beta_ref <- beta
  beta_ref[trt_cols] <- 0
  X <- as.matrix(dat[, names(beta), drop = FALSE])
  mu_ref <- intercept + as.vector(X %*% beta_ref)

  resid <- dat$y - dat$mu
  mu_new <- dat$mu
  disc <- logical(nrow(dat))

  key_d <- paste(dat$id, dat$time)
  key_m <- paste(mask$id, mask$time)
  idx <- match(key_d, key_m)
  # Same check as `apply_mask()`: unmatched rows would otherwise give
  # NA in `is_missing` and fail later at `if (!any(mi))` with
  # "missing value where TRUE/FALSE needed".
  if (anyNA(idx)) {
    stop("Mask does not cover every row of `dat`: ",
         sum(is.na(idx)), " unmatched id/time combinations.")
  }
  is_missing <- mask$missing[idx]

  for (i in unique(dat$id)) {
    k <- which(dat$id == i)
    k <- k[order(dat$time[k])]
    mi <- is_missing[k]
    if (!any(mi)) next
    d <- which(mi)[1]
    at <- k[d:length(k)]
    disc[at] <- TRUE
    if (method == "cr") at <- k
    prev <- if (d > 1L) k[d - 1L] else k[1]

    mu_new[at] <- switch(
      method,
      j2r  = mu_ref[at],
      cr   = mu_ref[at],
      lmcf = dat$mu[prev],
      cir  = dat$mu[prev] + (mu_ref[at] - mu_ref[prev])
    )
  }

  dat$mu <- mu_new
  dat$y <- mu_new + resid
  dat$discontinued <- disc
  dat
}
