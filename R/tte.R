# Time-to-event endpoints: single and recurrent.
#
# Time-to-event does not fit the covariance-based parameterization of
# `dgm_conditional()` / `dgm_marginal()`. Censoring, hazards and risk
# sets are a different structure, not a different marginal law, so this
# is a sibling generator behind the same contract rather than another
# `family` argument.
#
# Generation is by inversion of the cumulative hazard, which is exact
# rather than a discretization. Under proportional hazards the
# individual cumulative hazard is
#
#     H_i(t) = H_0(t) exp(lp_i),        lp_i = x_i'beta + b_i
#
# so if E ~ Exp(1) then T = H_0^{-1}(E / exp(lp_i)) has exactly the
# intended distribution. For recurrent events the same inversion is
# applied to the successive partial sums of i.i.d. Exp(1) variates,
# which generates the event times of a non-homogeneous Poisson process
# with intensity h_i(t). Both are exact.

#' Specify a time-to-event data-generating mechanism
#'
#' @param shape Numeric. Weibull shape. `1` gives an exponential
#'   baseline hazard, values above `1` an increasing hazard, below `1`
#'   a decreasing one. Ignored when `breaks` and `rates` are supplied.
#' @param scale Numeric. Weibull scale, in the time unit of the study.
#'   With `shape = 1` this is the mean time to event at `lp = 0`.
#' @param frailty_sd Numeric. Standard deviation of a Gaussian frailty
#'   on the log-hazard scale, `b_i ~ N(0, frailty_sd^2)`. Zero gives no
#'   frailty. For recurrent events the frailty is what induces
#'   within-subject dependence between gap times.
#' @param recurrent Logical. `FALSE` generates time to first event;
#'   `TRUE` generates all events in the follow-up window.
#' @param breaks Numeric vector of interior cut points for a
#'   piecewise-exponential baseline hazard, e.g. `c(6, 12)` for the
#'   intervals `[0, 6)`, `[6, 12)`, `[12, Inf)`. Supply with `rates`.
#' @param rates Numeric vector of baseline hazard rates, one per
#'   interval, so `length(breaks) + 1` of them.
#' @param effect_lag Numeric. Delay before the covariate effect
#'   begins. Zero gives proportional hazards throughout; a positive
#'   value gives a **non-proportional** hazard in which the arms are
#'   indistinguishable until `effect_lag` and diverge thereafter.
#' @param max_events Integer. Safety cap on events per subject.
#' @return An object of class `dgm_tte`.
#' @details
#' Two baseline hazards are available. The Weibull, with
#' `H_0(t) = (t / scale)^shape` and `H_0^{-1}(u) = scale * u^(1/shape)`;
#' and the piecewise exponential, whose cumulative hazard is piecewise
#' linear and therefore invertible in closed form on each interval. The
#' piecewise form is the usual choice in trial simulation because it
#' expresses an arbitrary hazard shape without committing to a
#' parametric family.
#'
#' `effect_lag` breaks proportional hazards deliberately. The intensity
#' becomes `h_i(t) = h_0(t) exp(lp_i * I(t > lag))`, so
#'
#'     H_i(t) = H_0(min(t, lag)) +
#'              [H_0(t) - H_0(min(t, lag))] * exp(lp_i),
#'
#' which is still invertible in closed form and so is generated
#' exactly. A delayed separation of this kind is standard in
#' immuno-oncology and is the shape anti-amyloid trials in Alzheimer's
#' disease are usually argued to show; a constant hazard ratio cannot
#' represent it, and fitting one to such data estimates a
#' time-averaged quantity that corresponds to no instant of the trial.
#'
#' A frailty makes the marginal distribution non-proportional even
#' though the conditional model is proportional hazards, which is the
#' usual reason a fitted hazard ratio attenuates toward one. That is a
#' property of the mechanism, not a defect, and it is why the estimand
#' must state whether it names the conditional or the marginal hazard
#' ratio.
#' @export
dgm_tte <- function(shape = 1, scale = 12, breaks = NULL,
                    rates = NULL, frailty_sd = 0, effect_lag = 0,
                    recurrent = FALSE, max_events = 100L) {
  stopifnot(frailty_sd >= 0, effect_lag >= 0,
            is.logical(recurrent), max_events >= 1)

  piecewise <- !is.null(breaks) || !is.null(rates)
  if (piecewise) {
    if (is.null(breaks) || is.null(rates)) {
      stop("`breaks` and `rates` must be supplied together.")
    }
    if (length(rates) != length(breaks) + 1L) {
      stop("`rates` must have one more entry than `breaks`: ",
           length(breaks), " breaks require ", length(breaks) + 1L,
           " rates, got ", length(rates), ".")
    }
    if (is.unsorted(breaks, strictly = TRUE) || any(breaks <= 0)) {
      stop("`breaks` must be strictly increasing and positive.")
    }
    if (any(rates <= 0)) stop("`rates` must be positive.")
    cuts <- c(0, breaks)
    # Cumulative hazard at each cut point, for inversion.
    widths <- diff(c(cuts, Inf))
    cumH <- c(0, cumsum(rates[-length(rates)] *
                          diff(cuts)))
  } else {
    stopifnot(shape > 0, scale > 0)
    cuts <- cumH <- NULL
  }

  structure(list(shape = shape, scale = scale,
                 breaks = breaks, rates = rates,
                 cuts = cuts, cumH = cumH, piecewise = piecewise,
                 frailty_sd = frailty_sd, effect_lag = effect_lag,
                 recurrent = recurrent,
                 max_events = as.integer(max_events)),
            class = c("dgm_tte", "dgm"))
}

#' @export
print.dgm_tte <- function(x, ...) {
  cat("<dgm_tte>",
      if (x$recurrent) "recurrent events" else "time to first event",
      "\n")
  if (isTRUE(x$piecewise)) {
    cat("  baseline: piecewise exponential, breaks [",
        paste(x$breaks, collapse = ", "), "], rates [",
        paste(signif(x$rates, 4), collapse = ", "), "]\n", sep = "")
  } else {
    cat("  baseline: Weibull(shape =", x$shape,
        ", scale =", x$scale, ")",
        if (x$shape == 1) " [exponential]" else "", "\n", sep = "")
  }
  cat("  frailty SD:", x$frailty_sd,
      if (x$effect_lag > 0)
        paste0("   effect lag: ", x$effect_lag,
               " (non-proportional hazards)") else "", "\n")
  invisible(x)
}

# Baseline cumulative hazard and its inverse, for either baseline.
# Both are exact: the Weibull inverts analytically and the piecewise
# form is piecewise linear, so inversion is a lookup plus one division.
.H0 <- function(t, dgm) {
  if (!isTRUE(dgm$piecewise)) return((t / dgm$scale)^dgm$shape)
  j <- findInterval(t, dgm$cuts)
  dgm$cumH[j] + dgm$rates[j] * (t - dgm$cuts[j])
}

.H0inv <- function(u, dgm) {
  if (!isTRUE(dgm$piecewise)) return(dgm$scale * u^(1 / dgm$shape))
  j <- findInterval(u, dgm$cumH)
  dgm$cuts[j] + (u - dgm$cumH[j]) / dgm$rates[j]
}

# Inverse of the individual cumulative hazard, allowing a delayed
# covariate effect. With lag L the effect multiplies the hazard only
# beyond L, so hazard accrued before L is unmodified:
#   H_i(t) = H0(min(t,L)) + [H0(t) - H0(min(t,L))] * r
# Inverting for t given a target e on the Exp(1) scale.
.Hinv_ind <- function(e, rate, dgm) {
  lag <- dgm$effect_lag
  if (lag <= 0) return(.H0inv(e / rate, dgm))
  Hlag <- .H0(lag, dgm)
  ifelse(e <= Hlag,
         .H0inv(e, dgm),
         .H0inv(Hlag + (e - Hlag) / rate, dgm))
}

#' Generate time-to-event outcomes
#'
#' @param subjects A data frame with one row per subject, containing
#'   `id` and any columns named in `beta`.
#' @param dgm A `dgm_tte`.
#' @param beta Named numeric vector of log-hazard ratios. Names must
#'   match columns of `subjects`.
#' @param followup Numeric. Administrative censoring time. Either a
#'   single value applied to every subject, or a vector of length
#'   `nrow(subjects)` giving each subject's own follow-up, as produced
#'   by staggered accrual with a common close-out.
#' @return For a single event, `subjects` with `time` and `status`
#'   added (`status = 1` event, `0` censored). For recurrent events, a
#'   counting-process data frame with one row per interval and columns
#'   `id`, `tstart`, `tstop`, `status`, `enum`, plus the subject
#'   columns; this is the format `survival::coxph()` expects for
#'   Andersen-Gill models.
#' @details
#' Administrative censoring from a common close-out is the natural
#' companion to [accrue()] and [close_out()]: pass
#' `followup = close - enroll` and each subject is censored at the
#' study end, with follow-up decreasing in enrolment time.
#' @examples
#' set.seed(1)
#' subj <- data.frame(id = 1:200,
#'                    trt = rep(0:1, each = 100))
#' g <- dgm_tte(shape = 1.2, scale = 18)
#' generate_tte(subj, g, beta = c(trt = log(0.7)), followup = 24)[1:3, ]
#' @export
generate_tte <- function(subjects, dgm, beta, followup) {
  stopifnot(inherits(dgm, "dgm_tte"), is.data.frame(subjects),
            "id" %in% names(subjects))
  if (is.null(names(beta)) || any(names(beta) == "")) {
    stop("`beta` must be fully named.")
  }
  miss <- setdiff(names(beta), names(subjects))
  if (length(miss)) {
    stop("`beta` names not found in `subjects`: ",
         paste(miss, collapse = ", "))
  }
  n <- nrow(subjects)
  if (length(followup) == 1L) followup <- rep(followup, n)
  stopifnot(length(followup) == n, all(followup > 0))

  X <- as.matrix(subjects[, names(beta), drop = FALSE])
  lp <- as.vector(X %*% beta)
  if (dgm$frailty_sd > 0) {
    lp <- lp + stats::rnorm(n, 0, dgm$frailty_sd)
  }
  rate <- exp(lp)

  if (!dgm$recurrent) {
    e <- stats::rexp(n, rate = 1)
    t_event <- .Hinv_ind(e, rate, dgm)
    out <- subjects
    out$time <- pmin(t_event, followup)
    out$status <- as.integer(t_event <= followup)
    attr(out, "lp") <- lp
    return(out)
  }

  rows <- vector("list", n)
  for (i in seq_len(n)) {
    fu <- followup[i]
    cum <- 0
    times <- numeric(0)
    for (k in seq_len(dgm$max_events)) {
      cum <- cum + stats::rexp(1, rate = 1)
      tk <- .Hinv_ind(cum, rate[i], dgm)
      if (tk > fu) break
      times <- c(times, tk)
    }
    starts <- c(0, times)
    stops <- c(times, fu)
    status <- c(rep(1L, length(times)), 0L)
    # Drop a final zero-length interval when an event lands exactly at
    # the censoring time.
    keep <- stops > starts + .Machine$double.eps^0.5
    if (!any(keep)) keep[1] <- TRUE
    rows[[i]] <- data.frame(
      id = subjects$id[i],
      tstart = starts[keep], tstop = stops[keep],
      status = status[keep],
      enum = seq_len(sum(keep)),
      subjects[i, setdiff(names(subjects), "id"), drop = FALSE],
      row.names = NULL
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "lp") <- lp
  out
}

#' Derive a time-to-event endpoint from a longitudinal trajectory
#'
#' Defines the event as the first visit at which a generated
#' longitudinal response crosses a threshold, so that a survival
#' endpoint and a continuous endpoint describe the *same* underlying
#' process rather than being drawn independently.
#'
#' @param dat Long-format data from [generate_outcomes()].
#' @param threshold Numeric. The crossing value.
#' @param direction Character. `"above"` treats crossing upward as the
#'   event, `"below"` downward.
#' @param interpolate Logical. Linearly interpolate the crossing time
#'   between the bracketing visits, rather than using the visit time
#'   itself. Interval censoring at visit times is the honest default;
#'   interpolation is offered for comparison.
#' @return One row per subject with `id`, `time`, `status`.
#' @details
#' This is what compendium 13 requires: a latent continuous severity
#' process from which time to threshold crossing is derived. Because
#' both endpoints come from one trajectory, a comparison of MMRM
#' against survival analysis is a comparison of *analyses*, not of two
#' unrelated data-generating mechanisms.
#'
#' Subjects never crossing are censored at their last observed visit,
#' so dropout and administrative censoring both propagate correctly.
#' @export
tte_from_trajectory <- function(dat, threshold,
                                direction = c("above", "below"),
                                interpolate = FALSE) {
  direction <- match.arg(direction)
  stopifnot(all(c("id", "time", "y") %in% names(dat)))
  ids <- unique(dat$id)
  rows <- lapply(ids, function(i) {
    k <- which(dat$id == i)
    k <- k[order(dat$time[k])]
    y <- dat$y[k]
    tt <- dat$time[k]
    ok <- !is.na(y)
    if (!any(ok)) {
      return(data.frame(id = i, time = max(tt), status = 0L))
    }
    y <- y[ok]; tt <- tt[ok]
    hit <- if (direction == "above") y >= threshold else y <= threshold
    if (!any(hit)) {
      return(data.frame(id = i, time = max(tt), status = 0L))
    }
    m <- which(hit)[1]
    t_evt <- tt[m]
    if (interpolate && m > 1L) {
      y0 <- y[m - 1L]; y1 <- y[m]
      if (!isTRUE(all.equal(y1, y0))) {
        frac <- (threshold - y0) / (y1 - y0)
        frac <- min(max(frac, 0), 1)
        t_evt <- tt[m - 1L] + frac * (tt[m] - tt[m - 1L])
      }
    }
    data.frame(id = i, time = t_evt, status = 1L)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
