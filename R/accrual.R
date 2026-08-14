# Staggered accrual and common close-out.
#
# `trial_schedule()` defines the schedule relative to a subject's own
# randomization. When subjects enrol at different calendar times and
# the study runs to a single close-out date, the realized schedule is
# subject specific: a subject enrolled early has more calendar time
# available before close-out and is therefore observed for longer than
# the nominal schedule, while the last enrolled subject is observed for
# exactly the nominal duration.
#
# Calendar time of subject i's visit j is  e_i + t_j,  where e_i is the
# enrolment time. A visit is observed when that falls at or before the
# close-out date C.

#' Generate enrolment times
#'
#' @param n Integer. Number of subjects.
#' @param period Numeric. Length of the accrual window, in the study
#'   time unit.
#' @param pattern Character. `"uniform"` spreads enrolment evenly at
#'   random over the window; `"linear"` places subjects at
#'   deterministic equally spaced times; `"poisson"` draws
#'   inter-arrival times from an exponential distribution with the rate
#'   implied by `n` and `period`.
#' @return Numeric vector of `n` enrolment times, sorted.
#' @details
#' `"uniform"` and `"linear"` return times within `[0, period]`.
#' `"poisson"` does not: it accumulates exponential inter-arrival
#' times at rate `n / period`, so `period` sets the *expected* span of
#' the last enrolment and any particular realization may overrun it.
#' That is the intended behavior of a Poisson process with a target
#' rate, but it means a downstream [close_out()] driven by the last
#' enrolment can fall after `period`.
#' @export
accrue <- function(n, period,
                   pattern = c("uniform", "linear", "poisson")) {
  pattern <- match.arg(pattern)
  stopifnot(n >= 1, period >= 0)
  e <- switch(
    pattern,
    uniform = stats::runif(n, 0, period),
    linear = seq(0, period, length.out = n),
    poisson = {
      gaps <- stats::rexp(n, rate = n / max(period, .Machine$double.eps))
      cumsum(gaps)
    }
  )
  sort(e)
}

#' Determine the study close-out date
#'
#' @param enroll Numeric vector of enrolment times.
#' @param schedule A `trial_schedule`.
#' @param rule Character. `"lslv"` (last subject, last visit) closes
#'   when the last enrolled subject completes the nominal schedule.
#'   `"fixed"` uses `at` directly.
#' @param at Numeric. Calendar close-out date, required when
#'   `rule = "fixed"`.
#' @return A single numeric calendar time.
#' @export
close_out <- function(enroll, schedule, rule = c("lslv", "fixed"),
                      at = NULL) {
  rule <- match.arg(rule)
  stopifnot(inherits(schedule, "trial_schedule"))
  if (rule == "fixed") {
    if (is.null(at)) stop("`at` is required when rule = 'fixed'.")
    return(as.numeric(at))
  }
  max(enroll) + max(schedule$time)
}

#' Realize a schedule under staggered accrual and a common close-out
#'
#' Expands a nominal `trial_schedule` across subjects with given
#' enrolment times, and marks which visits fall before the close-out
#' date. When `extend = TRUE` (the default) subjects enrolled before
#' the last one continue to be observed past the nominal final visit,
#' on the same spacing, up to close-out.
#'
#' @param schedule A `trial_schedule` with equal spacing. Extension
#'   requires a known `interval`, so schedules built from explicit
#'   `times` cannot be extended.
#' @param enroll Numeric vector of enrolment times, one per subject.
#' @param close Numeric. Close-out calendar time, from `close_out()`.
#' @param extend Logical. Continue observing early enrollees past the
#'   nominal final visit, up to close-out.
#' @return A data frame with one row per subject-visit actually
#'   observed, with columns `id`, `enroll`, `index`, `time` (study time
#'   since randomization), `calendar` (`enroll + time`), `phase`, `h`,
#'   and `on_treatment`. Visits after close-out are dropped.
#' @details
#' Under `rule = "lslv"` every subject completes at least the nominal
#' schedule, and follow-up duration is a decreasing function of
#' enrolment time. The resulting imbalance is informative: it is a
#' design feature of common close-out trials, and it is the mechanism
#' by which the information fraction at an interim analysis is
#' determined.
#'
#' Visits during the extension inherit the phase of the nominal final
#' visit. If the nominal schedule ends in a common close, extended
#' visits are also `common_close`; if it ends on treatment, they are
#' `treatment`.
#' @export
realize_schedule <- function(schedule, enroll, close,
                             extend = TRUE) {
  stopifnot(inherits(schedule, "trial_schedule"))
  interval <- attr(schedule, "interval")
  J1 <- attr(schedule, "J1")
  if (extend && (is.na(interval) || is.null(interval))) {
    stop("Extension requires an equally spaced schedule; ",
         "rebuild with `interval` rather than `times`.")
  }

  last_idx <- max(schedule$index)
  last_phase <- schedule$phase[schedule$index == last_idx]

  rows <- lapply(seq_along(enroll), function(i) {
    e <- enroll[i]
    idx <- schedule$index
    t_j <- schedule$time
    ph <- schedule$phase
    h <- schedule$h
    ot <- schedule$on_treatment

    if (extend) {
      avail <- close - e
      n_extra <- floor(avail / interval) - last_idx
      if (!is.na(n_extra) && n_extra >= 1) {
        k <- seq_len(n_extra)
        idx <- c(idx, last_idx + k)
        t_j <- c(t_j, (last_idx + k) * interval)
        ph <- c(as.character(ph), rep(as.character(last_phase),
                                      length(k)))
        h <- c(h, rep(1L, length(k)))
        ot <- c(ot, rep(FALSE, length(k)))
      } else {
        ph <- as.character(ph)
      }
    } else {
      ph <- as.character(ph)
    }

    cal <- e + t_j
    keep <- cal <= close + .Machine$double.eps^0.5
    data.frame(
      id = i, enroll = e,
      index = idx[keep], time = t_j[keep],
      calendar = cal[keep], phase = ph[keep],
      h = h[keep], on_treatment = ot[keep]
    )
  })

  out <- do.call(rbind, rows)
  out$phase <- factor(out$phase,
                      levels = c("run_in", "randomization",
                                 "treatment", "common_close"))
  rownames(out) <- NULL
  out
}
