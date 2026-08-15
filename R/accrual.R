# Staggered accrual and common close-out.
#
# `trial_schedule()` defines the schedule relative to a subject's own
# randomization. When subjects enroll at different calendar times and
# the study runs to a single close-out date, the realized schedule is
# subject specific: a subject enrolled early has more calendar time
# available before close-out and is therefore observed for longer than
# the nominal schedule, while the last enrolled subject is observed for
# exactly the nominal duration.
#
# Calendar time of subject i's visit j is  e_i + t_j,  where e_i is the
# enrollment time. A visit is observed when that falls at or before the
# close-out date C.

#' Generate enrollment times
#'
#' @param n Integer. Number of subjects.
#' @param period Numeric. Length of the accrual window, in the study
#'   time unit.
#' @param pattern Character. `"uniform"` spreads enrollment evenly at
#'   random over the window; `"linear"` places subjects at
#'   deterministic equally spaced times; `"poisson"` draws
#'   inter-arrival times from an exponential distribution with the rate
#'   implied by `n` and `period`.
#' @return Numeric vector of `n` enrollment times, sorted.
#' @details
#' `"uniform"` and `"linear"` return times within `[0, period]`.
#' `"poisson"` does not: it accumulates exponential inter-arrival
#' times at rate `n / period`, so `period` sets the *expected* span of
#' the last enrollment and any particular realization may overrun it.
#' That is the intended behavior of a Poisson process with a target
#' rate, but it means a downstream [close_out()] driven by the last
#' enrollment can fall after `period`.
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
#' @param enroll Numeric vector of enrollment times.
#' @param schedule A `trial_schedule`.
#' @param rule Character. `"lslv"` (last subject, last visit) closes
#'   when the last enrolled subject completes the nominal schedule.
#'   `"fixed"` uses `at` directly.
#' @param at Numeric. Calendar close-out date, required when
#'   `rule = "fixed"`.
#' @return A single numeric calendar time: `max(enroll) +
#'   max(schedule$time)` under `"lslv"`, or `as.numeric(at)` under
#'   `"fixed"`. Pass it to [realize_schedule()] as `close`.
#' @examples
#' set.seed(1)
#' s <- trial_schedule(treatment = 4, interval = 3)
#' e <- accrue(n = 20, period = 6, pattern = "uniform")
#'
#' # Last subject, last visit: close 12 months after the final
#' # enrollment, so every subject completes the nominal schedule.
#' close_out(e, s, rule = "lslv")
#'
#' # A fixed calendar date instead: subjects enrolled late are
#' # censored before completing.
#' close_out(e, s, rule = "fixed", at = 15)
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
#' enrollment times, and marks which visits fall before the close-out
#' date. When `extend = TRUE` (the default) subjects enrolled before
#' the last one continue to be observed past the nominal final visit,
#' on the same spacing, up to close-out.
#'
#' @param schedule A `trial_schedule` with equal spacing. Extension
#'   requires a known `interval`, so schedules built from explicit
#'   `times` cannot be extended.
#' @param enroll Numeric vector of enrollment times, one per subject.
#' @param close Numeric. Close-out calendar time, from `close_out()`.
#' @param extend Logical. Continue observing early enrollees past the
#'   nominal final visit, up to close-out.
#' @return A data frame with one row per subject-visit actually
#'   observed, and columns
#'   \describe{
#'     \item{id}{subject index, `1` to `length(enroll)`, in the order
#'       `enroll` was given}
#'     \item{enroll}{that subject's enrollment time, repeated across
#'       their visits}
#'     \item{index}{visit index `j`, continuing past the nominal final
#'       index when the visit is an extension}
#'     \item{time}{study time since randomization, `index * interval`}
#'     \item{calendar}{`enroll + time`, the calendar time of the visit}
#'     \item{phase}{factor with levels `run_in`, `randomization`,
#'       `treatment`, `common_close`; extended visits inherit the phase
#'       of the nominal final visit}
#'     \item{h}{phase indicator, `1` after randomization; `1` for all
#'       extended visits}
#'     \item{on_treatment}{logical; `FALSE` for all extended visits,
#'       since the nominal treatment period has ended}
#'   }
#'   Visits whose `calendar` exceeds `close` are dropped, so the number
#'   of rows per subject varies and the data frame is unbalanced by
#'   construction.
#' @details
#' Under `rule = "lslv"` every subject completes at least the nominal
#' schedule, and follow-up duration is a decreasing function of
#' enrollment time. The resulting imbalance is informative: it is a
#' design feature of common close-out trials, and it is the mechanism
#' by which the information fraction at an interim analysis is
#' determined.
#'
#' Visits during the extension inherit the phase of the nominal final
#' visit. If the nominal schedule ends in a common close, extended
#' visits are also `common_close`; if it ends on treatment, they are
#' `treatment`.
#' @examples
#' set.seed(1)
#' s <- trial_schedule(treatment = 4, interval = 3)
#'
#' # Accrue 20 subjects over 6 months, then close out when the last
#' # of them completes the nominal 12-month schedule.
#' e <- accrue(n = 20, period = 6, pattern = "uniform")
#' cl <- close_out(e, s, rule = "lslv")
#' rs <- realize_schedule(s, e, cl)
#' head(rs, 3)
#'
#' # Follow-up is longest for the first enrolled and exactly nominal
#' # for the last: the design feature that makes the data unbalanced.
#' range(tapply(rs$time, rs$id, max))
#'
#' # Without extension every subject gets the nominal schedule only.
#' fixed <- realize_schedule(s, e, cl, extend = FALSE)
#' range(table(fixed$id))
#' @export
realize_schedule <- function(schedule, enroll, close,
                             extend = TRUE) {
  if (!inherits(schedule, "trial_schedule")) {
    stop("`schedule` must be a `trial_schedule` from ",
         "`trial_schedule()`; got an object of class ",
         paste(class(schedule), collapse = "/"), ".")
  }
  if (!is.numeric(enroll) || !length(enroll) || anyNA(enroll)) {
    stop("`enroll` must be a non-empty numeric vector of enrollment ",
         "times with no missing values; see `accrue()`.")
  }
  if (!is.numeric(close) || length(close) != 1L || !is.finite(close)) {
    stop("`close` must be a single finite close-out time; see ",
         "`close_out()`.")
  }
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
    # A subject enrolled so late that even their baseline falls after
    # close-out keeps no visits at all. `data.frame()` would then try
    # to recycle the scalar `id` against zero-length columns and fail
    # with "arguments imply differing number of rows". Return NULL and
    # report the affected subjects once, below.
    if (!any(keep)) return(NULL)
    data.frame(
      id = i, enroll = e,
      index = idx[keep], time = t_j[keep],
      calendar = cal[keep], phase = ph[keep],
      h = h[keep], on_treatment = ot[keep]
    )
  })

  dropped <- which(vapply(rows, is.null, logical(1)))
  if (length(dropped)) {
    warning(length(dropped), " subject(s) enrolled after the ",
            "close-out date and contribute no visits (position(s) ",
            paste(utils::head(dropped, 5), collapse = ", "),
            if (length(dropped) > 5) ", ..." else "",
            "). They are dropped, so the realized schedule holds ",
            length(enroll) - length(dropped), " of ", length(enroll),
            " subjects and its `id` values are renumbered. Supply an ",
            "`arm` of the reduced length to `runin_design()`, or ",
            "extend `close`.")
    rows <- rows[-dropped]
    if (!length(rows)) {
      stop("No subject has any visit on or before the close-out ",
           "date; every enrollment falls after `close`.")
    }
    # Renumber so that ids remain 1..n, which is the contract
    # `runin_design()` matches `arm` against.
    for (k in seq_along(rows)) rows[[k]]$id <- k
  }

  out <- do.call(rbind, rows)
  out$phase <- factor(out$phase,
                      levels = c("run_in", "randomization",
                                 "treatment", "common_close"))
  rownames(out) <- NULL

  # The design-matrix builder needs the phase boundaries, which are a
  # property of the schedule rather than of any subject's realized
  # grid. Carrying them here is what lets `runin_design()` accept this
  # ragged frame directly, so that staggered accrual reaches the
  # generation stage instead of dead-ending at a data frame.
  class(out) <- c("realized_schedule", "data.frame")
  attr(out, "J0") <- attr(schedule, "J0")
  attr(out, "J1") <- J1
  attr(out, "interval") <- interval
  attr(out, "t_last") <- schedule$time[schedule$index == J1]
  out
}
