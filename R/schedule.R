# Visit schedules spanning run-in, treatment, and common-close phases.
#
# The three-phase design follows the formulation in the companion
# compendium `res/04-runin-power-analysis`, whose closed-form variance
# results are derived under exactly this parameterization:
#
#   run-in         j = -J0, ..., -1   all participants untreated
#   randomization  j = 0              baseline
#   treatment      j = 1, ..., J1     randomized to treatment or placebo
#   common close   j = J1+1, ..., J1+J2   all participants off-treatment
#
# with t_j = j * interval, so run-in times are negative and the
# randomization visit sits at t = 0.

#' Construct a trial visit schedule
#'
#' Defines the observation schedule for a longitudinal trial spanning up
#' to three phases. Either give the phase counts with a common spacing,
#' or give explicit visit times for an unequally spaced design.
#'
#' @param run_in Integer. Number of pre-randomization visits (`J0`).
#' @param treatment Integer. Number of post-randomization on-treatment
#'   visits (`J1`).
#' @param common_close Integer. Number of post-treatment visits during
#'   which all arms are off treatment (`J2`).
#' @param interval Numeric. Spacing between visits, in the time unit of
#'   the study. Ignored when `times` is supplied.
#' @param times Numeric vector. Explicit visit times for unequally
#'   spaced designs. Must contain `0`, the randomization visit. When
#'   supplied, `treatment_end` determines where the common close begins.
#' @param treatment_end Numeric. The time of the last on-treatment
#'   visit. Required with `times` when a common close is present;
#'   defaults to the final time, i.e. no common close.
#' @return An object of class `trial_schedule`: a data frame with one
#'   row per visit and columns
#'   \describe{
#'     \item{index}{visit index `j`, negative during run-in, `0` at
#'       randomization}
#'     \item{time}{visit time `t_j`, negative during run-in}
#'     \item{phase}{factor: `run_in`, `randomization`, `treatment`,
#'       `common_close`}
#'     \item{h}{phase indicator, `0` for `j <= 0` and `1` for `j > 0`}
#'     \item{on_treatment}{logical, `TRUE` for `1 <= j <= J1`}
#'   }
#'   with attributes `J0`, `J1`, `J2`, and `interval`.
#' @export
trial_schedule <- function(run_in = 0L,
                           treatment = 4L,
                           common_close = 0L,
                           interval = 1,
                           times = NULL,
                           treatment_end = NULL) {
  if (is.null(times)) {
    # Split rather than checked jointly: the four have different
    # admissible ranges and different fixes, and a joint failure that
    # named none of them left the caller to guess which one was wrong.
    .bad_scalar <- function(x) {
      !is.numeric(x) || length(x) != 1L || is.na(x)
    }
    .saw <- function(x) {
      if (.bad_scalar(x)) {
        paste0("it is ", paste(class(x), collapse = "/"),
               " of length ", length(x))
      } else {
        paste0("it is ", format(x))
      }
    }
    if (.bad_scalar(run_in) || run_in < 0) {
      stop("`run_in` must be a single count of pre-randomization ",
           "visits of 0 or more, but ", .saw(run_in),
           ". Use 0 for a trial with no run-in phase.",
           call. = FALSE)
    }
    if (.bad_scalar(treatment) || treatment < 1) {
      stop("`treatment` must be a single count of post-randomization ",
           "on-treatment visits of 1 or more, but ", .saw(treatment),
           ". A trial needs at least one visit after randomization ",
           "to estimate a treatment effect.",
           call. = FALSE)
    }
    if (.bad_scalar(common_close) || common_close < 0) {
      stop("`common_close` must be a single count of post-treatment ",
           "visits of 0 or more, but ", .saw(common_close),
           ". Use 0 for a trial with no common close phase.",
           call. = FALSE)
    }
    if (.bad_scalar(interval) || interval <= 0) {
      stop("`interval` must be a single positive visit spacing, but ",
           .saw(interval),
           ". Give the spacing in the study time unit, e.g. ",
           "interval = 3 for visits every three months, or supply ",
           "`times` for an unequally spaced design.",
           call. = FALSE)
    }
    # Visit counts are counts. A fractional value would otherwise be
    # truncated silently by `as.integer()` below, so that
    # `treatment = 2.5` quietly became 2.
    frac <- c(run_in = run_in, treatment = treatment,
              common_close = common_close)
    bad <- frac[frac != trunc(frac)]
    if (length(bad)) {
      stop("Visit counts must be whole numbers; got ",
           paste(names(bad), "=", format(bad), collapse = ", "),
           ". Use `times` to specify a non-integer visit grid.")
    }
    j <- seq.int(-run_in, treatment + common_close)
    t_j <- j * interval
    J0 <- as.integer(run_in)
    J1 <- as.integer(treatment)
    J2 <- as.integer(common_close)
  } else {
    if (!is.numeric(times) || length(times) < 1L || anyNA(times)) {
      stop("`times` must be a non-empty numeric vector of visit ",
           "times with no missing values, but it is ",
           paste(class(times), collapse = "/"), " of length ",
           length(times),
           if (anyNA(times)) " containing NA" else "",
           ". Give the times on the study time scale, negative during ",
           "run-in and 0 at randomization.",
           call. = FALSE)
    }
    if (is.unsorted(times)) {
      out_of_order <- which(diff(times) < 0)[1L]
      stop("`times` must be in increasing order, but element ",
           out_of_order + 1L, " (", format(times[out_of_order + 1L]),
           ") is less than element ", out_of_order, " (",
           format(times[out_of_order]),
           "). Visit indices are taken from position, so sort the ",
           "times first: times = sort(c(",
           paste(format(times), collapse = ", "), ")).",
           call. = FALSE)
    }
    if (!any(abs(times) < .Machine$double.eps^0.5)) {
      stop("`times` must include 0, the randomization visit.")
    }
    t_j <- times
    zero <- which.min(abs(times))
    j <- seq_along(times) - zero
    if (is.null(treatment_end)) treatment_end <- max(times)
    J0 <- sum(j < 0L)
    J1 <- sum(j > 0L & t_j <= treatment_end)
    J2 <- sum(t_j > treatment_end)
    interval <- NA_real_
  }

  phase <- factor(
    ifelse(j < 0L, "run_in",
           ifelse(j == 0L, "randomization",
                  ifelse(j <= J1, "treatment", "common_close"))),
    levels = c("run_in", "randomization", "treatment",
               "common_close")
  )

  out <- data.frame(
    index = as.integer(j),
    time = t_j,
    phase = phase,
    h = as.integer(j > 0L),
    on_treatment = j >= 1L & j <= J1
  )
  attr(out, "J0") <- J0
  attr(out, "J1") <- J1
  attr(out, "J2") <- J2
  attr(out, "interval") <- interval
  class(out) <- c("trial_schedule", "data.frame")
  out
}

#' @export
print.trial_schedule <- function(x, ...) {
  cat("<trial_schedule>  J0 =", attr(x, "J0"),
      " J1 =", attr(x, "J1"),
      " J2 =", attr(x, "J2"), "\n")
  cat("  visits:", nrow(x),
      " time range: [", min(x$time), ",", max(x$time), "]\n")
  print(as.data.frame(x), row.names = FALSE)
  invisible(x)
}

#' Build the fixed-effect design columns for a run-in trial
#'
#' Constructs the three slope columns of the model
#' `Y_ij = alpha + a_i + (delta (1 - h_j) + beta h_j + gamma h_j g_i
#' + b_i) t_j + w_ij`, expanded over subjects.
#'
#' The `common_close` argument selects how the treatment column behaves
#' after the last on-treatment visit. The two conventions are not
#' equivalent and give different estimands; see Details.
#'
#' @param schedule Either a `trial_schedule` from [trial_schedule()],
#'   giving every subject the same visit grid, or a
#'   `realized_schedule` from [realize_schedule()], giving each
#'   subject the grid their enrollment date and the common close-out
#'   allow. The second is the staggered-accrual path, and produces a
#'   ragged design with differing numbers of rows per subject.
#' @param arm Arm allocation, one entry per subject, of length `n`.
#'   Any vector `factor()` accepts: a factor is used as given, and
#'   anything else (integer, character, logical) is coerced with
#'   `factor()`. Arbitrarily many levels are supported, not just two,
#'   so a three-arm dose-ranging design is written by supplying three
#'   levels. The level ordering matters only through `reference`, which
#'   defaults to the first level; every other level gets its own
#'   treatment column.
#' @param hinge Logical, length 1 or one entry per arm level. Whether
#'   each arm's mean profile changes slope at randomization. `FALSE`
#'   for the reference arm constrains its post-randomization slope to
#'   equal the run-in slope (`beta = delta`), which makes the run-in
#'   observations directly informative about the control slope. See
#'   Details.
#' @param reference The arm level treated as the reference (control).
#'   Defaults to the first sorted level.
#' @param common_close Character. `"revert"` (the default) sets the
#'   treatment column to zero during the common close, following
#'   `res/04-runin-power-analysis`. `"retain"` holds the accumulated
#'   treatment contribution constant after the last on-treatment visit.
#' @return A data frame with one row per subject-visit and columns as
#'   below. From a `trial_schedule` there are `n * nrow(schedule)`
#'   rows, subjects varying slowest and visits fastest. From a
#'   `realized_schedule` the frame is ragged -- each subject keeps the
#'   visits their enrollment date allowed -- and the row order is that
#'   of the realized schedule.
#'
#'   Note when analyzing a ragged design: `fit_ancova()` defaults to
#'   `max(dat$time)`, which under staggered accrual only the earliest
#'   enrollees reach. Pass `visit_time` explicitly, normally the last
#'   nominal visit of the original schedule, or the fit will be based
#'   on a handful of subjects and may return [null_fit()].
#'   \describe{
#'     \item{id}{subject index. From a `trial_schedule`, `1` to `n` in
#'       the order `arm` was given; from a `realized_schedule`, the
#'       `id` of that schedule, with `arm` matched to it by position}
#'     \item{arm}{the arm factor, recycled across that subject's
#'       visits}
#'     \item{index, time, phase}{copied from `schedule`}
#'     \item{x_slope}{the common slope column, equal to `time` on every
#'       visit}
#'     \item{x_hinge}{`h * time`, the reference arm's post-randomization
#'       slope increment. Present only when the reference arm is hinged
#'       **and** the schedule has a run-in phase (`J0 > 0`). With
#'       `J0 = 0` the column would be aliased with `x_slope`, since
#'       `h` differs from `1` only at the randomization visit where
#'       `time` is zero, so it is omitted rather than handed over as a
#'       rank-deficient design.}
#'     \item{x_trt_<level>}{one column per non-reference arm level that
#'       is hinged, holding the treatment-phase time for subjects in
#'       that arm and zero otherwise. Non-hinged arms are skipped, so
#'       there are as many such columns as there are hinged
#'       non-reference levels, not one per non-reference level.}
#'   }
#' @details
#' Under `"revert"` the treatment indicator is replaced by
#' `g_i * I(1 <= j <= J1)`, so the mean for a treated subject returns to
#' the control trajectory at the start of the common close. This is a
#' discontinuity in the mean: the accumulated benefit `gamma * t_{J1}`
#' is removed at once. It is the convention under which the closed-form
#' variance results of the companion compendium were derived, and it
#' makes the common-close visits informative about `beta` and the
#' variance components while contributing nothing to `gamma`.
#'
#' Under `"retain"` the treatment column is held at its value at the
#' last on-treatment visit, so both arms progress with slope `beta`
#' thereafter while the treated arm keeps its accumulated advantage.
#' This is the more common clinical reading of an off-treatment
#' follow-up. It is offered because the choice is a modeling decision
#' rather than a detail, and it should be made explicitly.
#'
#' The slope for a subject in arm `k` is
#' \itemize{
#'   \item run-in (`j <= 0`): `delta`, common to all arms;
#'   \item treatment (`1 <= j <= J1`): `delta + theta_ref + theta_k`;
#'   \item common close (`j > J1`): `delta + theta_ref`, all arms.
#' }
#' `theta_ref` is the reference arm's slope change at randomization and
#' is dropped when `hinge` is `FALSE` for the reference, constraining
#' that arm to a single straight line through baseline. `theta_k` is
#' arm `k`'s increment over the reference and is the treatment effect
#' on the rate of change. This is a reparameterization of the
#' `delta`/`beta`/`gamma` form: `theta_ref = beta - delta` and
#' `theta_k = gamma`, spanning the same column space when every arm is
#' hinged.
#' @export
runin_design <- function(schedule,
                         arm,
                         hinge = TRUE,
                         reference = NULL,
                         common_close = c("revert", "retain")) {
  ragged <- inherits(schedule, "realized_schedule")
  if (!ragged && !inherits(schedule, "trial_schedule")) {
    stop("`schedule` must be a `trial_schedule` from ",
         "`trial_schedule()`, or a `realized_schedule` from ",
         "`realize_schedule()` for a staggered-accrual design; ",
         "got an object of class ",
         paste(class(schedule), collapse = "/"), ".")
  }
  common_close <- match.arg(common_close)

  arm_f <- if (is.factor(arm)) arm else factor(arm)
  levs <- levels(arm_f)
  if (is.null(reference)) reference <- levs[1]
  reference <- as.character(reference)
  if (!reference %in% levs) {
    stop("`reference` must be one of: ", paste(levs, collapse = ", "))
  }

  if (length(hinge) == 1L) {
    hinge <- stats::setNames(rep(hinge, length(levs)), levs)
  } else {
    # Names are required rather than inferred. Assigning `levs` to an
    # unnamed vector matched it to *sorted level* order, not to the
    # order the arms were written, so `hinge = c(TRUE, FALSE, TRUE)`
    # for arms low/hi/pbo was silently applied as hi/low/pbo -- model
    # misspecification with no diagnostic. It also errored obscurely
    # from `names<-` whenever the lengths disagreed, because that
    # assignment ran before the length check below.
    if (is.null(names(hinge))) {
      stop("`hinge` must be length 1, or named for every arm level (",
           paste(levs, collapse = ", "), "). An unnamed vector of ",
           "length ", length(hinge), " would be matched to sorted ",
           "level order rather than to the order the arms were ",
           "given.")
    }
    if (!all(levs %in% names(hinge))) {
      stop("`hinge` is missing entries for arm level(s): ",
           paste(setdiff(levs, names(hinge)), collapse = ", "),
           ". Supply one entry per level, or a single value for all.")
    }
  }

  n <- length(arm_f)
  J1 <- attr(schedule, "J1")
  J0 <- attr(schedule, "J0")

  if (ragged) {
    # Subjects already have their own grids, of differing length, so
    # the columns are taken row-wise and `arm` is broadcast per
    # subject rather than per row.
    #
    # `arm` is matched to subjects by IDENTITY, never by rank. Ranking
    # is only correct for ids that are the integers 1..n: sorting a
    # factor orders by level rather than label, and sorting character
    # ids orders "S10" before "S2", each of which silently assigns
    # arms to the wrong subjects. Requiring 1..n here keeps the
    # contract identical to the balanced branch, where id is
    # `rep(seq_len(n), each = p)` and `arm[i]` belongs to `id == i`.
    if (is.null(attr(schedule, "J0")) ||
        is.null(attr(schedule, "J1")) ||
        is.null(attr(schedule, "t_last"))) {
      stop("`schedule` is classed as a `realized_schedule` but has ",
           "lost the phase-boundary attributes `J0`, `J1` and ",
           "`t_last`. Subsetting columns drops them. Rebuild it with ",
           "`realize_schedule()`, or subset rows only.")
    }
    raw_id <- schedule$id
    if (is.factor(raw_id) || is.character(raw_id)) {
      stop("`id` in a realized schedule must be integer, but it is ",
           class(raw_id)[1], ". `arm` is matched to subject `id` ",
           "1..n, and sorting a ", class(raw_id)[1], " id would key ",
           "the allocation to level or lexicographic order instead. ",
           "Use `as.integer()` on `id`, keeping `arm` in the same ",
           "subject order.")
    }
    ids <- sort(unique(raw_id))
    if (n != length(ids)) {
      stop("`arm` has length ", n, " but the realized schedule holds ",
           length(ids), " subject(s); supply one arm per subject.")
    }
    if (!isTRUE(all.equal(as.numeric(ids), as.numeric(seq_len(n))))) {
      stop("`id` in a realized schedule must be the integers 1 to ",
           n, ", so that `arm[i]` belongs to subject `i`; got ",
           paste(utils::head(ids, 5), collapse = ", "),
           if (length(ids) > 5) ", ..." else "",
           ". Renumber the subjects, keeping `arm` in the same order.")
    }
    t_last <- attr(schedule, "t_last")
    subj_id <- schedule$id
    idx <- schedule$index
    t_j <- schedule$time
    h <- schedule$h
    ontrt <- schedule$on_treatment
    phase_v <- schedule$phase
    a <- arm_f[match(subj_id, ids)]
  } else {
    p <- nrow(schedule)
    t_last <- schedule$time[schedule$index == J1]
    subj_id <- rep(seq_len(n), each = p)
    idx <- rep(schedule$index, times = n)
    t_j <- rep(schedule$time, times = n)
    h <- rep(schedule$h, times = n)
    ontrt <- rep(schedule$on_treatment, times = n)
    phase_v <- rep(schedule$phase, times = n)
    a <- rep(arm_f, each = p)
  }

  # Treatment-phase time, respecting the common-close convention.
  t_trt <- if (common_close == "revert") {
    as.numeric(ontrt) * t_j
  } else {
    as.numeric(h) * pmin(t_j, t_last)
  }

  out <- data.frame(
    id = subj_id,
    arm = a,
    index = idx,
    time = t_j,
    phase = phase_v,
    x_slope = t_j
  )

  # Reference arm's slope change at randomization, applying to every
  # arm from randomization onward (including the common close).
  #
  # Only meaningful when there are pre-randomization visits. With
  # J0 = 0 every observation has h = 1 except the randomization visit
  # itself, where t = 0, so h * t and t coincide exactly and the two
  # columns are aliased: there is no pre-randomization slope for the
  # hinge to depart from. Emitting the column anyway would hand the
  # user a rank-deficient design whose fit silently drops a term.
  if (isTRUE(hinge[[reference]]) && J0 > 0L) {
    out$x_hinge <- h * t_j
  }

  # One treatment-effect column per non-reference, hinged arm.
  for (k in setdiff(levs, reference)) {
    if (!isTRUE(hinge[[k]])) next
    out[[paste0("x_trt_", k)]] <- as.numeric(a == k) * t_trt
  }

  out
}
