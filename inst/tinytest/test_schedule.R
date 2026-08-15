library(tinytest)

# Visit schedules: run-in, treatment, common close, staggered accrual.
#
# Semantics follow res/04-runin-power-analysis:
#   run-in        j = -J0..-1   all untreated, t_j < 0
#   randomization j = 0         t_j = 0
#   treatment     j = 1..J1
#   common close  j > J1        all arms off treatment

s <- trial_schedule(run_in = 2, treatment = 4, common_close = 2,
                    interval = 3)

# ---- schedule structure -------------------------------------------

expect_equal(nrow(s), 9L,
             info = 'J0 + 1 + J1 + J2 visits')
expect_equal(attr(s, "J0"), 2L, info = 'J0 recorded')
expect_equal(attr(s, "J1"), 4L, info = 'J1 recorded')
expect_equal(attr(s, "J2"), 2L, info = 'J2 recorded')
expect_equal(s$time[s$index == 0], 0,
             info = 'randomization visit sits at t = 0')
expect_true(all(s$time[s$index < 0] < 0),
            info = 'run-in times are negative')
expect_equal(sum(s$on_treatment), 4L,
             info = 'exactly J1 on-treatment visits')
expect_true(all(!s$on_treatment[s$phase == "common_close"]),
            info = 'common close is off treatment')
expect_true(all(s$h[s$index > 0] == 1L),
            info = 'phase indicator h = 1 after randomization')

# Explicit unequally spaced times.
su <- trial_schedule(times = c(-6, -3, 0, 2, 5, 9, 14),
                     treatment_end = 9)
expect_equal(attr(su, "J0"), 2L, info = 'J0 from explicit times')
expect_equal(attr(su, "J1"), 3L, info = 'J1 from treatment_end')
expect_equal(attr(su, "J2"), 1L, info = 'J2 from treatment_end')
expect_error(trial_schedule(times = c(1, 2, 3)),
             info = 'times must contain the randomization visit')

# ---- common-close convention --------------------------------------

arm <- factor(c("placebo", "active"), levels = c("placebo", "active"))

d_rev <- runin_design(s, arm, common_close = "revert")
g_rev <- d_rev$x_trt_active[d_rev$arm == "active"]
expect_equal(g_rev, c(0, 0, 0, 3, 6, 9, 12, 0, 0),
             info = 'revert: treatment column zeroes in common close')

d_ret <- runin_design(s, arm, common_close = "retain")
g_ret <- d_ret$x_trt_active[d_ret$arm == "active"]
expect_equal(g_ret, c(0, 0, 0, 3, 6, 9, 12, 12, 12),
             info = 'retain: treatment column holds at last on-trt value')

# The reference arm's post-randomization slope applies during the
# common close under both conventions: all arms revert to slope beta.
expect_equal(d_rev$x_hinge[d_rev$arm == "active"],
             c(0, 0, 0, 3, 6, 9, 12, 15, 18),
             info = 'reference hinge column continues through common close')

# ---- per-arm hinge -------------------------------------------------

d_h <- runin_design(s, arm, hinge = TRUE)
expect_true("x_hinge" %in% names(d_h),
            info = 'hinged reference contributes a slope-change column')

d_nh <- runin_design(s, arm,
                     hinge = c(placebo = FALSE, active = TRUE))
expect_false("x_hinge" %in% names(d_nh),
             info = 'unhinged reference drops the slope-change column')
expect_true("x_trt_active" %in% names(d_nh),
            info = 'treatment column survives an unhinged reference')

# Unhinged reference costs exactly one parameter.
np <- function(d) sum(grepl("^x_", names(d)))
expect_equal(np(d_h) - np(d_nh), 1L,
             info = 'unhinged reference is a one-parameter restriction')

# The hinged parameterization spans the same space as delta/beta/gamma.
h <- s$h; tt <- s$time; ot <- as.numeric(s$on_treatment)
orig <- cbind((1 - h) * tt, h * tt, ot * tt)
new <- as.matrix(d_h[d_h$arm == "active",
                     c("x_slope", "x_trt_active")])
expect_equal(qr(cbind(orig, new))$rank, qr(orig)$rank,
             info = paste('hinge parameterization is a reparameterization',
                          'of delta/beta/gamma, not a different model'))

# ---- staggered accrual and common close-out ------------------------

s2 <- trial_schedule(run_in = 1, treatment = 4, interval = 3)
e <- accrue(5, period = 12, pattern = "linear")
expect_equal(length(e), 5L, info = 'accrue returns one time per subject')
expect_false(is.unsorted(e), info = 'enrolment times are sorted')

cl <- close_out(e, s2, rule = "lslv")
expect_equal(cl, max(e) + max(s2$time),
             info = 'LSLV closes when the last enrollee completes')

r <- realize_schedule(s2, e, cl, extend = TRUE)
expect_true(all(r$calendar <= cl + 1e-8),
            info = 'no visit occurs after close-out')
expect_equal(sum(r$id == 5L), nrow(s2),
             info = 'last enrollee is observed for exactly the nominal schedule')
expect_true(sum(r$id == 1L) > nrow(s2),
            info = 'first enrollee is observed beyond the nominal schedule')

fu <- tapply(r$time, r$id, max)
expect_true(all(diff(as.numeric(fu)) <= 0),
            info = 'follow-up duration decreases with enrolment time')
expect_equal(as.numeric(fu[5]), max(s2$time),
             info = 'last enrollee gets exactly the nominal duration')

# Without extension every subject gets at most the nominal schedule.
r0 <- realize_schedule(s2, e, cl, extend = FALSE)
expect_true(all(table(r0$id) <= nrow(s2)),
            info = 'extend = FALSE caps at the nominal schedule')

# Explicit-times schedules cannot be extended: no spacing is defined.
expect_error(realize_schedule(su, e, cl, extend = TRUE),
             info = 'extension requires an equally spaced schedule')

# ---- staggered accrual reaches the generation stage -----------------
#
# `realize_schedule()` used to dead-end: it returned a plain ragged
# data frame that `runin_design()` refused, so the staggered-accrual
# branch advertised in DESCRIPTION had no route into a DGM. It now
# carries the phase boundaries forward and `runin_design()` accepts
# it, so this whole path must keep working end to end.

set.seed(7)
s_acc <- trial_schedule(run_in = 1, treatment = 4, interval = 3)
e_acc <- accrue(40, period = 12, pattern = "uniform")
cl_acc <- close_out(e_acc, s_acc, rule = "lslv")
rs_acc <- realize_schedule(s_acc, e_acc, cl_acc)

expect_inherits(rs_acc, "realized_schedule",
                info = 'realize_schedule is classed for runin_design')
expect_equal(attr(rs_acc, "J1"), attr(s_acc, "J1"),
             info = 'phase boundaries carried forward')

arm_acc <- factor(rep(c("placebo", "active"), each = 20),
                  levels = c("placebo", "active"))
d_acc <- runin_design(rs_acc, arm_acc, reference = "placebo")

expect_true(all(c("x_slope", "x_hinge", "x_trt_active") %in%
                  names(d_acc)),
            info = 'ragged design gets the same model columns')
expect_equal(length(unique(d_acc$id)), 40L,
             info = 'one subject per arm entry')
expect_true(length(unique(table(d_acc$id))) > 1L,
            info = 'design is genuinely ragged under staggered entry')
# `arm` is per subject, not per row: every row of a subject agrees.
expect_true(all(vapply(split(as.character(d_acc$arm), d_acc$id),
                       function(a) length(unique(a)) == 1L,
                       logical(1))),
            info = 'arm is constant within subject')

g_acc <- dgm_conditional(G = diag(c(9, 0.04)), sigma2 = 4)
dat_acc <- generate_outcomes(
  d_acc, g_acc,
  beta = c(x_slope = 0.5, x_hinge = -0.1, x_trt_active = -0.25),
  intercept = 20)
expect_equal(nrow(dat_acc), nrow(d_acc),
             info = 'generation accepts the ragged design')
expect_true(all(!is.na(dat_acc$y)),
            info = 'every ragged row gets an outcome')

# A non-schedule first argument is refused by name, not by a bare
# stopifnot echo.
expect_error(runin_design(data.frame(a = 1), arm_acc),
             pattern = "schedule",
             info = 'runin_design names the offending argument')

# `arm` maps to subjects by SORTED id, not by row order.
#
# Keying the mapping to first appearance would make a schedule sorted
# by calendar date -- the natural ordering for staggered accrual --
# silently assign arms to the wrong subjects.
arm_map <- function(x) {
  u <- unique(runin_design(x, arm_acc, reference = "placebo")[,
                                                    c("id", "arm")])
  as.character(u$arm[order(u$id)])
}
ref_map <- arm_map(rs_acc)

reord <- function(rows) {
  z <- rs_acc[rows, ]
  attributes(z)[c("J0", "J1", "interval", "t_last")] <-
    attributes(rs_acc)[c("J0", "J1", "interval", "t_last")]
  class(z) <- class(rs_acc)
  z
}
expect_equal(arm_map(reord(order(-rs_acc$id, rs_acc$time))), ref_map,
             info = 'reversed row order maps arms identically')
expect_equal(arm_map(reord(order(rs_acc$calendar))), ref_map,
             info = 'calendar-sorted rows map arms identically')
expect_equal(ref_map, as.character(arm_acc),
             info = 'arm i belongs to subject i, as in the balanced path')

# ---- ragged-path input validation -----------------------------------
#
# Each of these previously either crashed obscurely or, worse, ran
# silently with the treatment allocation wrong.

# A subject enrolled after close-out keeps no visits. Previously this
# crashed with "arguments imply differing number of rows: 1, 0".
e_late <- accrue(10, 6, "linear")
expect_warning(
  rs_late <- realize_schedule(s_acc, e_late,
                              close = close_out(e_late, s_acc,
                                                "fixed", at = 2)),
  info = 'subjects past close-out are dropped with a warning')
expect_true(nrow(rs_late) > 0L,
            info = 'the surviving subjects still yield a schedule')
expect_equal(sort(unique(rs_late$id)), seq_len(length(unique(rs_late$id))),
             info = 'ids are renumbered 1..n after dropping')

# `arm` is matched by identity, so non-integer ids are refused rather
# than ranked: sorting a factor uses level order and sorting character
# ids puts "S10" before "S2", each silently inverting the allocation.
rs_f <- rs_acc; rs_f$id <- factor(rs_f$id)
expect_error(runin_design(rs_f, arm_acc), pattern = "must be integer",
             info = 'factor id is refused, not ranked')
rs_c <- rs_acc; rs_c$id <- paste0("S", rs_c$id)
expect_error(runin_design(rs_c, arm_acc), pattern = "must be integer",
             info = 'character id is refused, not ranked')

# Column-subsetting keeps the class but drops the phase boundaries.
stripped <- rs_acc[, c("id", "enroll", "index", "time", "calendar",
                       "phase", "h", "on_treatment")]
expect_error(runin_design(stripped, arm_acc),
             pattern = "phase-boundary",
             info = 'a degraded realized_schedule is rejected by name')

# An unnamed multi-element `hinge` was matched to sorted level order.
expect_error(
  runin_design(trial_schedule(run_in = 2, treatment = 4, interval = 3),
               factor(rep(c("a", "b", "c"), each = 2)),
               hinge = c(TRUE, FALSE)),
  pattern = "named for every arm level",
  info = 'unnamed hinge is refused rather than silently reordered')

expect_error(realize_schedule(s_acc, numeric(0), close = 20),
             pattern = "non-empty numeric",
             info = 'empty enroll is refused')
expect_error(realize_schedule(s_acc, c(1, NA, 3), close = 20),
             pattern = "non-empty numeric",
             info = 'NA enroll is refused')
expect_error(realize_schedule(s_acc, c(1, 2), close = Inf),
             pattern = "single finite",
             info = 'non-finite close is refused')
