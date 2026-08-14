library(tinytest)

# Time-to-event: single and recurrent. Generation is by inversion of
# the cumulative hazard, so the distributional claims below are exact
# statements about the mechanism, not approximations.

has_surv <- requireNamespace("survival", quietly = TRUE)

# ---- single event ---------------------------------------------------

set.seed(4001)
n <- 4000L
subj <- data.frame(id = seq_len(2 * n),
                   trt = rep(0:1, each = n))

# Exponential baseline, no covariate effect: median = scale * log(2).
g_exp <- dgm_tte(shape = 1, scale = 20)
e0 <- generate_tte(subj, g_exp, beta = c(trt = 0), followup = 1e6)
expect_true(all(c("time", "status") %in% names(e0)),
            info = 'single event returns time and status')
expect_true(all(e0$status == 1L),
            info = 'no censoring when follow-up is effectively infinite')
expect_equal(stats::median(e0$time), 20 * log(2), tolerance = 0.05,
             info = 'exponential median is scale * log(2)')
expect_equal(mean(e0$time), 20, tolerance = 0.05,
             info = 'exponential mean is the scale')

# Weibull shape is recovered: -log S(t) should be linear in log t with
# slope equal to the shape.
g_wb <- dgm_tte(shape = 1.6, scale = 15)
e1 <- generate_tte(subj, g_wb, beta = c(trt = 0), followup = 1e6)
qs <- stats::quantile(e1$time, c(0.2, 0.4, 0.6, 0.8))
slope <- stats::coef(stats::lm(log(-log(1 - c(.2, .4, .6, .8))) ~
                                 log(as.numeric(qs))))[2]
expect_equal(as.numeric(slope), 1.6, tolerance = 0.05,
             info = 'Weibull shape recovered from the survival curve')

# Administrative censoring.
e2 <- generate_tte(subj, g_exp, beta = c(trt = 0), followup = 10)
expect_true(all(e2$time <= 10 + 1e-9),
            info = 'no time exceeds the follow-up window')
expect_equal(mean(e2$status), 1 - exp(-10 / 20), tolerance = 0.02,
             info = 'event rate matches 1 - S(followup)')

# Subject-specific follow-up, as produced by a common close-out.
fu <- runif(nrow(subj), 5, 30)
e3 <- generate_tte(subj, g_exp, beta = c(trt = 0), followup = fu)
expect_true(all(e3$time <= fu + 1e-9),
            info = 'per-subject follow-up respected')

# The hazard ratio is recovered.
if (has_surv) {
  hr <- 0.65
  e4 <- generate_tte(subj, g_wb, beta = c(trt = log(hr)),
                     followup = 30)
  fit <- survival::coxph(survival::Surv(time, status) ~ trt, data = e4)
  est <- unname(stats::coef(fit))
  se <- sqrt(unname(stats::vcov(fit)[1, 1]))
  expect_true(abs(est - log(hr)) < 4 * se,
              info = 'Cox model recovers the true log hazard ratio')
}

# ---- recurrent events -----------------------------------------------

set.seed(4002)
subj2 <- data.frame(id = 1:1500, trt = rep(0:1, each = 750))
g_rec <- dgm_tte(shape = 1, scale = 10, recurrent = TRUE)
r0 <- generate_tte(subj2, g_rec, beta = c(trt = 0), followup = 20)

expect_true(all(c("tstart", "tstop", "status", "enum") %in% names(r0)),
            info = 'recurrent returns counting-process format')
expect_true(all(r0$tstop > r0$tstart),
            info = 'every interval has positive length')

# Intervals for a subject must tile [0, followup] without gaps.
chk <- tapply(seq_len(nrow(r0)), r0$id, function(k) {
  z <- r0[k, ][order(r0$tstart[k]), ]
  isTRUE(all.equal(z$tstart[-1], z$tstop[-nrow(z)])) &&
    abs(z$tstart[1]) < 1e-9 &&
    abs(z$tstop[nrow(z)] - 20) < 1e-6
})
expect_true(all(unlist(chk)),
            info = 'intervals tile [0, followup] contiguously')

# Exactly one censoring row per subject, and it is the last.
last_status <- tapply(seq_len(nrow(r0)), r0$id, function(k) {
  z <- r0[k, ][order(r0$tstart[k]), ]
  z$status[nrow(z)]
})
expect_true(all(unlist(last_status) == 0L),
            info = 'each subject ends with a censored interval')

# A homogeneous Poisson process of rate 1/scale over [0, fu] has
# expected count fu / scale.
n_ev <- tapply(r0$status, r0$id, sum)
expect_equal(mean(unlist(n_ev)), 20 / 10, tolerance = 0.06,
             info = 'mean event count matches the Poisson rate')
expect_equal(stats::var(unlist(n_ev)), 20 / 10, tolerance = 0.15,
             info = 'without frailty, counts are Poisson (var = mean)')

# Frailty induces overdispersion, which is its purpose.
set.seed(4003)
g_fr <- dgm_tte(shape = 1, scale = 10, recurrent = TRUE,
                frailty_sd = 0.8)
r1 <- generate_tte(subj2, g_fr, beta = c(trt = 0), followup = 20)
n_ev_fr <- unlist(tapply(r1$status, r1$id, sum))
expect_true(stats::var(n_ev_fr) > 1.5 * mean(n_ev_fr),
            info = 'frailty produces overdispersed event counts')

# Recurrent-event rate ratio is recovered under Andersen-Gill.
if (has_surv) {
  set.seed(4004)
  rr <- 0.6
  r2 <- generate_tte(subj2, g_rec, beta = c(trt = log(rr)),
                     followup = 20)
  fit_ag <- survival::coxph(
    survival::Surv(tstart, tstop, status) ~ trt + survival::cluster(id),
    data = r2)
  est <- unname(stats::coef(fit_ag))[1]
  se <- sqrt(unname(stats::vcov(fit_ag)[1, 1]))
  expect_true(abs(est - log(rr)) < 4 * se,
              info = 'Andersen-Gill recovers the true rate ratio')
}

expect_error(generate_tte(subj, g_exp, beta = c(nope = 1),
                          followup = 10),
             info = 'unknown beta names are rejected')

# ---- endpoint derived from a longitudinal trajectory ----------------

set.seed(4005)
s <- trial_schedule(treatment = 6, interval = 3)
arm <- factor(rep(c("placebo", "active"), each = 300),
              levels = c("placebo", "active"))
d <- runin_design(s, arm, reference = "placebo")
gd <- dgm_conditional(G = diag(c(4, 0.02)), sigma2 = 2)
lo <- generate_outcomes(d, gd,
                        beta = c(x_slope = 0.35, x_trt_active = -0.20),
                        intercept = 10)

ev <- tte_from_trajectory(lo, threshold = 16, direction = "above")
expect_equal(nrow(ev), length(unique(lo$id)),
             info = 'one row per subject')
expect_true(all(ev$status %in% c(0L, 1L)),
            info = 'status is binary')
expect_true(all(ev$time <= max(s$time) + 1e-9),
            info = 'event times lie within the observed schedule')

# Non-crossers are censored at their last observed visit.
never <- ev$id[ev$status == 0L]
if (length(never)) {
  last_t <- tapply(lo$time[lo$id %in% never & !is.na(lo$y)],
                   lo$id[lo$id %in% never & !is.na(lo$y)], max)
  expect_equal(ev$time[match(as.integer(names(last_t)), ev$id)],
               as.numeric(last_t),
               info = 'non-crossers censored at last observed visit')
}

# The treated arm worsens more slowly, so it should cross less often.
#
# Note what is NOT asserted: that the mean crossing time among those
# who crossed is later in the treated arm. It is not, and the reason
# is instructive. Conditioning on having crossed selects the fastest
# progressors within each arm, and it selects more severely in the arm
# where crossing is rarer. Measured on this data set: 64% of placebo
# and 25% of active cross, a Cox hazard ratio of 0.29, yet the mean
# time among crossers is 12.36 against 12.16 -- indistinguishable.
#
# This is why a survival endpoint is analyzed through the hazard with
# censored observations retained, not through the mean event time of
# the subset that had events. It is also precisely the comparison
# compendium 13 exists to study.
arm_of <- arm[ev$id]
p_pl <- mean(ev$status[arm_of == "placebo"])
p_ac <- mean(ev$status[arm_of == "active"])
expect_true(p_ac < p_pl,
            info = 'treatment lowers the proportion crossing')

if (has_surv) {
  fit_x <- survival::coxph(survival::Surv(time, status) ~ arm_of,
                           data = ev)
  lhr <- unname(stats::coef(fit_x))
  se_x <- sqrt(unname(stats::vcov(fit_x)[1, 1]))
  expect_true(lhr < 0 && abs(lhr) > 2 * se_x,
              info = paste('treatment reduces the hazard of crossing,',
                           'which is the correct comparison here'))
}

# Interpolation moves crossing times earlier than the visit grid.
ev_i <- tte_from_trajectory(lo, threshold = 16, direction = "above",
                            interpolate = TRUE)
both <- ev$status == 1L & ev_i$status == 1L
expect_true(all(ev_i$time[both] <= ev$time[both] + 1e-9),
            info = 'interpolated crossing is never later than the visit')

# Dropout propagates: a subject with all values missing is censored.
lo2 <- lo
lo2$y[lo2$id == 1] <- NA_real_
ev2 <- tte_from_trajectory(lo2, threshold = 16)
expect_equal(ev2$status[ev2$id == 1], 0L,
             info = 'a subject with no observations is censored')

# ---- piecewise-exponential baseline ---------------------------------

g_pw <- dgm_tte(breaks = c(6, 12), rates = c(0.02, 0.06, 0.15))
expect_true(g_pw$piecewise, info = 'piecewise baseline recorded')

# H0 is piecewise linear with the stated slopes, and H0inv inverts it.
expect_equal(zzrctsim:::.H0(0, g_pw), 0, info = 'H0(0) = 0')
expect_equal(zzrctsim:::.H0(6, g_pw), 6 * 0.02, tolerance = 1e-12,
             info = 'first interval accrues at its own rate')
expect_equal(zzrctsim:::.H0(12, g_pw), 6 * 0.02 + 6 * 0.06,
             tolerance = 1e-12,
             info = 'second interval accrues at its own rate')
expect_equal(zzrctsim:::.H0(20, g_pw), 6 * 0.02 + 6 * 0.06 + 8 * 0.15,
             tolerance = 1e-12, info = 'final open interval')
tt_grid <- c(0.5, 3, 6, 9, 12, 18, 30)
expect_equal(zzrctsim:::.H0inv(zzrctsim:::.H0(tt_grid, g_pw), g_pw),
             tt_grid,
             tolerance = 1e-10,
             info = 'H0inv is the exact inverse of H0')

# The realized hazard in each interval matches the specified rate.
set.seed(4101)
subj3 <- data.frame(id = 1:40000, trt = 0)
e_pw <- generate_tte(subj3, g_pw, beta = c(trt = 0), followup = 20)
# Occurrence/exposure rate on an interval: events divided by
# person-time actually spent in it. `tolerance` in tinytest is
# relative, so 0.05 means "within 5%", which is the right order for a
# Monte Carlo estimate at this sample size.
emp_rate <- function(lo, hi) {
  ev <- sum(e_pw$status == 1L & e_pw$time > lo & e_pw$time <= hi)
  pt <- sum(pmax(0, pmin(e_pw$time, hi) - lo))
  ev / pt
}
expect_equal(emp_rate(0, 6), 0.02, tolerance = 0.05,
             info = 'realized hazard matches rate on interval 1')
expect_equal(emp_rate(6, 12), 0.06, tolerance = 0.05,
             info = 'realized hazard matches rate on interval 2')
expect_equal(emp_rate(12, 20), 0.15, tolerance = 0.05,
             info = 'realized hazard matches rate on interval 3')

expect_error(dgm_tte(breaks = c(6, 12), rates = c(0.02, 0.06)),
             info = 'rates must be one longer than breaks')
expect_error(dgm_tte(breaks = c(12, 6), rates = c(1, 2, 3)),
             info = 'breaks must be increasing')
expect_error(dgm_tte(breaks = c(6)),
             info = 'breaks without rates is an error')

# ---- delayed treatment effect (non-proportional hazards) ------------

set.seed(4102)
subj4 <- data.frame(id = 1:12000, trt = rep(0:1, each = 6000))
lag <- 8
g_lag <- dgm_tte(shape = 1, scale = 20, effect_lag = lag)
e_lag <- generate_tte(subj4, g_lag, beta = c(trt = log(0.4)),
                      followup = 40)

# Before the lag the arms are identical by construction.
pre0 <- mean(e_lag$time[e_lag$trt == 0] <= lag &
               e_lag$status[e_lag$trt == 0] == 1L)
pre1 <- mean(e_lag$time[e_lag$trt == 1] <= lag &
               e_lag$status[e_lag$trt == 1] == 1L)
expect_true(abs(pre0 - pre1) < 0.02,
            info = 'no arm separation before the effect lag')

# After the lag the treated arm has the lower hazard.
post0 <- mean(e_lag$status[e_lag$trt == 0 & e_lag$time > lag] == 1L)
post1 <- mean(e_lag$status[e_lag$trt == 1 & e_lag$time > lag] == 1L)
expect_true(post1 < post0,
            info = 'treated hazard is lower after the effect lag')

if (has_surv) {
  # A Cox model fitted to lagged data violates proportional hazards,
  # and the test should detect it. This is the point of the option:
  # it generates data a constant hazard ratio cannot represent.
  fz <- survival::coxph(survival::Surv(time, status) ~ trt,
                        data = e_lag)
  ph <- survival::cox.zph(fz)
  expect_true(ph$table[1, "p"] < 0.01,
              info = 'delayed effect is detected as non-proportional')

  # And the fitted constant HR is attenuated toward 1 relative to the
  # true post-lag hazard ratio of 0.4.
  expect_true(exp(unname(stats::coef(fz))) > 0.4,
              info = 'a constant HR understates the post-lag effect')

  # With no lag, proportionality holds.
  set.seed(4103)
  g_nolag <- dgm_tte(shape = 1, scale = 20)
  e_nolag <- generate_tte(subj4, g_nolag, beta = c(trt = log(0.4)),
                          followup = 40)
  fz0 <- survival::coxph(survival::Surv(time, status) ~ trt,
                         data = e_nolag)
  expect_true(survival::cox.zph(fz0)$table[1, "p"] > 0.01,
              info = 'no lag: proportional hazards not rejected')
}

# Piecewise baseline and delayed effect compose.
set.seed(4104)
g_both <- dgm_tte(breaks = c(6, 12), rates = c(0.02, 0.06, 0.15),
                  effect_lag = 6, recurrent = TRUE)
r_both <- generate_tte(subj2, g_both, beta = c(trt = log(0.5)),
                       followup = 24)
expect_true(all(r_both$tstop > r_both$tstart),
            info = 'piecewise + lag + recurrent produces valid intervals')
expect_true(nrow(r_both) > nrow(subj2),
            info = 'recurrent events generated under a piecewise baseline')
