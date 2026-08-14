library(tinytest)

# Power by simulation: the single-size estimate, the curve, and the
# inversion of the curve for a target power. Every run below uses a
# fixed master seed, so each assertion is deterministic rather than
# probabilistic: `sim_power()` draws its replicate streams from
# `sim_streams(B, seed)`, and the same seed reproduces the same
# replicates exactly. Sizes, replicate counts and visit grids are kept
# deliberately small to keep the file fast.

s <- trial_schedule(treatment = 2, interval = 6)
gd <- dgm_conditional(G = diag(c(9, 0.04)), sigma2 = 4)

mk_design <- function(n) {
  arm <- factor(rep(c("placebo", "active"), each = n),
                levels = c("placebo", "active"))
  runin_design(s, arm, reference = "placebo")
}

b <- c(x_slope = 0.5, x_trt_active = -0.25)
b_null <- c(x_slope = 0.5, x_trt_active = 0)
tmax <- max(s$time)
es <- estimand("final-visit contrast", b[["x_trt_active"]] * tmax)
es0 <- estimand("null contrast", 0)
anc <- list(ancova = function(z) fit_ancova(z, reference = "placebo"))

# ---- sim_power -------------------------------------------------------

p20 <- sim_power(20L, mk_design, gd, b, es, anc, B = 100L,
                 intercept = 20, seed = 1L)

expect_true(is.data.frame(p20), info = 'sim_power returns a data frame')
expect_equal(nrow(p20), 1L,
             info = 'one row per method for a single method')
expect_equal(names(p20), c("n_per_arm", "method", "power", "mcse"),
             info = 'documented column contract')
expect_equal(p20$n_per_arm, 20L,
             info = 'the requested size is carried through')
expect_equal(p20$method, "ancova",
             info = 'the fitter name labels the row')
expect_true(p20$power >= 0 && p20$power <= 1,
            info = 'power is a proportion')
expect_true(p20$mcse > 0,
            info = 'a non-degenerate rejection rate has a positive MCSE')

# Two fitters give two rows, one per method.
p_two <- sim_power(20L, mk_design, gd, b, es,
                   analyse = list(
                     final = function(z) fit_ancova(z, reference = "placebo"),
                     mid = function(z) fit_ancova(z, visit_time = 6,
                                                  reference = "placebo")),
                   B = 50L, intercept = 20, seed = 1L)
expect_equal(nrow(p_two), 2L, info = 'one row per method')
expect_equal(sort(p_two$method), c("final", "mid"),
             info = 'both methods are reported')

# Reproducibility: the same seed gives the same answer.
expect_equal(sim_power(20L, mk_design, gd, b, es, anc, B = 100L,
                       intercept = 20, seed = 1L)$power,
             p20$power,
             info = 'sim_power is reproducible from its seed')

# `dropout` is plumbed through sim_power -> dropout_mask -> apply_mask
# -> run_simulation. Heavy dropout must both change the result and
# still yield a valid row.
p40 <- sim_power(40L, mk_design, gd, b, es, anc, B = 100L,
                 intercept = 20, seed = 1L)
p40_drop <- sim_power(40L, mk_design, gd, b, es, anc, B = 100L,
                      intercept = 20, seed = 1L,
                      dropout = list(target = 0.6, by = "arm"))

expect_equal(nrow(p40_drop), 1L,
             info = 'dropout run still returns a single valid row')
expect_true(p40_drop$power >= 0 && p40_drop$power <= 1,
            info = 'power under dropout is still a proportion')
expect_false(isTRUE(all.equal(p40_drop$power, p40$power)),
             info = '`dropout` changes the result: it is plumbed through')
expect_true(p40_drop$power < p40$power,
            info = 'heavy dropout reduces power relative to none')

# ---- power_curve -----------------------------------------------------

n_grid <- c(10L, 20L, 40L, 80L)
pc <- power_curve(n_grid, design_fn = mk_design, dgm = gd, beta = b,
                  estimand = es, analyse = anc, B = 100L,
                  intercept = 20, seed = 1L)

expect_inherits(pc, "power_curve", info = 'power_curve S3 class')
expect_inherits(pc, "data.frame", info = 'power_curve is a data frame')
expect_equal(nrow(pc), length(n_grid),
             info = 'one row per grid point for a single method')
expect_equal(pc$n_per_arm, n_grid,
             info = 'grid values appear in the supplied order')
expect_equal(rownames(pc), as.character(seq_len(nrow(pc))),
             info = 'row names are reset after rbind')
expect_true(all(pc$power >= 0 & pc$power <= 1),
            info = 'every curve point is a proportion')

# Monotonicity. This is not a flaky assertion: the seed is fixed, so
# the curve is a deterministic function of the inputs, and the effect
# and grid are wide enough (power rises from well under a half to one)
# that the increase is unambiguous. The tolerance still allows a
# single-point dip of one Monte Carlo standard error, so a change of
# RNG stream would not break the test spuriously.
expect_true(all(diff(pc$power) > -max(pc$mcse)),
            info = 'power is non-decreasing in n up to Monte Carlo error')
expect_true(pc$power[nrow(pc)] > pc$power[1],
            info = 'the largest size has more power than the smallest')

expect_true(is.character(capture.output(print(pc))[1]),
            info = 'print.power_curve runs')

# ---- sample_size: interpolation branch -------------------------------

ss <- sample_size(target = 0.80, n_grid = n_grid, confirm_B = 100L,
                  design_fn = mk_design, dgm = gd, beta = b,
                  estimand = es, analyse = anc, B = 100L,
                  intercept = 20, seed = 1L)

expect_inherits(ss, "sample_size", info = 'sample_size S3 class')
expect_equal(sort(names(ss)),
             sort(c("n_per_arm", "target", "curve", "confirmation")),
             info = 'documented element contract')
expect_equal(ss$target, 0.80, info = 'target is retained')
expect_true(is.integer(ss$n_per_arm) && ss$n_per_arm >= 1L,
            info = 'selected size is a positive integer')
expect_inherits(ss$curve, "power_curve",
                info = 'the curve is returned with the answer')
expect_true(ss$n_per_arm > min(n_grid) && ss$n_per_arm <= max(n_grid),
            info = 'the interpolated size lies inside the grid')

# The interpolated size brackets the target power on the curve.
below <- ss$curve$power < 0.80
expect_true(ss$n_per_arm > max(ss$curve$n_per_arm[below]),
            info = 'selection exceeds every grid size short of target')

expect_true(is.data.frame(ss$confirmation) &&
              nrow(ss$confirmation) == 1L,
            info = 'confirmation is a single sim_power row')
expect_equal(ss$confirmation$n_per_arm, ss$n_per_arm,
             info = 'confirmation is run at the selected size')

expect_true(is.character(capture.output(print(ss))[1]),
            info = 'print.sample_size runs')

# ---- sample_size: the confirmation uses a separate RNG stream --------

# The curve seed is offset by one for the confirmation run, so the
# confirmed power is an independent check of the selected size rather
# than a re-reading of the replicates that selected it. This is
# verified exactly: the confirmation must equal a `sim_power()` call at
# seed + 1, and must not equal one at the curve's own seed. Both are
# deterministic given the seeds, so neither assertion is flaky.
conf_offset <- sim_power(ss$n_per_arm, mk_design, gd, b, es, anc,
                         B = 100L, intercept = 20, seed = 2L)
conf_same <- sim_power(ss$n_per_arm, mk_design, gd, b, es, anc,
                       B = 100L, intercept = 20, seed = 1L)

expect_equal(ss$confirmation$power, conf_offset$power,
             info = 'confirmation runs at the curve seed offset by one')
expect_false(isTRUE(all.equal(conf_offset$power, conf_same$power)),
             info = 'the offset stream is genuinely a different sample')

# Reproducibility of the whole search.
ss_again <- sample_size(target = 0.80, n_grid = n_grid,
                        confirm_B = 100L, design_fn = mk_design,
                        dgm = gd, beta = b, estimand = es,
                        analyse = anc, B = 100L, intercept = 20,
                        seed = 1L)
expect_equal(ss_again$n_per_arm, ss$n_per_arm,
             info = 'the search is reproducible from its seed')
expect_equal(ss_again$confirmation$power, ss$confirmation$power,
             info = 'the confirmation is reproducible from its seed')

# ---- sample_size: confirm_B = 0 skips confirmation --------------------

ss0 <- sample_size(target = 0.80, n_grid = n_grid, confirm_B = 0L,
                   design_fn = mk_design, dgm = gd, beta = b,
                   estimand = es, analyse = anc, B = 100L,
                   intercept = 20, seed = 1L)
expect_true(is.null(ss0$confirmation),
            info = 'confirm_B = 0 skips the confirmation run')
expect_equal(ss0$n_per_arm, ss$n_per_arm,
             info = 'skipping confirmation does not change the answer')

# ---- sample_size: target not reached on the grid ---------------------

expect_error(
  sample_size(target = 0.90, n_grid = c(10L, 20L), confirm_B = 0L,
              design_fn = mk_design, dgm = gd, beta = b_null,
              estimand = es0, analyse = anc, B = 50L,
              intercept = 20, seed = 1L),
  pattern = "Extend `n_grid` upward",
  info = 'unreachable target is an error naming the remedy')

# ---- sample_size: target already exceeded at the smallest size -------

expect_warning(
  sample_size(target = 0.30, n_grid = c(40L, 80L), confirm_B = 0L,
              design_fn = mk_design, dgm = gd, beta = b,
              estimand = es, analyse = anc, B = 50L,
              intercept = 20, seed = 1L),
  pattern = "upper bound",
  info = 'target below the whole curve warns that n is an upper bound')

ss_ub <- suppressWarnings(
  sample_size(target = 0.30, n_grid = c(40L, 80L), confirm_B = 0L,
              design_fn = mk_design, dgm = gd, beta = b,
              estimand = es, analyse = anc, B = 50L,
              intercept = 20, seed = 1L))
expect_equal(ss_ub$n_per_arm, 40L,
             info = 'the upper bound returned is the smallest grid size')

# ---- sample_size: more than one method --------------------------------

expect_error(
  sample_size(target = 0.80, n_grid = c(10L, 20L), confirm_B = 0L,
              design_fn = mk_design, dgm = gd, beta = b, estimand = es,
              analyse = list(
                final = function(z) fit_ancova(z, reference = "placebo"),
                mid = function(z) fit_ancova(z, visit_time = 6,
                                             reference = "placebo")),
              B = 20L, intercept = 20, seed = 1L),
  pattern = "expects a single method",
  info = 'sample_size refuses to invert a multi-method curve')

# ---- sample_size: argument validation --------------------------------

expect_error(
  sample_size(target = 1.2, n_grid = c(10L, 20L), confirm_B = 0L,
              design_fn = mk_design, dgm = gd, beta = b, estimand = es,
              analyse = anc, B = 20L, intercept = 20, seed = 1L),
  info = 'target outside (0, 1) is rejected')
expect_error(
  sample_size(target = 0.80, n_grid = 10L, confirm_B = 0L,
              design_fn = mk_design, dgm = gd, beta = b, estimand = es,
              analyse = anc, B = 20L, intercept = 20, seed = 1L),
  info = 'a grid of one point cannot bracket the answer')
