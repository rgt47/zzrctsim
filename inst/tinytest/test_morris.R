library(tinytest)

# ADEMP conformance: RNG discipline, estimand, fitter contract,
# replicate driver, and performance measures with Monte Carlo SEs.

# ---- RNG discipline --------------------------------------------------

st <- sim_streams(5, seed = 42)
expect_equal(length(st), 5L, info = 'one stream per replicate')
expect_true(all(vapply(st, function(s) s[1] %% 100 == 7L, logical(1))),
            info = 'streams are L\'Ecuyer-CMRG (RNGkind code 7)')
expect_false(identical(st[[1]], st[[2]]),
             info = 'streams are distinct')

# Streams are reproducible from the master seed.
expect_identical(st, sim_streams(5, seed = 42),
                 info = 'same master seed gives the same streams')
expect_false(identical(st, sim_streams(5, seed = 43)),
             info = 'a different master seed gives different streams')

# with_rng_state replays a replicate exactly and restores the caller.
set.seed(999)
before <- .Random.seed
a <- with_rng_state(st[[3]], stats::rnorm(4))
b <- with_rng_state(st[[3]], stats::rnorm(4))
expect_equal(a, b, info = 'a stored stream replays exactly')
expect_identical(.Random.seed, before,
                 info = 'caller RNG state is restored')

# ---- estimand and fitter contract -----------------------------------

es <- estimand("slope difference at 15 months", -0.75)
expect_inherits(es, "estimand", info = 'estimand constructor')
expect_equal(es$value, -0.75, info = 'estimand value retained')

fr <- fit_result(estimate = 2, se = 0.5, p_value = 0.01, df = 100)
expect_true(fr$ci_lower < 2 && fr$ci_upper > 2,
            info = 'CI derived from estimate, se and df')
expect_equal(fr$ci_upper - fr$ci_lower,
             2 * stats::qt(0.975, 100) * 0.5, tolerance = 1e-12,
             info = 'CI width uses the t reference when df is finite')
nf <- null_fit()
expect_false(nf$converged, info = 'null_fit is non-converged')
expect_true(is.na(nf$estimate), info = 'null_fit carries NA')

# ---- driver, end to end ----------------------------------------------

s <- trial_schedule(treatment = 4, interval = 3)
arm <- factor(rep(c("placebo", "active"), each = 60),
              levels = c("placebo", "active"))
d <- runin_design(s, arm, reference = "placebo")
gd <- dgm_conditional(G = diag(c(9, 0.04)), sigma2 = 4)

# True estimand: change-from-baseline contrast at the final visit.
tmax <- max(s$time)
b_true <- c(x_slope = 0.5, x_hinge = -0.1, x_trt_active = -0.25)
theta <- b_true[["x_trt_active"]] * tmax

gen <- function() {
  z <- generate_outcomes(d, gd, beta = b_true, intercept = 20)
  z$arm <- z$arm
  z
}
res <- run_simulation(
  B = 200,
  generate = gen,
  analyse = list(ancova = function(z) fit_ancova(z, reference = "placebo")),
  estimand = estimand("final-visit contrast", theta),
  seed = 7L
)

expect_inherits(res, "sim_results", info = 'driver returns sim_results')
expect_equal(nrow(res), 200L, info = 'one row per replicate per method')
expect_true(all(res$converged), info = 'all replicates converged')
expect_equal(length(attr(res, "streams")), 200L,
             info = 'RNG state stored for every replicate')
expect_equal(attr(res, "estimand")$value, theta,
             info = 'estimand carried on the results')

# The stored stream reproduces its replicate exactly.
rep1 <- with_rng_state(attr(res, "streams")[[1]], gen())
rep1_fit <- fit_ancova(rep1, reference = "placebo")
expect_equal(rep1_fit$estimate, res$estimate[res$rep == 1],
             tolerance = 1e-10,
             info = 'a stored stream reproduces its replicate exactly')

# Re-running the whole study with the same seed is deterministic.
res2 <- run_simulation(200, gen,
                       list(ancova = function(z)
                         fit_ancova(z, reference = "placebo")),
                       estimand("final-visit contrast", theta), seed = 7L)
expect_equal(res$estimate, res2$estimate, tolerance = 1e-12,
             info = 'the study is reproducible from its seed')

# ---- performance measures --------------------------------------------

perf <- compute_performance(res)
need <- c("n_converged", "bias", "emp_se", "mod_se",
          "rel_error_mod_se", "mse", "coverage",
          "bias_elim_coverage", "rejection")
expect_true(all(need %in% perf$measure),
            info = 'Morris Table 6 measures are all present')

get <- function(m, col = "estimate") {
  perf[[col]][perf$measure == m]
}
expect_true(all(!is.na(perf$mcse[perf$measure %in%
                                   setdiff(need, "n_converged")])),
            info = 'every measure carries a Monte Carlo SE')

# The fitter is unbiased for this estimand.
expect_true(abs(get("bias")) < 3 * get("bias", "mcse"),
            info = 'ANCOVA is unbiased within Monte Carlo error')

# Model-based and empirical SE agree when the model is correct.
expect_true(abs(get("rel_error_mod_se")) < 10,
            info = 'model SE within 10% of empirical SE')

# Coverage is nominal.
expect_true(abs(get("coverage") - 0.95) < 3 * get("coverage", "mcse"),
            info = 'coverage is nominal within Monte Carlo error')

# MCSE formulas have the right order of magnitude.
B <- get("n_converged")
expect_equal(get("bias", "mcse"), get("emp_se") / sqrt(B),
             tolerance = 1e-8,
             info = 'MCSE(bias) = empirical SE / sqrt(B)')
expect_equal(get("coverage", "mcse"),
             sqrt(get("coverage") * (1 - get("coverage")) / B),
             tolerance = 1e-8,
             info = 'MCSE(coverage) is the binomial form')

# ---- null DGM: rejection rate is the type I error rate ---------------

b_null <- c(x_slope = 0.5, x_hinge = -0.1, x_trt_active = 0)
res0 <- run_simulation(
  B = 400,
  generate = function() generate_outcomes(d, gd, beta = b_null,
                                          intercept = 20),
  analyse = list(ancova = function(z)
    fit_ancova(z, reference = "placebo")),
  estimand = estimand("null contrast", 0), seed = 11L
)
p0 <- compute_performance(res0)
t1 <- p0$estimate[p0$measure == "rejection"]
t1_mcse <- p0$mcse[p0$measure == "rejection"]
expect_true(abs(t1 - 0.05) < 3 * t1_mcse,
            info = 'type I error is nominal under the null DGM')

# ---- errors are trapped, not propagated ------------------------------

bad <- run_simulation(
  B = 10,
  generate = function() stop("deliberate failure"),
  analyse = list(m = function(z) fit_ancova(z)),
  estimand = estimand("x", 0), seed = 1L
)
expect_equal(sum(bad$converged), 0L,
             info = 'failing replicates recorded as non-converged')
expect_equal(length(attr(bad, "errors")), 10L,
             info = 'error messages retained')

# ---- common random numbers -------------------------------------------

# Two methods on the same data must see identical data sets, so their
# estimates are paired.
res_cmp <- run_simulation(
  B = 100, generate = gen,
  analyse = list(
    full = function(z) fit_ancova(z, reference = "placebo"),
    mid = function(z) fit_ancova(z, visit_time = 6,
                                 reference = "placebo")),
  estimand = estimand("final-visit contrast", theta), seed = 3L
)
e_full <- res_cmp$estimate[res_cmp$method == "full"]
e_mid <- res_cmp$estimate[res_cmp$method == "mid"]
expect_equal(length(e_full), 100L, info = 'both methods run on all reps')
expect_true(stats::cor(e_full, e_mid) > 0.3,
            info = paste('methods sharing data are positively',
                         'correlated: common random numbers in effect'))

# ---- nsim justification ----------------------------------------------

expect_equal(nsim_for_mcse("proportion", mcse = 0.005, p = 0.05),
             ceiling(0.05 * 0.95 / 0.005^2),
             info = 'nsim for a proportion')
expect_equal(nsim_for_mcse("proportion", mcse = 0.016, p = 0.5), 977L,
             info = 'worst-case proportion at MCSE 0.016')
