# Power curves and sample size determination by simulation.
#
# Sample size is obtained by root-finding on a simulated power curve
# rather than from a closed-form expression, so it accommodates the
# design features that have no analytic power formula: run-in and
# common-close phases, staggered accrual with a common close-out,
# per-arm differential dropout, reference-based post-discontinuation
# trajectories, and analysis models whose variance is not available in
# closed form.
#
# Closed-form calculation remains the right first step. `longpower`
# gives a starting value in a fraction of the time; the role of
# simulation is to check it under the design actually planned.

#' Simulated power at a given sample size
#'
#' @param n_per_arm Integer. Subjects per arm.
#' @param design_fn Function of `n_per_arm` returning a design data
#'   frame, as from [runin_design()].
#' @param dgm A `dgm`.
#' @param beta Named fixed-effect vector, as for [generate_outcomes()].
#' @param analyze A fitter, or named list of fitters.
#' @param estimand An [estimand()].
#' @param B Integer. Replicates.
#' @param alpha Numeric. Two-sided significance level.
#' @param intercept Numeric. The overall mean at time zero, passed to
#'   [generate_outcomes()]. It shifts every arm equally and so does not
#'   affect power for a contrast, but it does set the scale on which a
#'   response-dependent `dropout` hazard is centered, so it should be
#'   given a realistic value whenever `dropout` is used.
#' @param dropout Optional list of arguments passed to
#'   [dropout_mask()], applied to each replicate.
#' @param seed Integer. Master seed.
#' @return A data frame with one row per method in `analyze` and
#'   columns
#'   \describe{
#'     \item{n_per_arm}{the sample size evaluated, repeated down the
#'       rows}
#'     \item{method}{the fitter name; `"method"` when `analyze` was a
#'       bare function rather than a named list}
#'     \item{power}{the proportion of replicates with `p < alpha`, that
#'       is the `rejection` measure from [compute_performance()]. It is
#'       power only under an alternative; under the null it is the
#'       type I error rate.}
#'     \item{mcse}{the Monte Carlo standard error of `power`, the
#'       binomial `sqrt(power * (1 - power) / m)` over the `m`
#'       replicates that converged with a non-missing p-value}
#'   }
#' @examples
#' design_fn <- function(n) {
#'   a <- factor(rep(c("placebo", "active"), each = n),
#'               levels = c("placebo", "active"))
#'   runin_design(trial_schedule(treatment = 4, interval = 3), a,
#'                reference = "placebo")
#' }
#' g <- dgm_conditional(G = diag(c(9, 0.04)), sigma2 = 4)
#' bt <- c(x_slope = 0.5, x_trt_active = -0.25)
#' es <- estimand("difference in change at month 12", -3)
#'
#' # B is far too small for a real study; see `nsim_for_mcse()`.
#' sim_power(n_per_arm = 30, design_fn = design_fn, dgm = g,
#'           beta = bt, estimand = es, analyze = fit_ancova,
#'           B = 20, intercept = 20)
#' @export
sim_power <- function(n_per_arm, design_fn, dgm, beta,
                      analyze, estimand, B = 500L, alpha = 0.05,
                      intercept = 0, dropout = NULL,
                      seed = 20260810L) {
  # Argument order follows the ADEMP flow that `run_simulation()` uses
  # -- generate, then analyze, then the estimand they are judged
  # against -- so the two functions read the same way.
  d <- design_fn(n_per_arm)
  gen <- function() {
    z <- generate_outcomes(d, dgm, beta = beta, intercept = intercept)
    if (!is.null(dropout)) {
      mk <- do.call(dropout_mask, c(list(dat = z), dropout))
      z <- apply_mask(z, mk)
    }
    z
  }
  res <- run_simulation(B, gen, analyze, estimand, seed = seed)
  perf <- compute_performance(res, alpha = alpha)
  keep <- perf$measure == "rejection"
  data.frame(n_per_arm = n_per_arm,
             method = perf$method[keep],
             power = perf$estimate[keep],
             mcse = perf$mcse[keep],
             stringsAsFactors = FALSE)
}

#' Simulated power curve
#'
#' @param n_grid Integer vector of per-arm sample sizes.
#' @param ... Passed to [sim_power()].
#' @return An object of class `c("power_curve", "data.frame")`: the
#'   [sim_power()] results row-bound over `n_grid`, so
#'   `length(n_grid) * length(analyze)` rows with the same four
#'   columns
#'   \describe{
#'     \item{n_per_arm}{the per-arm sample size for the row}
#'     \item{method}{the fitter name}
#'     \item{power}{simulated rejection proportion at that size}
#'     \item{mcse}{its Monte Carlo standard error}
#'   }
#'   Every point on the curve is simulated under the same `seed`, so
#'   the sizes share random numbers and the curve is smoother than
#'   independent runs would give. That correlation is deliberate but it
#'   means the points are not independent estimates.
#' @examples
#' design_fn <- function(n) {
#'   a <- factor(rep(c("placebo", "active"), each = n),
#'               levels = c("placebo", "active"))
#'   runin_design(trial_schedule(treatment = 4, interval = 3), a,
#'                reference = "placebo")
#' }
#' g <- dgm_conditional(G = diag(c(9, 0.04)), sigma2 = 4)
#'
#' power_curve(n_grid = c(30, 60), design_fn = design_fn, dgm = g,
#'             beta = c(x_slope = 0.5, x_trt_active = -0.25),
#'             estimand = estimand("change difference", -3),
#'             analyze = fit_ancova, B = 20)
#' @export
power_curve <- function(n_grid, ...) {
  out <- do.call(rbind, lapply(n_grid, function(n) {
    sim_power(n_per_arm = n, ...)
  }))
  rownames(out) <- NULL
  class(out) <- c("power_curve", "data.frame")
  out
}

#' @export
print.power_curve <- function(x, ...) {
  cat("<power_curve>\n")
  print(as.data.frame(x), row.names = FALSE)
  invisible(x)
}

#' Sample size for a target power, by simulation
#'
#' Searches for the smallest per-arm sample size attaining `power`.
#' A monotone interpolation is fitted through a simulated power curve
#' and inverted, then the answer is confirmed by a longer run at the
#' selected size.
#'
#' @param target Numeric. Desired power.
#' @param n_grid Integer vector. Sample sizes at which to evaluate the
#'   curve. Should bracket the answer.
#' @param confirm_B Integer. Replicates for the confirmation run at the
#'   selected size. `0` skips confirmation.
#' @param ... Passed to [sim_power()].
#' @return An object of class `sample_size`, a list with
#'   \describe{
#'     \item{n_per_arm}{integer. The selected per-arm sample size:
#'       `ceiling()` of the interpolated size, or the smallest value on
#'       `n_grid` when `target` was already exceeded there, in which
#'       case a warning states that it is only an upper bound.}
#'     \item{target}{the requested power, as supplied}
#'     \item{curve}{the `power_curve` the selection was interpolated
#'       from}
#'     \item{confirmation}{a one-row [sim_power()] data frame at
#'       `n_per_arm`, run with `confirm_B` replicates on an independent
#'       RNG stream, or `NULL` when `confirm_B = 0`. Its `power` and
#'       `mcse` are the figures to quote.}
#'   }
#' @details
#' The returned sample size inherits Monte Carlo error from the curve.
#' The confirmation run reports the achieved power with its MCSE at the
#' chosen size, which is the number to quote. If the confirmation MCSE
#' is large relative to the distance from `target`, raise `B`; see
#' [nsim_for_mcse()].
#'
#' The confirmation run uses a different RNG stream from the curve
#' (the curve's seed, offset by one), so that it is an independent
#' check of the selected size rather than a re-reading of the
#' replicates that selected it.
#' @export
sample_size <- function(target = 0.80, n_grid, confirm_B = 2000L, ...) {
  stopifnot(target > 0, target < 1, length(n_grid) >= 2)
  dots <- list(...)
  # Checked before simulating, not after: inverting a curve is only
  # meaningful for one method, and running the whole grid first would
  # make the user pay for the entire study before being told the call
  # was invalid.
  if (is.list(dots$analyze) && length(dots$analyze) > 1L) {
    stop("`sample_size()` expects a single method; `analyze` has ",
         length(dots$analyze), ": ",
         paste(names(dots$analyze), collapse = ", "),
         ". Call it once per method.")
  }
  curve <- do.call(power_curve, c(list(n_grid = n_grid), dots))
  meths <- unique(curve$method)
  if (length(meths) > 1L) {
    stop("`sample_size()` expects a single method; got: ",
         paste(meths, collapse = ", "))
  }

  if (max(curve$power) < target) {
    stop("Target power ", target, " not reached on the grid; ",
         "the largest size gave ", format(max(curve$power), digits = 3),
         ". Extend `n_grid` upward.")
  }
  if (min(curve$power) > target) {
    # `approx()` would return NA here, since the target lies outside
    # the range of observed powers. Report the smallest grid value and
    # say plainly that it is an upper bound, rather than interpolating
    # off the end of the curve.
    smallest <- which.min(curve$n_per_arm)
    n_sel <- as.integer(curve$n_per_arm[smallest])
    # The quoted power must come from the row actually being returned,
    # not from `min(curve$power)`: on a curve made non-monotone by
    # Monte Carlo error those are different rows, and the warning
    # would then attribute one size's power to another.
    warning("Target power ", target, " is already exceeded at the ",
            "smallest size on the grid (", n_sel, ", power ",
            format(curve$power[smallest], digits = 3), "). The ",
            "returned size is an upper bound; extend `n_grid` ",
            "downward for a sharper answer.")
  } else {
    # A simulated power curve is not guaranteed monotone: Monte Carlo
    # error can put a dip in it, and `approx()` will interpolate
    # straight through that dip and return a size from a region where
    # power decreases in n. Say so rather than answering silently.
    by_n <- curve[order(curve$n_per_arm), ]
    if (is.unsorted(by_n$power)) {
      warning("The simulated power curve is not monotone in ",
              "`n_per_arm` (powers ",
              paste(format(by_n$power, digits = 3), collapse = ", "),
              " at sizes ",
              paste(by_n$n_per_arm, collapse = ", "),
              "). This is Monte Carlo error; the interpolated size is ",
              "correspondingly uncertain. Raise `B`, or widen the ",
              "spacing of `n_grid`.")
    }

    o <- order(curve$power)
    n_hat <- stats::approx(x = curve$power[o], y = curve$n_per_arm[o],
                           xout = target, ties = "ordered")$y
    n_sel <- as.integer(ceiling(n_hat))

    # `approx()` with tied powers returns the LAST tie, so a curve
    # that saturates (power 1.000 repeated across sizes, which is
    # routine) yields the largest size rather than the smallest.
    # The documented contract is the smallest size attaining `target`,
    # so take that directly whenever the grid already contains one at
    # or below the interpolated value.
    attains <- by_n$n_per_arm[by_n$power >= target]
    if (length(attains)) {
      n_sel <- as.integer(min(n_sel, min(attains)))
    }
  }
  if (!is.finite(n_sel) || n_sel < 1L) {
    stop("Could not determine a sample size from the supplied grid; ",
         "the simulated power curve was ",
         paste(format(curve$power, digits = 3), collapse = ", "), ".")
  }

  confirmation <- NULL
  if (confirm_B > 0L) {
    # `B` may already be present in `...` for the curve; the
    # confirmation run overrides it.
    dots$B <- confirm_B
    # The confirmation must also use a different RNG stream. Reusing
    # the curve's seed would draw the same replicates that produced
    # the interpolated `n_sel`, so the "confirmed" power would share
    # its Monte Carlo error with the selection rather than providing
    # an independent check of it.
    curve_seed <- if (is.null(dots$seed)) {
      eval(formals(sim_power)$seed)
    } else {
      dots$seed
    }
    dots$seed <- as.integer(curve_seed) + 1L
    confirmation <- do.call(sim_power,
                            c(list(n_per_arm = n_sel), dots))
  }

  structure(list(n_per_arm = n_sel, target = target,
                 curve = curve, confirmation = confirmation),
            class = "sample_size")
}

#' @export
print.sample_size <- function(x, ...) {
  cat("<sample_size> target power", x$target, "\n")
  cat("  n per arm:", x$n_per_arm,
      " (total", 2 * x$n_per_arm, "for two arms)\n")
  if (!is.null(x$confirmation)) {
    cat("  confirmed power:",
        format(x$confirmation$power, digits = 4),
        " (MCSE", format(x$confirmation$mcse, digits = 3), ")\n")
  }
  cat("\n")
  print(as.data.frame(x$curve), row.names = FALSE)
  invisible(x)
}
