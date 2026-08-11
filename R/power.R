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
#' @param estimand An [estimand()].
#' @param analyse A fitter, or named list of fitters.
#' @param B Integer. Replicates.
#' @param alpha Numeric. Two-sided significance level.
#' @param intercept Numeric.
#' @param dropout Optional list of arguments passed to
#'   [dropout_mask()], applied to each replicate.
#' @param seed Integer. Master seed.
#' @return A data frame with `n_per_arm`, `method`, `power` and
#'   `mcse`.
#' @export
sim_power <- function(n_per_arm, design_fn, dgm, beta, estimand,
                      analyse, B = 500L, alpha = 0.05,
                      intercept = 0, dropout = NULL,
                      seed = 20260810L) {
  d <- design_fn(n_per_arm)
  gen <- function() {
    z <- generate_outcomes(d, dgm, beta = beta, intercept = intercept)
    if (!is.null(dropout)) {
      mk <- do.call(dropout_mask, c(list(dat = z), dropout))
      z <- apply_mask(z, mk)
    }
    z
  }
  res <- run_simulation(B, gen, analyse, estimand, seed = seed)
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
#' @return A data frame of `n_per_arm`, `method`, `power`, `mcse`,
#'   with class `power_curve`.
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
#' @return An object of class `sample_size`: a list with `n_per_arm`,
#'   `target`, the `curve`, and the `confirmation` row.
#' @details
#' The returned sample size inherits Monte Carlo error from the curve.
#' The confirmation run reports the achieved power with its MCSE at the
#' chosen size, which is the number to quote. If the confirmation MCSE
#' is large relative to the distance from `target`, raise `B`; see
#' [nsim_for_mcse()].
#' @export
sample_size <- function(target = 0.80, n_grid, confirm_B = 2000L, ...) {
  stopifnot(target > 0, target < 1, length(n_grid) >= 2)
  dots <- list(...)
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
    n_sel <- as.integer(min(curve$n_per_arm))
    warning("Target power ", target, " is already exceeded at the ",
            "smallest size on the grid (", n_sel, ", power ",
            format(min(curve$power), digits = 3), "). The returned ",
            "size is an upper bound; extend `n_grid` downward for a ",
            "sharper answer.")
  } else {
    # Monotone interpolation of n as a function of power.
    o <- order(curve$power)
    n_hat <- stats::approx(x = curve$power[o], y = curve$n_per_arm[o],
                           xout = target, ties = "ordered")$y
    n_sel <- as.integer(ceiling(n_hat))
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
