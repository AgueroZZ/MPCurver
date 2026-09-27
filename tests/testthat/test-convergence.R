convergence_data <- function() {
  set.seed(326)
  t <- seq(-1, 1, length.out = 40)
  cbind(t, t^2, sin(3 * t), rev(t)^2, cos(3 * t), -t) +
    matrix(rnorm(240, sd = 0.2), 40, 6)
}

test_that("normalized changes use N times D and ignore objective offsets", {
  change <- .mpcurve_elbo_change
  expect_equal(change(2, -30, 40, 6, "normalized", 1), 2 / 240)
  expect_equal(change(2, 1e9, 40, 6, "normalized", 1), 2 / 240)
  expect_equal(change(2, -30, 40, 6, "relative", 1), 2 / 31)
  expect_equal(change(2, -30, 40, 6, "relative", 1e-12), 2 / (30 + 1e-12))
})

test_that("single-ordering stopping matches the first eligible trace increment", {
  X <- convergence_data()
  for (S in list(NULL, rep(0.2, ncol(X)))) {
    baseline <- fit_mpcurve(X, K = 6, iter = 60, tol = 0, S = S)
    expect_identical(baseline$fit$control$convergence, "normalized")
    trace <- baseline$elbo_trace
    delta <- diff(trace)
    previous <- head(trace, -1L)
    for (rule in c("normalized", "relative")) {
      metric <- abs(delta) / if (rule == "normalized") length(X) else (abs(previous) + 1)
      threshold <- median(metric[delta >= 0])
      expected <- which(delta >= 0 & metric < threshold)[1]
      expect_true(is.finite(expected))
      fit <- fit_mpcurve(X, K = 6, iter = 60, tol = threshold,
                        convergence = rule, S = S)
      expect_true(fit$converged)
      expect_equal(fit$elbo_trace, trace[seq_len(expected + 1L)], tolerance = 1e-10)
      expect_identical(fit$fit$control$convergence, rule)
    }
  }
})

test_that("partition stopping uses the full matrix and only post-annealing increments", {
  X <- convergence_data()
  for (S in list(NULL, rep(0.2, ncol(X)))) {
    args <- list(X = X, K = 6, intrinsic_dim = 2, n_outer = 2,
                 inner_iter = 1, T_start = 2, T_end = 1,
                 max_converge_iter = 40, S = S)
    baseline <- do.call(fit_mpcurve, c(args, list(tol_outer = 0)))
    trace <- baseline$elbo_trace
    # History: initial state, two annealing steps, then exact T = 1 sweeps.
    delta <- diff(trace)[-(1:2)]
    previous <- head(trace, -1L)[-(1:2)]
    eligible <- delta >= -1e-8 * (abs(previous) + 1)
    for (rule in c("normalized", "relative")) {
      metric <- abs(delta) / if (rule == "normalized") length(X) else (abs(previous) + 1e-12)
      threshold <- median(metric[eligible])
      expected <- which(eligible & metric < threshold)[1]
      fit <- do.call(fit_mpcurve, c(args, list(tol_outer = threshold, convergence = rule)))
      expect_true(fit$converged)
      expect_equal(fit$elbo_trace, trace[seq_len(expected + 3L)], tolerance = 1e-10)
      expect_identical(fit$fit$control$convergence, rule)
      expect_equal(fit$fit$convergence_info$n_observations, length(X))
    }
  }
})

test_that("continuation inherits, overrides, and migrates stopping rules", {
  X <- convergence_data()
  for (M in 1:2) {
    for (rule in c("normalized", "relative")) {
      fit <- fit_mpcurve(X, K = 6, intrinsic_dim = M, iter = 3,
                         n_outer = 1, T_start = 1, T_end = 1,
                         tol = 0, tol_outer = 0, convergence = rule)
      continued <- do_mpcurve(fit, iter = 3)
      expect_identical(continued$fit$control$convergence, rule)
      expect_equal(head(continued$elbo_trace, length(fit$elbo_trace)), fit$elbo_trace)
      override <- do_mpcurve(fit, iter = 3, convergence = "normalized")
      expect_identical(override$fit$control$convergence, "normalized")
      fit$fit$control$convergence <- NULL
      legacy <- do_mpcurve(fit, iter = 3)
      expect_identical(legacy$fit$control$convergence, "relative")
      expect_equal(legacy$elbo_trace, continued$elbo_trace, tolerance = 1e-10)
      expect_error(do_mpcurve(fit, convergence = "invalid"), "arg")
    }
  }
})

test_that("public fitting validates the rule and shares the default tolerance", {
  X <- convergence_data()
  expect_error(fit_mpcurve(X, convergence = "invalid"), "arg")
  expect_equal(formals(fit_mpcurve)$tol, 1e-6)
  expect_equal(formals(fit_mpcurve)$tol_outer, 1e-6)
  fit <- soft_two_trajectory_cavi(X, K = 6, n_outer = 1,
                                 max_converge_iter = 1, verbose = FALSE)
  expect_identical(fit$control$convergence, "normalized")
})

test_that("continued fits stop at the chosen normalized or relative increment", {
  X <- convergence_data()
  for (M in 1:2) {
    start <- fit_mpcurve(X, K = 6, intrinsic_dim = M, iter = 2,
                        n_outer = 1, T_start = 1, T_end = 1,
                        tol = 0, tol_outer = 0)
    baseline <- do_mpcurve(start, iter = 40)
    old_length <- length(start$elbo_trace)
    tail_trace <- baseline$elbo_trace[old_length:length(baseline$elbo_trace)]
    delta <- diff(tail_trace)
    previous <- head(tail_trace, -1L)
    eligible <- if (M == 1) delta >= 0 else delta >= -1e-8 * (abs(previous) + 1)
    for (rule in c("normalized", "relative")) {
      offset <- if (M == 1) 1 else 1e-12
      metric <- abs(delta) / if (rule == "normalized") length(X) else (abs(previous) + offset)
      threshold <- median(metric[eligible])
      expected <- which(eligible & metric < threshold)[1]
      fit <- do_mpcurve(start, iter = 40, tol = threshold, tol_outer = threshold,
                        convergence = rule)
      expect_true(fit$converged)
      expect_equal(fit$elbo_trace, head(baseline$elbo_trace, old_length + expected),
                   tolerance = 1e-10)
    }
  }
})

test_that("dimension selection passes the rule to every candidate", {
  X <- convergence_data()
  original <- fit_mpcurve
  rules <- character()
  testthat::local_mocked_bindings(fit_mpcurve = function(...) {
    fit <- original(...)
    rules <<- c(rules, fit$fit$control$convergence)
    fit
  })
  for (rule in c("normalized", "relative")) {
    rules <- character()
    fit <- select_mpcurve_dimension(X, max_intrinsic_dim = 2, K = 6,
                                    iter = 2, n_outer = 1,
                                    T_start = 1, T_end = 1, convergence = rule)
    expect_equal(length(rules), 3L)
    expect_true(all(rules == rule))
    expect_identical(fit$fit$control$convergence, rule)
  }
})
