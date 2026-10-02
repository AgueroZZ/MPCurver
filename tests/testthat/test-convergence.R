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

test_that("default partition fitting starts at T = 1 without extra sweeps", {
  X <- convergence_data()
  expect_identical(mpcurve_control()$anneal_steps, 0L)
  expect_error(mpcurve_control(anneal_steps = -1), "integer >= 0")
  expect_error(mpcurve_control(anneal_steps = 0.5), "integer >= 0")
  for (S in list(NULL, rep(0.2, ncol(X)))) {
    initial <- fit_mpcurve(X, S = S, num_bins = 6, intrinsic_dim = 2,
                          max_iter = 0, tol = 0)
    expect_equal(initial$fit$iter, 0)
    expect_equal(initial$temperature_history, 1)
    expect_length(initial$objective_history, 1)
    expect_equal(initial$fit$n_anneal, 1)
    fit <- fit_mpcurve(X, S = S, num_bins = 6, intrinsic_dim = 2,
                      max_iter = 3, tol = 0)
    expect_equal(fit$fit$iter, 3)
    expect_equal(fit$temperature_history, rep(1, 4))
    expect_equal(fit$objective_history[1], initial$objective_history)
    expect_identical(fit$fit$control$n_outer, 0L)
    continued <- do_mpcurve(fit, max_iter = 2, tol = 0)
    uninterrupted <- fit_mpcurve(X, S = S, num_bins = 6, intrinsic_dim = 2,
                                max_iter = 5, tol = 0)
    expect_equal(continued$temperature_history, rep(1, 6))
    expect_equal(continued$objective_history, uninterrupted$objective_history,
                 tolerance = 1e-10)
  }
})

test_that("disabled annealing ignores its temperature and sweep settings", {
  X <- convergence_data()
  baseline <- fit_mpcurve(X, num_bins = 6, intrinsic_dim = 2, max_iter = 2, tol = 0)
  inactive <- fit_mpcurve(X, num_bins = 6, intrinsic_dim = 2, max_iter = 2, tol = 0,
    control = mpcurve_control(anneal_steps = 0, anneal_start = 20, anneal_sweeps = 5))
  expect_identical(inactive$objective_history, baseline$objective_history)
  expect_identical(inactive$partition$pi_weights, baseline$partition$pi_weights)
  expect_identical(inactive$fit$gamma, baseline$fit$gamma)
})

test_that("explicit annealing precedes the temperature-1 fitting budget", {
  fit <- fit_mpcurve(convergence_data(), num_bins = 6, intrinsic_dim = 2,
    max_iter = 3, tol = 0,
    control = mpcurve_control(anneal_steps = 3, anneal_start = 4, anneal_sweeps = 2))
  expect_equal(fit$temperature_history, c(4, 4, 4, 2, 2, 1, 1, 1, 1, 1))
  expect_equal(fit$fit$iter, 9)
  expect_length(fit$objective_history, 10)
  expect_equal(fit$fit$n_anneal, 7)
})

test_that("single-ordering stopping matches the first eligible trace increment", {
  X <- convergence_data()
  for (S in list(NULL, rep(0.2, ncol(X)))) {
    baseline <- fit_mpcurve(
    X,
    num_bins = 6,
    max_iter = 60,
    tol = 0,
    S = S
  )
    expect_identical(baseline$fit$control$convergence, "normalized")
    trace <- baseline$elbo_trace
    delta <- diff(trace)
    previous <- head(trace, -1L)
    for (rule in c("normalized", "relative")) {
      metric <- abs(delta) / if (rule == "normalized") length(X) else (abs(previous) + 1)
      threshold <- median(metric[delta >= 0])
      expected <- which(delta >= 0 & metric < threshold)[1]
      expect_true(is.finite(expected))
      fit <- fit_mpcurve(
    X,
    num_bins = 6,
    max_iter = 60,
    tol = threshold,
    S = S,
    control = mpcurve_control(convergence = rule)
  )
      expect_true(fit$converged)
      expect_equal(fit$elbo_trace, trace[seq_len(expected + 1L)], tolerance = 1e-10)
      expect_identical(fit$fit$control$convergence, rule)
    }
  }
})

test_that("partition stopping uses the full matrix and only post-annealing increments", {
  X <- convergence_data()
  for (S in list(NULL, rep(0.2, ncol(X)))) {
    args <- list(
    X = X,
    num_bins = 6,
    intrinsic_dim = 2,
    max_iter = 40,
    S = S,
    control = mpcurve_control(anneal_steps = 2, anneal_sweeps = 1, anneal_start = 2)
  )
    baseline <- do.call(fit_mpcurve, c(args, list(tol = 0)))
    trace <- baseline$elbo_trace
    # History: initial state, two annealing steps, then exact T = 1 sweeps.
    delta <- diff(trace)[-(1:2)]
    previous <- head(trace, -1L)[-(1:2)]
    eligible <- delta >= -1e-8 * (abs(previous) + 1)
    for (rule in c("normalized", "relative")) {
      metric <- abs(delta) / if (rule == "normalized") length(X) else (abs(previous) + 1e-12)
      threshold <- median(metric[eligible])
      expected <- which(eligible & metric < threshold)[1]
      args$control$convergence <- rule
      fit <- do.call(fit_mpcurve, c(args, list(tol = threshold)))
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
      fit <- fit_mpcurve(
    X,
    num_bins = 6,
    intrinsic_dim = M,
    max_iter = 3,
    tol = 0,
    control = mpcurve_control(anneal_steps = 1, anneal_start = 1, convergence = rule)
  )
      continued <- do_mpcurve(fit, max_iter = 3)
      expect_identical(continued$fit$control$convergence, rule)
      expect_equal(head(continued$elbo_trace, length(fit$elbo_trace)), fit$elbo_trace)
      override <- do_mpcurve(fit, max_iter = 3, control = mpcurve_continue_control(convergence = "normalized"))
      expect_identical(override$fit$control$convergence, "normalized")
      fit$fit$control$convergence <- NULL
      legacy <- do_mpcurve(fit, max_iter = 3)
      expect_identical(legacy$fit$control$convergence, "relative")
      expect_equal(legacy$elbo_trace, continued$elbo_trace, tolerance = 1e-10)
      expect_error(do_mpcurve(fit, control = mpcurve_continue_control(convergence = "invalid")), "arg")
    }
  }
})

test_that("public fitting validates the rule and shares the default tolerance", {
  X <- convergence_data()
  expect_error(fit_mpcurve(
    X,
    control = mpcurve_control(convergence = "invalid")
  ), "arg")
  expect_equal(formals(fit_mpcurve)$tol, 1e-6)
  expect_false("tol_outer" %in% names(formals(fit_mpcurve)))
  fit <- soft_two_trajectory_cavi(X, K = 6, n_outer = 1,
                                 max_converge_iter = 1, verbose = FALSE)
  expect_identical(fit$control$convergence, "normalized")
})

test_that("continued fits stop at the chosen normalized or relative increment", {
  X <- convergence_data()
  for (M in 1:2) {
    start <- fit_mpcurve(
    X,
    num_bins = 6,
    intrinsic_dim = M,
    max_iter = 2,
    tol = 0,
    control = mpcurve_control(anneal_steps = 1, anneal_start = 1)
  )
    baseline <- do_mpcurve(start, max_iter = 40)
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
      fit <- do_mpcurve(
        start,
        max_iter = 40,
        tol = threshold,
        control = mpcurve_continue_control(convergence = rule)
      )
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
    fit <- select_mpcurve_dimension(
  X,
  max_intrinsic_dim = 2,
    num_bins = 6,
    max_iter = 2,
    control = mpcurve_control(anneal_steps = 1, anneal_start = 1, convergence = rule)
  )
    expect_equal(length(rules), 2L)
    expect_true(all(rules == rule))
    expect_identical(fit$fit$control$convergence, rule)
  }
})
