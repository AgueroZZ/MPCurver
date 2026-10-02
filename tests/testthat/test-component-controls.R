test_that("continuation inherits settings and consistently controls precision updates", {
  set.seed(790)
  X <- matrix(rnorm(160), 40, 4)
  for (M in 1:2) {
    original <- fit_mpcurve(X, num_bins = 5, intrinsic_dim = M, max_iter = 2,
      tol = 0, position_prior = "fixed", control = mpcurve_control(
        lambda_init = 2, fix_lambda = TRUE, lambda_bounds = c(1, 10),
        lambda_sd_prior_rate = 0.7, convergence = "relative", anneal_steps = 1))
    precision <- function(fit) if (M == 1) fit$fit$lambda_vec else fit$lambda_mat
    inherited <- do_mpcurve(original, max_iter = 2)
    expect_true(all(precision(inherited) == 2))
    expect_equal(inherited$fit$control$lambda_sd_prior_rate, 0.7)
    expect_equal(inherited$fit$control$convergence, "relative")
    expect_equal(inherited$fit$control$position_prior, "fixed")
    expect_equal(inherited$data, X)
    expect_equal(head(inherited$elbo_trace, length(original$elbo_trace)), original$elbo_trace)

    reset <- expect_no_warning(do_mpcurve(original, max_iter = 2, control =
      mpcurve_continue_control(lambda_init = 5, fix_lambda = TRUE)))
    expect_true(all(precision(reset) == 5))
    clipped <- expect_no_warning(do_mpcurve(reset, max_iter = 2,
      control = list(lambda_bounds = c(1, 3))))
    expect_true(all(precision(clipped) == 3))
    free <- expect_no_warning(do_mpcurve(reset, max_iter = 2, control =
      mpcurve_continue_control(fix_lambda = FALSE, lambda_sd_prior_rate = 0,
        lambda_bounds = c(1, 3), sigma2_bounds = c(0.2, 0.4))))
    expect_false(free$fit$control$fix_lambda)
    expect_null(free$fit$control$lambda_sd_prior_rate)
    expect_true(all(precision(free) >= 1 & precision(free) <= 3))
    expect_true(all(free$params$sigma2 >= 0.2 & free$params$sigma2 <= 0.4))
    expect_length(free$continuation_history, 2)
    expect_equal(free$continuation_history[[2]]$control$fix_lambda, FALSE)
    expect_true(all(precision(original) == 2))

    # A reset alone must keep empirical-Bayes updates enabled on both paths.
    estimated <- fit_mpcurve(X, num_bins = 5, intrinsic_dim = M, max_iter = 2,
      tol = 0, control = mpcurve_control(anneal_steps = 1))
    updated <- expect_no_warning(do_mpcurve(estimated, max_iter = 2,
      control = list(lambda_init = 5)))
    expect_false(updated$fit$control$fix_lambda)
    expect_false(all(precision(updated) == 5))
  }
})

test_that("continuation preserves known measurement errors and initialization provenance", {
  set.seed(791)
  X <- matrix(rnorm(160), 40, 4)
  for (M in 1:2) {
    fit <- fit_mpcurve(X, S = rep(0.3, 4), num_bins = 5, intrinsic_dim = M,
      max_iter = 1, init_control = mpcurve_init_control(method_args = list(scale = TRUE)),
      control = mpcurve_control(anneal_steps = 1))
    more <- do_mpcurve(fit, max_iter = 2)
    expect_equal(more$measurement_sd, fit$measurement_sd)
    expect_null(more$params$sigma2)
    expect_equal(more$fit$control$method_args, fit$fit$control$method_args)
  }
})

test_that("retired and irrelevant controls fail before fitting", {
  X <- matrix(0, 20, 4)
  expect_error(do_mpcurve(structure(list(), class = "mpcurve"), max_iter = 1.5), "integer")
  expect_error(mpcurve_continue_control(lambda_bounds = c(5, 1)), "ordered")
  expect_error(mpcurve_continue_control(fix_lambda = NA), "TRUE or FALSE")
  for (M in list(1L, 2L, "auto")) {
    expect_error(fit_mpcurve(X, intrinsic_dim = M, initial_method = "isomap",
      init_control = list(method_args = list(control = list(eps = 1e-8)))),
      "Unknown isomap control")
  }
  expect_error(isomap_ordering(X, control = list(num_landmarks = 1.5)), "integer")
  expect_error(fiedler_ordering(X, control = list(weight = "unknown")), "arg")
  expect_error(PCA_ordering(X, component = 1.5), "integer")
  expect_error(simulate_mpcurve(control = list(K = 5)), "Unknown control")
  expect_error(simulate_dual_trajectory(control = list(latent_positions = matrix(0, 10, 2))),
               "Unknown control")
  expect_error(simulate_intrinsic_trajectories(n = 2.5), "integer")
  expect_error(simulate_intrinsic_trajectories(d_signal = c(2, 1.5)), "integers")
  expect_error(simulate_two_order_gp_dataset(control = list(permute_cols = NA)), "TRUE or FALSE")
})

test_that("simulation controls set interpretable generating quantities", {
  sim <- simulate_mpcurve(n = 40, d = 3, num_bins = 5, seed = 793, control = list(
    rw_order = 3, lambda_range = c(5, 5), noise_sd_range = c(0.2, 0.2),
    position_weights = c(1, 0, 0, 0, 0)))
  expect_equal(sim$rw_q, 3L)
  expect_equal(sim$lambda_vec, rep(5, 3))
  expect_equal(sim$sigma2, rep(0.04, 3))
  expect_true(all(sim$z == 1L))
  latent <- cbind(seq(0, 1, length.out = 40), seq(1, 0, length.out = 40))
  sim <- simulate_intrinsic_trajectories(n = 40, d_signal = c(2, 2), d_noise = 0,
    noise_sd = 0, seed = 793, trajectory_family = "linear", control = list(
      latent_positions = latent, signal_range = c(1, 1),
      linear_slope_range = c(2, 2), intercept_sd = 0))
  expect_equal(sim$latent_positions, latent)
  for (j in seq_len(ncol(sim$X))) {
    m <- match(sim$true_assign[j], sim$ordering_labels)
    expect_equal(abs(diff(range(sim$X[, j]))), 2)
    expect_equal(abs(cor(sim$X[, j], latent[, m])), 1)
  }
})
