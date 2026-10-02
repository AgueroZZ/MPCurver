test_that("advanced controls set RW3 and fixed precision across fitting paths", {
  set.seed(612)
  X <- matrix(rnorm(200), 40, 5)
  for (M in list(1L, 2L, "auto")) {
    fit <- fit_mpcurve(X, num_bins = 6, intrinsic_dim = M, max_iter = 2,
      init_control = mpcurve_init_control(max_intrinsic_dim = 2),
      control = mpcurve_control(rw_order = 3, lambda_init = 5,
                               fix_lambda = TRUE, anneal_steps = 1))
    expect_equal(fit$fit$rw_q, 3L)
    precision <- if (fit$model_intrinsic_dim == 1L) fit$fit$lambda_vec else fit$lambda_mat
    expect_length(as.numeric(precision), ncol(X) * fit$model_intrinsic_dim)
    expect_true(all(precision == 5))
    continued <- do_mpcurve(fit, max_iter = 2, tol = 0)
    precision <- if (fit$model_intrinsic_dim == 1L) continued$fit$lambda_vec else continued$lambda_mat
    expect_length(as.numeric(precision), ncol(X) * fit$model_intrinsic_dim)
    expect_true(all(precision == 5))
  }
  estimated <- fit_mpcurve(X, num_bins = 6, max_iter = 2,
                          control = mpcurve_control(lambda_init = 5))
  expect_false(all(estimated$fit$lambda_vec == 5))
})

test_that("bounds and a shared tolerance reach both fitting engines", {
  set.seed(613)
  X <- matrix(rnorm(200), 40, 5)
  for (M in c(1L, 2L)) {
    fit <- fit_mpcurve(X, num_bins = 6, intrinsic_dim = M, max_iter = 2, tol = 0.02,
      control = mpcurve_control(lambda_bounds = c(2, 3),
                               sigma2_bounds = c(0.7, 0.8), anneal_steps = 1))
    precision <- if (M == 1L) fit$fit$lambda_vec else fit$lambda_mat
    expect_true(all(precision >= 2 & precision <= 3))
    expect_true(all(fit$params$sigma2 >= 0.7 & fit$params$sigma2 <= 0.8))
    expect_equal(if (M == 1L) fit$fit$control$tol else fit$fit$control$tol_outer, 0.02)
  }
})

test_that("initialization helper settings reach all grouping paths", {
  set.seed(614)
  X <- matrix(rnorm(200), 40, 5)
  helper <- PCA_ordering
  received <- list()
  local_mocked_bindings(PCA_ordering = function(X, center = TRUE, scale = FALSE,
                                               scale01 = TRUE, component = 1L) {
    received[[length(received) + 1L]] <<- scale
    helper(X, center = center, scale = scale, scale01 = scale01, component = component)
  }, .package = "MPCurver")
  for (M in list(1L, 2L, "auto")) {
    received <- list()
    fit <- fit_mpcurve(X, num_bins = 6, intrinsic_dim = M, max_iter = 0,
      init_control = mpcurve_init_control(method_args = list(scale = TRUE),
                                         max_intrinsic_dim = 2),
      control = mpcurve_control(anneal_steps = 1))
    expect_true(length(received) >= 1)
    expect_true(all(unlist(received)))
  }
  received <- list()
  fit_mpcurve(X, num_bins = 6, intrinsic_dim = "auto", max_iter = 0,
    init_control = mpcurve_init_control(method_args = list(scale = TRUE),
                                       min_cluster_size = 100))
  expect_identical(unlist(received), TRUE)
})

test_that("invalid control settings fail before initialization", {
  X <- matrix(0, 20, 4)
  expect_error(mpcurve_control(lambda_bounds = c(2, 1)), "ordered")
  expect_error(mpcurve_control(sigma2_bounds = c(0, 1)), "positive")
  expect_error(mpcurve_control(anneal_sweeps = 1.5), "integer")
  expect_error(mpcurve_control(fix_lambda = NA), "TRUE or FALSE")
  expect_error(fit_mpcurve(X, control = list(unknown = 1)), "Unknown control")
  expect_error(fit_mpcurve(X, init_control = list(unknown = 1)), "Unknown init_control")
  for (M in list(1L, 2L, "auto")) {
    expect_error(fit_mpcurve(X, intrinsic_dim = M,
      init_control = mpcurve_init_control(method_args = list(k = 3))),
      "Unknown initialization argument")
  }
  expect_error(fit_mpcurve(X, initial_method = "isomap",
    init_control = mpcurve_init_control(pca_components = 1)), "requires initial_method")
  expect_error(fit_mpcurve(X, S = rep(1, 4),
    control = mpcurve_control(sigma2_init = 1)), "when S is non-NULL")
})
