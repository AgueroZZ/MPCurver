test_that("structural conditional trajectory update is unweighted and analytical", {
  set.seed(101)
  X <- matrix(rnorm(48), nrow = 12, ncol = 4)
  R <- matrix(runif(60), nrow = 12, ncol = 5)
  R <- R / rowSums(R)
  sigma2 <- c(0.7, 1.1, 1.4, 0.9)
  lambda <- c(0.5, 1, 2, 3)
  Q <- make_random_walk_precision(K = 5, d = 1, q = 2, ridge = 0)

  got <- .structural_partition_q_u(
    X = X,
    R = R,
    sigma2 = sigma2,
    measurement_sd = NULL,
    lambda_vec = lambda,
    Q_K = Q,
    ordering_label = "A"
  )

  Nk <- colSums(R)
  for (j in seq_len(ncol(X))) {
    A <- lambda[j] * Q + diag(Nk / sigma2[j], 5)
    b <- as.numeric(crossprod(X[, j], R)) / sigma2[j]
    expect_equal(got$S_list[[j]], solve(A), tolerance = 1e-9)
    expect_equal(got$m_mat[j, ], as.numeric(solve(A, b)), tolerance = 1e-9)
  }

  single_ordering <- .cavi_update_q_u_weighted(
    X = X,
    R = R,
    sigma2 = sigma2,
    lambda_vec = lambda,
    Q_K = Q,
    measurement_sd = NULL,
    feature_weights = rep(1, ncol(X)),
    rw_q = 2
  )
  expect_equal(got$m_mat, single_ordering$m_mat, tolerance = 1e-10)
  expect_equal(got$S_list, single_ordering$S_list, tolerance = 1e-10)

  # q(U_j | Z_j = m) has no w_jm argument, so exact zero assignment mass
  # cannot alter its precision or posterior moments.
  got_again <- .structural_partition_q_u(
    X = X,
    R = R,
    sigma2 = sigma2,
    measurement_sd = NULL,
    lambda_vec = lambda,
    Q_K = Q,
    ordering_label = "A"
  )
  expect_equal(got_again$m_mat, got$m_mat, tolerance = 0)
  expect_equal(got_again$S_list, got$S_list, tolerance = 0)
})

test_that("structural known-noise trajectory update uses full likelihood precision", {
  set.seed(102)
  X <- matrix(rnorm(36), nrow = 12, ncol = 3)
  R <- matrix(runif(48), nrow = 12, ncol = 4)
  R <- R / rowSums(R)
  S <- matrix(runif(36, 0.5, 1.5), nrow = 12, ncol = 3)
  lambda <- c(0.5, 1, 2)
  Q <- make_random_walk_precision(K = 4, d = 1, q = 2, ridge = 0)

  got <- .structural_partition_q_u(
    X = X,
    R = R,
    sigma2 = NULL,
    measurement_sd = S,
    lambda_vec = lambda,
    Q_K = Q,
    ordering_label = "A"
  )
  inv_v <- 1 / S[, 2]^2
  A <- lambda[2] * Q + diag(as.numeric(crossprod(inv_v, R)), 4)
  b <- as.numeric(crossprod(inv_v * X[, 2], R))
  expect_equal(got$S_list[[2]], solve(A), tolerance = 1e-9)
  expect_equal(got$m_mat[2, ], as.numeric(solve(A, b)), tolerance = 1e-9)
})

test_that("shared sigma update pools conditional residuals with q(Z) weights", {
  set.seed(103)
  X <- matrix(rnorm(40), nrow = 10, ncol = 4)
  R1 <- matrix(runif(30), 10, 3)
  R2 <- matrix(runif(30), 10, 3)
  R1 <- R1 / rowSums(R1)
  R2 <- R2 / rowSums(R2)
  Q <- make_random_walk_precision(K = 3, d = 1, q = 1, ridge = 0)
  q1 <- .structural_partition_q_u(X, R1, rep(1, 4), NULL, rep(1, 4), Q, "A")
  q2 <- .structural_partition_q_u(X, R2, rep(1, 4), NULL, rep(1, 4), Q, "B")
  weights <- cbind(c(0, 0.2, 0.6, 1), c(1, 0.8, 0.4, 0))

  got <- .structural_partition_update_sigma2(
    X,
    gamma = list(R1, R2),
    q_u = list(q1, q2),
    weights = weights,
    sigma_min = 1e-12,
    sigma_max = 1e12
  )
  E1 <- .structural_partition_expected_residuals(X, R1, q1)
  E2 <- .structural_partition_expected_residuals(X, R2, q2)
  expected <- (weights[, 1] * E1 + weights[, 2] * E2) / nrow(X)
  expect_equal(got$sigma2, expected, tolerance = 1e-10)
})

test_that("structural local blocks, assignment update, and objective match formulas", {
  set.seed(107)
  X <- matrix(rnorm(24), nrow = 8, ncol = 3)
  R1 <- matrix(runif(24), 8, 3)
  R2 <- matrix(runif(24), 8, 3)
  R1 <- R1 / rowSums(R1)
  R2 <- R2 / rowSums(R2)
  sigma2 <- c(0.8, 1.2, 1.5)
  lambda_mat <- cbind(c(0.7, 1.1, 1.4), c(1.3, 0.9, 0.6))
  Q <- make_random_walk_precision(K = 3, d = 1, q = 1, ridge = 0)
  q1 <- .structural_partition_q_u(X, R1, sigma2, NULL, lambda_mat[, 1], Q, "A")
  q2 <- .structural_partition_q_u(X, R2, sigma2, NULL, lambda_mat[, 2], Q, "B")
  local <- .structural_partition_local_blocks(
    X = X,
    gamma = list(R1, R2),
    q_u = list(q1, q2),
    sigma2 = sigma2,
    measurement_sd = NULL,
    lambda_mat = lambda_mat,
    Q_K = Q,
    rw_q = 1,
    lambda_sd_prior_rate = NULL
  )

  E11 <- .structural_partition_expected_residuals(X, R1, q1)[1]
  meta <- .rw_precision_metadata(Q, rw_q = 1)
  expected_B11 <- -0.5 * nrow(X) * log(2 * pi * sigma2[1]) -
    0.5 * E11 / sigma2[1] +
    0.5 * (meta$rank * (log(lambda_mat[1, 1]) - log(2 * pi)) +
             meta$logdet - lambda_mat[1, 1] * q1$eq_quad[1]) +
    0.5 * (ncol(R1) * log(2 * pi * exp(1)) + q1$logdetS[1])
  expect_equal(local[1, 1], expected_B11, tolerance = 1e-9)

  weights0 <- matrix(0.5, nrow = ncol(X), ncol = 2)
  control <- list(
    partition_prior = "fixed",
    partition_prior_init = NULL,
    assignment_prior = NULL,
    ordering_alpha = NA_real_
  )
  T_now <- 1.5
  update <- .structural_partition_update_weights(local, weights0, T_now, control)
  expected_w <- exp(local / T_now - matrixStats::rowLogSumExps(local / T_now))
  expect_equal(update$weights, expected_w, tolerance = 1e-10)

  pi1 <- colMeans(R1)
  pi2 <- colMeans(R2)
  objective <- .structural_partition_objective(
    gamma = list(R1, R2),
    position_pi = list(pi1, pi2),
    weights = update$weights,
    local_blocks = local,
    assignment_info = update$assignment_info
  )
  cell <- .cavi_cell_terms_from_state(R1, pi1) +
    .cavi_cell_terms_from_state(R2, pi2)
  entropy_z <- -sum(update$weights * log(update$weights))
  expected_objective <- cell + sum(update$weights * local) -
    ncol(X) * log(2) + T_now * entropy_z
  expect_equal(objective$objective, expected_objective, tolerance = 1e-9)
})

test_that("fixed-M structural fits use canonical shared state with intrinsic ridge zero", {
  set.seed(104)
  X <- matrix(rnorm(144), nrow = 24, ncol = 6)

  fit2 <- suppressWarnings(fit_mpcurve(
    X,
    intrinsic_dim = 2,
    K = 5,
    partition_init = "ordering_methods",
    method = c("PCA", "PCA"),
    ridge = 0,
    T_start = 1,
    T_end = 1,
    n_outer = 1,
    max_converge_iter = 3,
    verbose = FALSE
  ))
  expect_s3_class(fit2, "mpcurve")
  expect_identical(fit2$variational_family, "structured")
  expect_equal(fit2$intrinsic_dim, 2L)
  expect_equal(length(fit2$params$sigma2), ncol(X))
  expect_equal(dim(fit2$lambda_mat), c(ncol(X), 2L))
  expect_equal(dim(fit2$partition$pi_weights), c(ncol(X), 2L))
  expect_equal(rowSums(fit2$partition$pi_weights), rep(1, ncol(X)), tolerance = 1e-10)
  expect_null(fit2$fit$active_orderings)
  expect_null(fit2$fit$active_feature_pairs)
  expect_null(fit2$fit$effective_pi_weights)
  expect_null(fit2$fit$frozen_orderings)
  expect_true(all(vapply(fit2$fits, function(child) {
    isTRUE(all.equal(child$params$sigma2, fit2$params$sigma2, tolerance = 0))
  }, logical(1))))

  t1 <- fit2$temperature_history == 1
  expect_true(all(diff(fit2$objective_history[t1]) >= -1e-6))

  fit3 <- suppressWarnings(fit_mpcurve(
    X,
    intrinsic_dim = 3,
    K = 5,
    partition_init = "ordering_methods",
    method = rep("PCA", 3),
    ridge = 0,
    T_start = 1,
    T_end = 1,
    n_outer = 1,
    max_converge_iter = 2,
    verbose = FALSE
  ))
  expect_equal(fit3$intrinsic_dim, 3L)
  expect_equal(dim(fit3$partition$pi_weights), c(ncol(X), 3L))
  expect_equal(length(fit3$conditional_posterior$mean), 3L)
})

test_that("known S and structural continuation preserve canonical schema", {
  set.seed(105)
  X <- matrix(rnorm(120), nrow = 24, ncol = 5)
  S_vec <- runif(5, 0.6, 1.2)
  fit <- suppressWarnings(fit_mpcurve(
    X,
    S = S_vec,
    intrinsic_dim = 2,
    K = 5,
    partition_init = "ordering_methods",
    method = c("PCA", "PCA"),
    ridge = 0,
    T_start = 1,
    T_end = 1,
    n_outer = 1,
    max_converge_iter = 1,
    verbose = FALSE
  ))
  expect_null(fit$params$sigma2)
  before <- length(fit$objective_history)
  continued <- suppressWarnings(do_mpcurve(fit, iter = 2, S = S_vec))
  expect_gt(length(continued$objective_history), before)
  expect_null(continued$params$sigma2)
  expect_identical(continued$variational_family, "structured")

  S_mat <- matrix(runif(length(X), 0.5, 1.5), nrow(X), ncol(X))
  fit_mat <- suppressWarnings(fit_mpcurve(
    X,
    S = S_mat,
    intrinsic_dim = 2,
    K = 5,
    partition_init = "ordering_methods",
    method = c("PCA", "PCA"),
    ridge = 0,
    T_start = 1,
    T_end = 1,
    n_outer = 1,
    max_converge_iter = 1,
    verbose = FALSE
  ))
  expect_null(fit_mat$params$sigma2)
  expect_equal(fit_mat$measurement_sd, S_mat)
})

test_that("transition controls warn and greedy is disabled", {
  set.seed(106)
  X <- matrix(rnorm(80), nrow = 20, ncol = 4)
  expect_warning(
    fit_mpcurve(
      X,
      intrinsic_dim = 2,
      K = 4,
      partition_init = "ordering_methods",
      method = c("PCA", "PCA"),
      n_outer = 1,
      max_converge_iter = 0,
      freeze_feature = TRUE,
      verbose = FALSE
    ),
    "deprecated and ignored"
  )
  expect_error(
    fit_mpcurve(X, intrinsic_dim = 3, greedy = "forward"),
    "temporarily unavailable"
  )
})
