missing_objective <- function(fit) {
  if (fit$intrinsic_dim == 1L) fit$fit$elbo_trace else fit$objective_history
}

missing_example <- function() {
  set.seed(4218)
  t <- seq(-1, 1, length.out = 30)
  X <- cbind(t, t^2, sin(3 * t), cos(3 * t)) + matrix(rnorm(120, sd = 0.15), 30, 4)
  X[matrix(runif(120) < 0.2, 30, 4)] <- NA_real_
  X[1, ] <- c(-1, 1, -0.1, -1)
  X[rowSums(!is.na(X)) == 0, 1] <- 0
  dimnames(X) <- list(paste0("s", 1:30), paste0("g", 1:4))
  X
}

# This oracle evaluates scalar Gaussian densities, without masked helper calls
# or the implementation's expanded quadratic forms.
scalar_missing_blocks <- function(X, R, q_u, variance, lambda, Q) {
  d <- ncol(X)
  K <- ncol(R)
  likelihood <- numeric(d)
  logits <- matrix(0, nrow(X), K)
  residual <- numeric(d)
  for (j in seq_len(d)) for (i in seq_len(nrow(X))) {
    if (is.na(X[i, j])) next
    for (k in seq_len(K)) {
      expected <- dnorm(X[i, j], q_u$m_mat[j, k], sqrt(variance[i, j]), log = TRUE) -
        q_u$S_list[[j]][k, k] / (2 * variance[i, j])
      likelihood[j] <- likelihood[j] + R[i, k] * expected
      logits[i, k] <- logits[i, k] + expected
      residual[j] <- residual[j] + R[i, k] *
        ((X[i, j] - q_u$m_mat[j, k])^2 + q_u$S_list[[j]][k, k])
    }
  }
  eigenvalues <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
  positive <- eigenvalues[eigenvalues > 1e-8 * max(eigenvalues)]
  local <- likelihood
  for (j in seq_len(d)) {
    mu <- q_u$m_mat[j, ]
    V <- q_u$S_list[[j]]
    local[j] <- local[j] + 0.5 * (length(positive) * log(lambda[j] / (2 * pi)) +
      sum(log(positive)) - lambda[j] * (sum(mu * (Q %*% mu)) + sum(diag(Q %*% V)))) +
      0.5 * (K * log(2 * pi * exp(1)) + as.numeric(determinant(V, logarithm = TRUE)$modulus))
  }
  list(likelihood = likelihood, logits = logits, residual = residual, local = local)
}

test_that("observed Gaussian systems, density constants, and residuals match scalar calculations", {
  X <- missing_example()[1:8, 1:3]
  X[-1, 3] <- NA_real_
  K <- 4L
  set.seed(420)
  R <- matrix(runif(nrow(X) * K, 0.2, 1), nrow(X), K)
  R <- R / rowSums(R)
  lambda <- c(0.7, 1.3, 2.1)
  sigma2 <- c(0.4, 0.8, 1.2)
  for (ridge in c(0, 0.2)) for (known in c(FALSE, TRUE)) {
    Q <- make_random_walk_precision(K, d = 1, q = 2, ridge = ridge)
    variance <- if (known) matrix(seq(0.2, 1.3, length.out = length(X)), nrow(X)) else
      matrix(rep(sigma2, each = nrow(X)), nrow(X))
    S <- if (known) sqrt(variance) else NULL
    if (known) S[is.na(X)] <- NA_real_
    got <- .structural_partition_q_u(X, R, if (known) NULL else sigma2,
      S, lambda, Q)
    for (j in seq_len(ncol(X))) {
      precision <- lambda[j] * Q
      rhs <- numeric(K)
      for (i in seq_len(nrow(X))) if (!is.na(X[i, j])) for (k in seq_len(K)) {
        precision[k, k] <- precision[k, k] + R[i, k] / variance[i, j]
        rhs[k] <- rhs[k] + R[i, k] * X[i, j] / variance[i, j]
      }
      expect_equal(got$S_list[[j]], solve(precision), tolerance = 1e-10)
      expect_equal(got$m_mat[j, ], drop(solve(precision, rhs)), tolerance = 1e-10)
    }
    oracle <- scalar_missing_blocks(X, R, got, variance, lambda, Q)
    scores <- .cavi_feature_scores_from_state(X, R, got,
      if (known) NULL else sigma2, lambda, S, Q, include_prior = FALSE)
    expect_equal(scores, oracle$likelihood, tolerance = 1e-11)
    expect_equal(.structural_partition_expected_residuals(X, R, got),
      oracle$residual, tolerance = 1e-11)
    blocks <- .structural_partition_local_blocks(X, list(R), list(got),
      if (known) NULL else sigma2, S, matrix(lambda, ncol = 1), Q, 2, NULL)
    expect_equal(drop(blocks), oracle$local, tolerance = 1e-10)
    pi_vec <- c(0.1, 0.2, 0.3, 0.4)
    expected_logits <- sweep(oracle$logits, 2, log(pi_vec), "+")
    expected_R <- exp(expected_logits - matrixStats::rowLogSumExps(expected_logits))
    updated_R <- .structural_partition_update_gamma(X, got, pi_vec,
      if (known) NULL else sigma2, S, rep(1, ncol(X)))
    expect_equal(updated_R, expected_R, tolerance = 1e-10)
    cell <- sum(R * rep(log(pi_vec), each = nrow(X))) - sum(R * log(R))
    expect_equal(.cavi_observed_elbo(X, R, pi_vec, got,
      if (known) NULL else sigma2, S, lambda, Q, 2, NULL),
      cell + sum(oracle$local), tolerance = 1e-10)
    # The Gaussian coordinate maximizes its complete conditional local block.
    perturbed <- got
    perturbed$m_mat[1, ] <- perturbed$m_mat[1, ] + c(0.1, -0.2, 0.05, 0.1)
    changed <- scalar_missing_blocks(X, R, perturbed, variance, lambda, Q)
    expect_lt(changed$local[1], oracle$local[1])
  }
})

test_that("shared variance uses feature observation counts and assignment-weighted residuals", {
  X <- missing_example()[1:8, 1:3]
  X[-1, 3] <- NA_real_
  set.seed(4201)
  Rs <- lapply(1:2, function(m) {
    R <- matrix(runif(nrow(X) * 4, 0.2, 1), nrow(X), 4)
    R / rowSums(R)
  })
  Q <- make_random_walk_precision(4, d = 1, q = 2, ridge = 0.1)
  variance <- matrix(0.8, nrow(X), ncol(X))
  states <- lapply(Rs, function(R) .structural_partition_q_u(
    X, R, rep(0.8, 3), NULL, rep(1, 3), Q))
  residual <- vapply(1:2, function(m) scalar_missing_blocks(
    X, Rs[[m]], states[[m]], variance, rep(1, 3), Q)$residual, numeric(3))
  weights <- rbind(c(0.1, 0.9), c(0.7, 0.3), c(1, 0))
  expected <- rowSums(weights * residual) / colSums(!is.na(X))
  actual <- .structural_partition_update_sigma2(X, Rs, states, weights, 1e-10, 1e10)
  expect_equal(unname(actual$sigma2), unname(expected), tolerance = 1e-11)
  # Differentiate the independent scalar variance objective at its update.
  for (j in 1:3) {
    derivative <- -sum(!is.na(X[, j])) / (2 * actual$sigma2[j]) +
      sum(weights[j, ] * residual[j, ]) / (2 * actual$sigma2[j]^2)
    expect_equal(unname(derivative), 0, tolerance = 1e-10)
  }
})

test_that("missing-data kernels reduce to existing complete-data kernels", {
  X <- .mpcurve_initialization_data(missing_example())
  R <- matrix(1 / 5, nrow(X), 5)
  Q <- make_random_walk_precision(5, d = 1, q = 2, ridge = 0.1)
  sigma2 <- c(0.2, 0.4, 0.6, 0.8)
  lambda <- c(0.7, 1, 1.3, 2)
  for (S in list(NULL, rep(0.5, 4), matrix(seq(0.3, 0.8, length.out = length(X)), nrow(X)))) {
    variance <- .cavi_resolve_measurement_sd(S, X)$var_mat
    complete <- .structural_partition_q_u(X, R, sigma2, S, lambda, Q)
    observed <- .cavi_observed_q_u(X, R, sigma2, lambda, Q, variance)
    expect_equal(observed, complete, tolerance = 1e-10)
    expect_equal(.cavi_observed_residuals(X, R, observed),
      .structural_partition_expected_residuals(X, R, complete), tolerance = 1e-10)
  }
})

test_that("masked kernels match scalar sums under large offsets and feature weights", {
  # A common offset of 1e4 would cancel about nine digits in uncentered
  # expanded squares; the kernels center each feature at its observed mean.
  X <- missing_example()[1:12, ] + 1e4
  K <- 5L
  set.seed(4202)
  R <- matrix(runif(nrow(X) * K, 0.1, 1), nrow(X), K)
  R <- R / rowSums(R)
  Q <- make_random_walk_precision(K, d = 1, q = 2, ridge = 0.1)
  lambda <- c(0.5, 1, 2, 4)
  weights <- c(0.2, 1, 0.6, 0.9)
  for (known in c(FALSE, TRUE)) {
    variance <- if (known) matrix(seq(0.3, 1.1, length.out = length(X)), nrow(X)) else
      matrix(rep(c(0.3, 0.5, 0.7, 0.9), each = nrow(X)), nrow(X))
    sigma2 <- if (known) NULL else variance[1, ]
    S <- if (known) sqrt(variance) else NULL
    if (known) S[is.na(X)] <- NA_real_
    q_u <- .structural_partition_q_u(X, R, sigma2, S, lambda, Q)
    for (j in seq_len(ncol(X))) {
      precision <- lambda[j] * Q
      rhs <- numeric(K)
      for (i in which(!is.na(X[, j]))) {
        precision <- precision + diag(R[i, ] / variance[i, j])
        rhs <- rhs + R[i, ] * X[i, j] / variance[i, j]
      }
      expect_equal(q_u$m_mat[j, ], drop(solve(precision, rhs)), tolerance = 1e-10)
    }
    oracle <- scalar_missing_blocks(X, R, q_u, variance, lambda, Q)
    expect_equal(.cavi_observed_residuals(X, R, q_u), oracle$residual, tolerance = 1e-10)
    expect_equal(.cavi_feature_scores_from_state(X, R, q_u, sigma2, lambda, S, Q,
      include_prior = FALSE), oracle$likelihood, tolerance = 1e-10)
    weighted <- matrix(0, nrow(X), K)
    for (i in seq_len(nrow(X))) for (j in which(!is.na(X[i, ]))) for (k in seq_len(K)) {
      weighted[i, k] <- weighted[i, k] + weights[j] *
        (dnorm(X[i, j], q_u$m_mat[j, k], sqrt(variance[i, j]), log = TRUE) -
           q_u$S_list[[j]][k, k] / (2 * variance[i, j]))
    }
    var_mat <- if (known) replace(variance, is.na(X), NA_real_) else NULL
    expect_equal(.cavi_observed_logdensity(X, q_u, sigma2, var_mat, weights),
      weighted, tolerance = 1e-10)
  }
})

test_that("intrinsic precision requires observed support on its null space", {
  X <- matrix(c(1, NA, NA, NA, 2, 3, 4, 5), 4, 2)
  R <- matrix(0, 4, 4)
  R[, 1] <- 1
  Q <- make_random_walk_precision(4, d = 1, q = 2, ridge = 0)
  expect_error(.cavi_observed_q_u(X, R, c(1, 1), c(1, 1), Q),
               "feature 1 is singular: its 1 observed value.*occupy 1 position bin")
  R[1, ] <- c(0.5, 0, 0, 0.5)
  expect_true(all(is.finite(.cavi_observed_q_u(X, R, c(1, 1), c(1, 1), Q)$m_mat)))
  # A zero observed value still supplies likelihood precision.
  X[1, 1] <- 0
  expect_true(all(is.finite(.cavi_observed_q_u(X, R, c(1, 1), c(1, 1), Q)$m_mat)))
  # Mass no larger than the numerical probability floor carries no data, so
  # a bin holding only such mass leaves the system singular.
  R[1, ] <- c(1 - 3e-12, rep(1e-12, 3))
  expect_error(.cavi_observed_q_u(X, R, c(1, 1), c(1, 1), Q),
               "feature 1 is singular: its 1 observed value.*occupy 1 position bin")
  # Supported but numerically unreliable systems are rejected separately.
  R[1, ] <- c(1 - 2e-7, 1e-7, 1e-7, 0)
  expect_error(.cavi_observed_q_u(X, R, c(1, 1), c(1e6, 1), Q),
               "too ill-conditioned.*feature 1")
  expect_true(all(is.finite(.cavi_observed_q_u(
    X, R, c(1, 1), c(1, 1), Q + diag(0.1, 4))$m_mat)))
})

test_that("empty subset rows have prior-only positions and hyperpriors remain valid", {
  X <- missing_example()
  X[2, ] <- NA_real_
  R <- matrix(1 / 5, nrow(X), 5)
  pi_vec <- c(0.1, 0.2, 0.3, 0.25, 0.15)
  Q <- make_random_walk_precision(5, d = 1, q = 2, ridge = 0.1)
  q_u <- .cavi_observed_q_u(X, R, rep(1, 4), rep(1, 4), Q)
  gamma <- .structural_partition_update_gamma(
    X, q_u, pi_vec = pi_vec, sigma2 = rep(1, 4),
    measurement_sd = NULL, weights = rep(1, 4))
  expect_equal(as.numeric(gamma[2, ]), pi_vec, tolerance = 1e-10)
  X <- missing_example()
  for (M in 1:2) {
    fit <- fit_mpcurve(X, num_bins = 5, intrinsic_dim = M, max_iter = 6, tol = 0,
      control = mpcurve_control(rw_order = 3, lambda_sd_prior_rate = 0.5))
    expect_true(all(is.finite(missing_objective(fit))))
    expect_true(all(diff(missing_objective(fit)) >= -1e-7))
  }
})

test_that("partial NA inputs fit and continue under all current noise/prior modes", {
  X <- missing_example()
  S_matrix <- matrix(seq(0.2, 0.5, length.out = length(X)), nrow(X))
  S_matrix[is.na(X)] <- NA_real_
  for (M in 1:2) for (position in c("fixed", "adaptive")) {
    for (S in list(NULL, rep(0.3, ncol(X)), S_matrix)) for (ridge in c(0, 0.1)) {
      args <- list(X = X, S = S, num_bins = 5, intrinsic_dim = M,
        position_prior = position, partition_prior = position, tol = 0,
        control = mpcurve_control(ridge = ridge, lambda_init = 1))
      fit <- do.call(fit_mpcurve, c(args, list(max_iter = 4)))
      continued <- do_mpcurve(fit, max_iter = 2, tol = 0)
      whole <- do.call(fit_mpcurve, c(args, list(max_iter = 6)))
      expect_identical(fit$data, X)
      expect_equal(missing_objective(continued), missing_objective(whole), tolerance = 1e-8)
      expect_equal(fitted_trajectories(continued), fitted_trajectories(whole), tolerance = 1e-8)
      expect_true(all(is.finite(missing_objective(fit))))
      expect_true(all(diff(missing_objective(fit)) >= -1e-7))
      expect_equal(fit$observations$n_observed, sum(!is.na(X)))
      expect_true(all(is.finite(fitted_positions(fit))))
      expect_true(all(is.finite(fitted_trajectories(fit, "sd"))))
      expect_identical(dim(fitted_assignments(fit)), c(4L, as.integer(M)))
    }
  }
})

test_that("unused measurement errors cannot affect missing-data inference", {
  X <- missing_example()
  S1 <- matrix(0.3, nrow(X), ncol(X))
  S2 <- S1
  S3 <- S1
  S1[is.na(X)] <- NA_real_
  S2[is.na(X)] <- 1e5
  S3[is.na(X)] <- 0
  for (M in 1:2) {
    fits <- lapply(list(S1, S2, S3), function(S) fit_mpcurve(
      X, S = S, num_bins = 5, intrinsic_dim = M, max_iter = 3, tol = 0))
    for (j in 2:3) {
      expect_equal(missing_objective(fits[[j]]), missing_objective(fits[[1]]), tolerance = 1e-10)
      expect_equal(fitted_trajectories(fits[[j]]), fitted_trajectories(fits[[1]]), tolerance = 1e-10)
    }
  }
})

test_that("missing input contracts reject uninformative margins and invalid SDs", {
  X <- missing_example()
  Y <- X
  Y[1, ] <- NA_real_
  expect_error(fit_mpcurve(Y), "entirely missing sample")
  expect_error(select_mpcurve_dimension(Y, max_intrinsic_dim = 2), "entirely missing sample")
  Y <- X
  Y[, 1] <- NA_real_
  expect_error(fit_mpcurve(Y), "entirely missing feature")
  for (bad in c(NaN, Inf, -Inf)) {
    Y <- X
    Y[1, 1] <- bad
    expect_error(fit_mpcurve(Y), "NaN or Inf")
  }
  S <- matrix(1, nrow(X), ncol(X))
  S[1, 1] <- NA_real_
  expect_error(fit_mpcurve(X, S = S), "at observed X entries")
  expect_error(fit_mpcurve(X, S = c(NA, 1, 1, 1)), "at observed X entries")
  expect_equal(.mpcurve_elbo_change(1, 100, 30, 4, "normalized", 1,
                                   sum(!is.na(X))), 1 / sum(!is.na(X)))
})

test_that("missing inputs work through initialization, dimension selection, and plots", {
  X <- missing_example()
  for (method in c("PCA", "isomap", "fiedler", "random")) {
    set.seed(42)
    fit <- fit_mpcurve(X, initial_method = method, num_bins = 5, max_iter = 2)
    expect_true(all(is.finite(missing_objective(fit))))
    expect_identical(fit$data, X)
  }
  fit <- fit_mpcurve(X, intrinsic_dim = "auto", num_bins = 5, max_iter = 2,
    init_control = mpcurve_init_control(max_intrinsic_dim = 2))
  expect_identical(fit$data, X)
  selected <- select_mpcurve_dimension(X, max_intrinsic_dim = 2, num_bins = 5,
                                       max_iter = 2)
  expect_identical(selected$data, X)
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  for (M in 1:2) {
    fit <- fit_mpcurve(X, intrinsic_dim = M, num_bins = 5, max_iter = 2)
    expect_no_warning(plot(fit, dims = 1))
    expect_no_warning(plot(fit, dims = c(1, 2)))
    expect_no_warning(plot(fit, plot_type = "mu", dims = 1))
  }
})

test_that("precision initialization uses only bins that carry data", {
  K <- 8L
  Q <- as.matrix(make_random_walk_precision(K, d = 1, q = 2))
  mu <- c(0.2, NA, 1.1, NA, NA, 0.4, NA, 2)
  rough <- .cavi_occupied_roughness(mu, Q)
  # Brute force: smoothest completion of the empty bins.
  occupied <- !is.na(mu)
  fill <- function(u_empty) { u <- mu; u[!occupied] <- u_empty; u }
  best <- solve(Q[!occupied, !occupied], -Q[!occupied, occupied] %*% mu[occupied])
  expect_equal(rough$quad, drop(crossprod(fill(best), Q %*% fill(best))), tolerance = 1e-10)
  expect_equal(rough$rank, sum(occupied) - 2)
  lambda <- .cavi_lambda_init_from_mu_raw(rbind(mu, c(1:K)), Q, r_rank = K - 2)
  expect_equal(lambda[1], (sum(occupied) - 2) / rough$quad, tolerance = 1e-10)
  # A linear raw trajectory has no curvature on its occupied bins.
  expect_equal(lambda[2], 1e10)
  # One occupied bin leaves no curvature information: start from 1.
  expect_equal(.cavi_lambda_init_from_mu_raw(rbind(c(3, rep(NA, K - 1))), Q, K - 2), 1)
  # Raw moments leave bins without data undefined rather than zero.
  X <- matrix(rnorm(40), 20, 2)
  R <- matrix(0, 20, K); R[1:10, 2] <- 1; R[11:20, 5] <- 1
  raw <- .cavi_raw_moments_from_responsibilities(X, R)
  expect_true(all(is.na(raw$mu_raw[, -c(2, 5)])))
  expect_equal(unname(raw$mu_raw[, 2]), colMeans(X[1:10, ]))
})

test_that("singular trajectory precision stops in complete-data fits too", {
  set.seed(2)
  X <- matrix(rnorm(60), 20, 3)
  collapsed <- matrix(0, 20, 5); collapsed[, 2] <- 1
  expect_error(cavi(X, K = 5, responsibilities_init = collapsed, max_iter = 1,
                    verbose = FALSE), "every feature is singular: the data occupy 1 position bin")
  two_bins <- collapsed; two_bins[11:20, 2] <- 0; two_bins[11:20, 4] <- 1
  fit <- cavi(X, K = 5, responsibilities_init = two_bins, max_iter = 0, verbose = FALSE)
  expect_true(all(is.finite(fit$posterior$mean)))
  # A proper prior identifies every direction.
  fit <- cavi(X, K = 5, responsibilities_init = collapsed, ridge = 0.1, max_iter = 1,
              verbose = FALSE)
  expect_true(all(is.finite(fit$posterior$mean)))
})

test_that("feature similarities use pairwise-complete observations", {
  X <- missing_example()
  observed <- !is.na(X)
  for (metric in c("spearman", "pearson")) {
    S <- .compute_same_ordering_similarity(X, metric = metric)$S
    pairwise <- abs(cor(X, method = metric, use = "pairwise.complete.obs"))
    diag(pairwise) <- 1
    expect_equal(unname(S), unname(pairwise), tolerance = 1e-12)
    expect_equal(S[1, 2], abs(cor(X[observed[, 1] & observed[, 2], 1],
      X[observed[, 1] & observed[, 2], 2], method = metric)), ignore_attr = TRUE)
  }
  # Spline R^2: order the jointly observed samples by the predictor.
  df <- 3L
  S <- .compute_same_ordering_similarity(X, metric = "spline_r2", spline_r2_df = df)$S
  directional <- function(p, j) {
    rows <- which(observed[, p] & observed[, j])
    y <- X[rows[order(X[rows, p])], j]
    position <- (seq_along(rows) - 0.5) / length(rows)
    summary(lm(y ~ splines::ns(position, df = df)))$r.squared
  }
  expect_equal(S[1, 3], max(directional(1, 3), directional(3, 1)), tolerance = 1e-10,
               ignore_attr = TRUE)
  # Smooth-fit evidence uses the same samples, rescaled per observation.
  sm <- .compute_same_ordering_similarity(X, metric = "smooth_fit", K = 5L)
  rows <- which(observed[, 1] & observed[, 2])
  pair <- .cavi_fit_directional_smooth_evidence(x_cov = X[rows, 1], y = X[rows, 2], K = 5L)
  expect_equal(sm$directional_delta[1, 2], pair$delta * nrow(X) / length(rows),
               tolerance = 1e-10, ignore_attr = TRUE)
  # Values of a matrix S at missing entries cannot change the similarity.
  S_known <- matrix(0.3, nrow(X), ncol(X)); S_alt <- S_known
  S_known[!observed] <- NA_real_; S_alt[!observed] <- 50
  a <- .compute_same_ordering_similarity(X, S = S_known, metric = "smooth_fit", K = 5L)$S
  b <- .compute_same_ordering_similarity(X, S = S_alt, metric = "smooth_fit", K = 5L)$S
  expect_equal(a, b)
})

test_that("exported ordering helpers fill missing values with a warning", {
  X <- missing_example()
  filled <- .mpcurve_initialization_data(X)
  expect_warning(got <- PCA_ordering(X), "PCA_ordering\\(\\) does not accept missing values")
  expect_equal(got, PCA_ordering(filled))
  for (helper in list(fiedler_ordering, isomap_ordering, pcurve_ordering)) {
    expect_warning(helper(X), "does not accept missing values")
  }
  expect_silent(PCA_ordering(filled))
})

test_that("missing-data position updates match the complete-data limit", {
  # A missing entry equals an observed entry with infinite measurement SD; with
  # no probability floor the two paths agree to rounding.
  X <- missing_example()
  filled <- .mpcurve_initialization_data(X)
  S <- matrix(0.3, nrow(X), ncol(X)); S_na <- S; S_na[is.na(X)] <- NA_real_
  S_inf <- S; S_inf[is.na(X)] <- 1e8
  R0 <- matrix(1, nrow(X), 5); R0[cbind(seq_len(nrow(X)), rep(1:5, length.out = nrow(X)))] <- 6
  R0 <- R0 / rowSums(R0)
  a <- cavi(X, K = 5, S = S_na, responsibilities_init = R0, lambda_init = 1, max_iter = 15, tol = 0)
  b <- cavi(filled, K = 5, S = S_inf, responsibilities_init = R0, lambda_init = 1, max_iter = 15, tol = 0)
  expect_equal(a$gamma, b$gamma, tolerance = 1e-8)
  expect_equal(a$posterior$mean, b$posterior$mean, tolerance = 1e-8)
})

test_that("singular updates under missing data suggest a fixed position prior", {
  set.seed(7)
  X <- matrix(rnorm(40 * 3), 40, 3); X[sample(length(X), 30)] <- NA
  X[1, ] <- 0
  collapsed <- matrix(0, nrow(X), 5); collapsed[, 2] <- 1
  expect_error(fit_mpcurve(X, num_bins = 5, max_iter = 1,
    init_control = mpcurve_init_control(responsibilities_init = collapsed)),
    "singular.*consider position_prior = \"fixed\"")
  expect_error(fit_mpcurve(X, num_bins = 5, max_iter = 1, position_prior = "fixed",
    init_control = mpcurve_init_control(responsibilities_init = collapsed)),
    "singular")
})

test_that("similarity initialization reports one-feature clusters", {
  set.seed(11)
  t <- runif(80)
  X <- cbind(sin(2 * pi * t), cos(2 * pi * t), t, t^2, rnorm(80))
  X <- X + matrix(rnorm(length(X), sd = 0.05), nrow(X))
  colnames(X) <- paste0("g", 1:5)
  expect_message(fit_mpcurve(X, intrinsic_dim = 2, max_iter = 2),
                 "alone in ordering .*consider a smaller intrinsic_dim or another similarity_metric")
})
