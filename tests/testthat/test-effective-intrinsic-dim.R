test_that("effective ordering count ignores numerical prior dust", {
  count_effective <- getFromNamespace(
    ".mpcurve_effective_intrinsic_dim", "MPCurver"
  )
  priors <- list(partition = list(
    mode = "adaptive",
    omega = c(A = 0.5, B = 0.5, C = 8.89318162514244e-323, D = 0)
  ))

  expect_identical(count_effective(priors, 4L, 1e-12), 2L)
  expect_identical(count_effective(priors, 4L, 0), 3L)
  priors$partition$mode <- "fixed"
  expect_identical(count_effective(priors, 4L, 1e-12), 4L)
})

test_that("fitted and continued models report effective intrinsic dimension", {
  set.seed(117)
  X <- matrix(rnorm(144), nrow = 24, ncol = 6)
  common <- list(
    X = X,
    intrinsic_dim = 2L,
    K = 5L,
    method = c("PCA", "PCA"),
    partition_init = "ordering_methods",
    T_start = 1,
    T_end = 1,
    n_outer = 1L,
    max_converge_iter = 1L,
    effective_weight_tol = 0.05,
    verbose = FALSE
  )

  adaptive <- suppressWarnings(do.call(fit_mpcurve, common))
  omega <- fitted_prior(adaptive, type = "partition")$omega
  expect_identical(adaptive$effective_intrinsic_dim,
                   as.integer(sum(omega > 0.05)))
  expect_identical(adaptive$effective_weight_tol, 0.05)
  expect_identical(adaptive$fit$control$effective_weight_tol, 0.05)
  expect_identical(adaptive$active_intrinsic_dim, 2L)
  expect_identical(summary(adaptive)$effective_intrinsic_dim,
                   adaptive$effective_intrinsic_dim)
  expect_match(paste(capture.output(print(adaptive)), collapse = "\n"),
               "Effective dim", fixed = TRUE)

  continued <- suppressWarnings(do_mpcurve(adaptive, iter = 1L))
  omega_continued <- fitted_prior(continued, type = "partition")$omega
  expect_identical(continued$effective_weight_tol, 0.05)
  expect_identical(continued$effective_intrinsic_dim,
                   as.integer(sum(omega_continued > 0.05)))

  fixed <- suppressWarnings(do.call(
    fit_mpcurve, c(common, list(partition_prior = "fixed"))
  ))
  expect_identical(fixed$effective_intrinsic_dim, 2L)
  expect_identical(summary(fixed)$effective_intrinsic_dim, 2L)

  single <- suppressWarnings(fit_mpcurve(X, K = 5L, iter = 1L))
  expect_identical(single$effective_intrinsic_dim, 1L)
  expect_identical(summary(single)$effective_intrinsic_dim, 1L)
})

test_that("effective ordering threshold is validated before fitting", {
  X <- matrix(rnorm(24), nrow = 6, ncol = 4)
  expect_error(
    fit_mpcurve(X, intrinsic_dim = 2L, effective_weight_tol = -1),
    "effective_weight_tol"
  )
  expect_error(
    fit_mpcurve(X, intrinsic_dim = 2L, effective_weight_tol = 0.5),
    "effective_weight_tol"
  )
  expect_error(
    fit_mpcurve(X, intrinsic_dim = 2L, effective_weight_tol = NA_real_),
    "effective_weight_tol"
  )
})
