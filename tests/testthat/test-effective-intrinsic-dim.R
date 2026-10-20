test_that("explicit fitted and continued models retain the specified dimension", {
  set.seed(117)
  X <- matrix(rnorm(144), nrow = 24, ncol = 6)
  common <- list(X = X, intrinsic_dim = 2L, num_bins = 5L,
                 initial_method = "PCA", max_iter = 1L,
                 init_control = mpcurve_init_control(),
                 control = mpcurve_control(anneal_start = 1, anneal_steps = 1L,
                                          effective_count_tol = 0.05), verbose = FALSE)

  adaptive <- suppressWarnings(do.call(fit_mpcurve, common))
  expect_identical(adaptive$intrinsic_dim, 2L)
  expect_identical(adaptive$effective_count_tol, 0.05)
  expect_identical(adaptive$fit$control$effective_count_tol, 0.05)
  expect_identical(summary(adaptive)$intrinsic_dim,
                   adaptive$intrinsic_dim)
  expect_match(paste(capture.output(print(adaptive)), collapse = "\n"),
               "Intrinsic dim", fixed = TRUE)

  continued <- suppressWarnings(do_mpcurve(adaptive, max_iter = 1L))
  expect_identical(continued$effective_count_tol, 0.05)
  expect_identical(continued$intrinsic_dim, 2L)

  fixed <- suppressWarnings(do.call(
    fit_mpcurve, c(common, list(partition_prior = "fixed"))
  ))
  expect_identical(fixed$intrinsic_dim, 2L)
  expect_identical(summary(fixed)$intrinsic_dim, 2L)

  single <- suppressWarnings(fit_mpcurve(
    X,
    num_bins = 5L,
    max_iter = 1L
  ))
  expect_identical(single$intrinsic_dim, 1L)
  expect_identical(summary(single)$intrinsic_dim, 1L)
})

test_that("effective ordering threshold is validated before fitting", {
  X <- matrix(rnorm(24), nrow = 6, ncol = 4)
  expect_error(
    fit_mpcurve(
    X,
    intrinsic_dim = 2L,
    control = mpcurve_control(effective_count_tol = -1)
  ),
    "effective_count_tol"
  )
  expect_error(
    fit_mpcurve(
    X,
    intrinsic_dim = 2L,
    control = mpcurve_control(effective_count_tol = 2)
  ),
    "effective_count_tol"
  )
  expect_error(
    fit_mpcurve(
    X,
    intrinsic_dim = 2L,
    control = mpcurve_control(effective_count_tol = NA_real_)
  ),
    "effective_count_tol"
  )
})

test_that("adaptive EB weights agree with total soft assignment counts", {
  update <- getFromNamespace(".structural_partition_assignment_info", "MPCurver")
  weights <- rbind(c(0.9, 0.1), c(0.25, 0.75), c(0.8, 0.2))
  info <- update(weights, 1, list(partition_prior = "adaptive"))
  expect_equal(nrow(weights) * info$omega, colSums(weights))
  expect_identical(mpcurve_control()$effective_count_tol, 1e-8)
  expect_error(mpcurve_control(effective_weight_tol = 1e-12), "unused argument")
})
