test_that("PCA ordering supports explicit components", {
  set.seed(1)
  X <- matrix(rnorm(100), 20, 5)
  ord1 <- PCA_ordering(X, component = 1)
  ord2 <- PCA_ordering(X, component = 2)
  expect_length(ord1$t, nrow(X))
  expect_length(ord2$t, nrow(X))
  expect_false(isTRUE(all.equal(ord1$t, ord2$t)))
})

test_that("soft_partition_cavi returns canonical structural state for M=2", {
  set.seed(2)
  X <- matrix(rnorm(120), 24, 5)
  res <- suppressWarnings(soft_partition_cavi(
    X,
    M = 2,
    init_methods = c("PCA", "PCA"),
    pca_components = 1:2,
    partition_init = "ordering_methods",
    K = 5,
    T_start = 1,
    T_end = 1,
    n_outer = 1,
    max_converge_iter = 2,
    ridge = 0,
    verbose = FALSE
  ))
  expect_s3_class(res, "soft_partition_cavi")
  expect_identical(res$variational_family, "structured")
  expect_equal(dim(res$pi_weights), c(ncol(X), 2L))
  expect_equal(rowSums(res$pi_weights), rep(1, ncol(X)), tolerance = 1e-10)
  expect_equal(length(res$params$sigma2), ncol(X))
  expect_equal(dim(res$lambda_mat), c(ncol(X), 2L))
  expect_length(res$gamma, 2L)
  expect_length(res$conditional_posterior$mean, 2L)
  expect_null(res$active_orderings)
  expect_null(res$active_feature_pairs)
  expect_null(res$effective_pi_weights)
  expect_true(all(diff(res$objective_history) >= -1e-6))
})

test_that("fit_mpcurve supports fixed structural M=2 and M=3", {
  set.seed(3)
  X <- matrix(rnorm(168), 28, 6)
  fit2 <- suppressWarnings(fit_mpcurve(
    X,
    intrinsic_dim = 2,
    method = c("PCA", "PCA"),
    partition_init = "ordering_methods",
    K = 5,
    T_start = 1,
    T_end = 1,
    n_outer = 1,
    max_converge_iter = 1,
    verbose = FALSE
  ))
  fit3 <- suppressWarnings(fit_mpcurve(
    X,
    intrinsic_dim = 3,
    method = rep("PCA", 3),
    partition_init = "ordering_methods",
    K = 5,
    T_start = 1,
    T_end = 1,
    n_outer = 1,
    max_converge_iter = 1,
    verbose = FALSE
  ))
  expect_s3_class(fit2, "mpcurve")
  expect_s3_class(fit3, "mpcurve")
  expect_equal(fit2$intrinsic_dim, 2L)
  expect_equal(fit3$intrinsic_dim, 3L)
  expect_equal(fit2$params$sigma2, fit2$fit$params$sigma2)
  expect_equal(fit3$params$sigma2, fit3$fit$params$sigma2)
  expect_equal(length(fit2$fits), 2L)
  expect_equal(length(fit3$fits), 3L)
})

test_that("partition lambda and sigma initialization are shared and explicit", {
  set.seed(4)
  X <- matrix(rnorm(120), 24, 5)
  sigma_init <- seq(0.5, 1.5, length.out = ncol(X))
  fit <- suppressWarnings(fit_mpcurve(
    X,
    intrinsic_dim = 2,
    method = c("PCA", "PCA"),
    partition_init = "ordering_methods",
    K = 5,
    lambda = 2,
    fix_lambda = TRUE,
    sigma2_init = sigma_init,
    T_start = 1,
    T_end = 1,
    n_outer = 1,
    max_converge_iter = 0,
    verbose = FALSE
  ))
  expect_equal(unname(fit$lambda_mat), matrix(2, ncol(X), 2), tolerance = 0)
  expect_equal(fit$sigma2_trace[[1]], sigma_init, tolerance = 0)
  expect_true(all(vapply(fit$fits, function(child) {
    isTRUE(all.equal(child$params$sigma2, fit$params$sigma2, tolerance = 0))
  }, logical(1))))
})

test_that("do_mpcurve continues only from canonical structural state", {
  set.seed(5)
  X <- matrix(rnorm(100), 20, 5)
  fit <- suppressWarnings(fit_mpcurve(
    X,
    intrinsic_dim = 2,
    method = c("PCA", "PCA"),
    partition_init = "ordering_methods",
    K = 5,
    T_start = 1,
    T_end = 1,
    n_outer = 1,
    max_converge_iter = 1,
    verbose = FALSE
  ))
  before <- length(fit$objective_history)
  continued <- suppressWarnings(do_mpcurve(fit, iter = 2))
  expect_gt(length(continued$objective_history), before)
  expect_identical(continued$variational_family, "structured")

  legacy <- fit
  legacy$fit$variational_family <- NULL
  legacy$fit$control$variational_family <- NULL
  expect_error(do_mpcurve(legacy, iter = 1), "read-only")
})

test_that("partition initialization retains a common K", {
  set.seed(6)
  X <- cbind(
    rep(c(-1, 1), each = 20),
    rep(c(-2, 2), each = 20),
    rnorm(40),
    rnorm(40)
  )
  fit <- suppressWarnings(fit_mpcurve(
    X,
    intrinsic_dim = 2,
    method = c("PCA", "PCA"),
    partition_init = "ordering_methods",
    K = 6,
    discretization = "quantile",
    n_outer = 1,
    max_converge_iter = 0,
    verbose = FALSE
  ))
  expect_equal(fit$K, 6L)
  expect_true(all(vapply(fit$gamma, ncol, integer(1)) == 6L))
})

test_that("soft two-trajectory wrapper delegates to structural VI", {
  set.seed(7)
  X <- matrix(rnorm(100), 20, 5)
  fit <- suppressWarnings(soft_two_trajectory_cavi(
    X,
    K = 5,
    n_outer = 1,
    max_converge_iter = 1,
    verbose = FALSE
  ))
  expect_s3_class(fit, "soft_partition_cavi")
  expect_identical(fit$variational_family, "structured")
  expect_equal(fit$M, 2L)
})

test_that("plot.mpcurve uses canonical structural means for trajectory plots", {
  set.seed(8)
  X <- matrix(rnorm(100), 20, 5)
  fit <- suppressWarnings(fit_mpcurve(
    X,
    intrinsic_dim = 2,
    method = c("PCA", "PCA"),
    partition_init = "ordering_methods",
    K = 5,
    n_outer = 1,
    max_converge_iter = 0,
    verbose = FALSE
  ))

  expect_identical(
    names(formals(plot.mpcurve)),
    c("x", "plot_type", "dims", "data", "pal", "add_legend", "...")
  )

  fit_before <- serialize(fit, NULL)

  # A structural trajectory plot must not read the derived compatibility fits.
  fit_for_mu <- fit
  fit_for_mu$fits <- lapply(fit_for_mu$fits, function(child) {
    child$params$mu <- NULL
    child$posterior$mean <- NULL
    child
  })
  mu_fit_before <- serialize(fit_for_mu, NULL)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(fit))
  expect_invisible(plot(fit_for_mu, plot_type = "mu", dims = 1L))
  expect_invisible(plot(fit_for_mu, plot_type = "mu", dims = c(1L, 2L)))
  expect_invisible(plot(fit_for_mu, plot_type = "elbo"))
  expect_identical(serialize(fit, NULL), fit_before)
  expect_identical(serialize(fit_for_mu, NULL), mu_fit_before)

  missing_canonical_means <- fit
  missing_canonical_means$conditional_posterior$mean <- NULL
  expect_error(
    plot(missing_canonical_means, plot_type = "mu", dims = 1L),
    "conditional_posterior"
  )

  legacy_view <- fit
  legacy_view$variational_family <- NULL
  legacy_view$control$variational_family <- NULL
  legacy_view$conditional_posterior <- NULL
  expect_invisible(plot(legacy_view, plot_type = "mu", dims = 1L))
})

test_that("single-ordering trajectory plots honor the public dims argument", {
  set.seed(81)
  X <- matrix(rnorm(80), 20, 4)
  fit <- fit_mpcurve(X, intrinsic_dim = 1, K = 5, iter = 2, verbose = FALSE)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(fit, plot_type = "mu", dims = 4L))
  expect_error(
    plot(fit, plot_type = "mu", dims = 5L),
    "dims out of range"
  )
})

test_that("intrinsic_dim=1 remains the ordinary cavi path", {
  set.seed(9)
  X <- matrix(rnorm(80), 20, 4)
  fit <- fit_mpcurve(X, intrinsic_dim = 1, K = 5, iter = 2, verbose = FALSE)
  expect_s3_class(fit$fit, "cavi")
  expect_null(fit$partition)
  expect_equal(fit$intrinsic_dim, 1L)
})

test_that("non-cavi public partition requests still fail", {
  X <- matrix(rnorm(60), 15, 4)
  expect_error(
    fit_mpcurve(X, algorithm = "csmooth_em", intrinsic_dim = 2),
    "CAVI-only"
  )
})
