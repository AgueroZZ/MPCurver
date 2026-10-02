named_interface_data <- function() {
  set.seed(614)
  X <- matrix(rnorm(160), 40, 4)
  dimnames(X) <- list(paste0("sample_", seq_len(nrow(X))),
                     paste0("gene_", seq_len(ncol(X))))
  X
}

test_that("result extraction retains dimensions, identities, and posterior values", {
  X <- named_interface_data()
  for (M in c(1L, 2L)) {
    fit <- fit_mpcurve(X, intrinsic_dim = M, num_bins = 5, max_iter = 2)
    for (result in list(fit, do_mpcurve(fit, max_iter = 1))) {
      pos <- fitted_positions(result)
      prob <- fitted_positions(result, "probability")
      sd <- fitted_positions(result, "sd")
      mu <- fitted_trajectories(result)
      traj_sd <- fitted_trajectories(result, "sd")
      cov <- fitted_trajectories(result, "covariance")
      weights <- fitted_assignments(result)
      expect_identical(dim(pos), c(40L, M))
      expect_identical(dim(prob), c(40L, 5L, M))
      expect_identical(dim(mu), c(4L, 5L, M))
      expect_identical(dim(traj_sd), dim(mu))
      expect_identical(dim(cov), c(5L, 5L, 4L, M))
      expect_identical(dim(weights), c(4L, M))
      expect_identical(rownames(pos), rownames(X))
      expect_identical(rownames(weights), colnames(X))
      expect_identical(dimnames(mu)[[1]], colnames(X))
      expect_identical(colnames(pos), colnames(weights))
      expect_identical(colnames(pos), dimnames(mu)[[3]])
      expect_equal(unname(rowSums(weights)), rep(1, 4))
      grid <- seq(0, 1, length.out = 5)
      for (m in seq_len(M)) {
        expect_equal(unname(rowSums(prob[, , m])), rep(1, 40))
        expect_equal(unname(pos[, m]), as.numeric(prob[, , m] %*% grid))
        expect_equal(unname(sd[, m]^2),
                     as.numeric(prob[, , m] %*% grid^2) - unname(pos[, m]^2),
                     tolerance = 1e-12)
        expect_equal(unname(fitted_positions(result, "map")[, m]),
                     grid[max.col(prob[, , m], ties.method = "first")])
        means <- if (M == 1L) result$params$mu else result$params$mu[[m]]
        loc <- if (M == 1L) result$locations else result$locations[[m]]
        expect_equal(unname(mu[, , m]), unname(means))
        expect_identical(rownames(means), colnames(X))
        expect_identical(names(loc$mean$pseudotime), rownames(X))
        for (j in seq_len(4)) {
          expect_equal(unname(traj_sd[j, , m]^2), unname(diag(cov[, , j, m])))
          raw_cov <- if (M == 1L) result$fit$posterior$cov[[j]] else
            result$fit$conditional_posterior$cov[[m]][[j]]
          expect_equal(unname(cov[, , j, m]), unname(raw_cov))
        }
      }
    }
  }
})

test_that("automatic fits expose the same result contract after pruning", {
  set.seed(2)
  X <- matrix(rnorm(30 * 8), 30, 8)
  dimnames(X) <- list(paste0("s", 1:30), paste0("g", 1:8))
  fit <- fit_mpcurve(X, intrinsic_dim = "auto", num_bins = 5, max_iter = 200,
    init_control = mpcurve_init_control(max_intrinsic_dim = 3))
  expect_identical(fit$intrinsic_dim, 1L)
  expect_gt(fit$dimension_estimation$initial_intrinsic_dim, 1)
  for (result in list(fit, do_mpcurve(fit, max_iter = 1))) {
    expect_identical(dim(fitted_positions(result)), c(30L, 1L))
    expect_identical(dim(fitted_assignments(result)), c(8L, 1L))
    expect_identical(dim(fitted_trajectories(result)), c(8L, 5L, 1L))
    expect_identical(rownames(fitted_positions(result)), rownames(X))
    expect_identical(colnames(fitted_assignments(result)),
                     names(result$dimension_estimation$expected_feature_counts))
    expect_equal(fitted_prior(result, ordering = colnames(fitted_positions(result))),
                 fitted_prior(result))
  }
})

test_that("invalid data are rejected before initialization or candidate fitting", {
  X <- named_interface_data()
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    Y <- X
    Y[7, 3] <- bad
    expect_error(fit_mpcurve(Y, initial_method = "random"),
                 "sample 7 .*sample_7.*feature 3 .*gene_3")
    expect_error(select_mpcurve_dimension(Y, max_intrinsic_dim = 2),
                 "X must contain only finite values")
  }
  expect_error(fit_mpcurve(matrix("text", 10, 4)), "X must be numeric")
  expect_error(fit_mpcurve(X + 1i), "real-valued")
  Y <- as.data.frame(X)
  Y$gene_3 <- factor(rep("a", 40))
  expect_error(fit_mpcurve(Y), "feature 3 .*gene_3.*not numeric")
  expect_error(fit_mpcurve(X, S = rep("1", 4)), "S must contain only finite")
  expect_error(fit_mpcurve(X, S = rep(1i, 4)), "S must contain only finite")
  expect_error(fit_mpcurve(X, S = rep(1, 3)), "length-d feature vector")
  expect_error(fit_mpcurve(X, S = matrix(1, 4, 40)), "matching X")
  expect_error(fit_mpcurve(X, S = c(1, 1, -1, 1)), "nonnegative")
  fit <- fit_mpcurve(as.data.frame(X), num_bins = 5, max_iter = 1)
  expect_equal(fit$data, X)
})

test_that("named plotting and one-feature defaults work across plot types", {
  X <- named_interface_data()
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  for (M in c(1L, 2L)) {
    fit <- fit_mpcurve(X, intrinsic_dim = M, num_bins = 5, max_iter = 2)
    for (type in c("scatterplot", "mu")) {
      expect_no_warning(plot(fit, plot_type = type, dims = "gene_3"))
      expect_no_warning(plot(fit, plot_type = type, dims = c("gene_4", "gene_1")))
      expect_error(plot(fit, plot_type = type, dims = "absent"), "Unknown feature name")
      expect_error(plot(fit, plot_type = type, dims = NA), "without missing values")
      expect_error(plot(fit, plot_type = type, dims = 1.5), "integer feature indices")
    }
    expect_identical(.mpcurve_feature_label(fit, 3), "gene_3")
  }
  fit <- fit_mpcurve(X[, 1, drop = FALSE], num_bins = 5, max_iter = 2)
  expect_identical(dim(fitted_trajectories(fit)), c(1L, 5L, 1L))
  expect_no_warning(plot(fit))
  expect_no_warning(plot(fit, plot_type = "mu"))
  colnames(fit$data) <- rep("duplicate", ncol(fit$data))
  expect_identical(.mpcurve_plot_dimensions(NULL, 1L), 1L)
  expect_error(.mpcurve_plot_dimensions("duplicate", 2L, rep("duplicate", 2)), "Ambiguous")
})
