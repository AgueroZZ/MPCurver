automatic_dimension_fixture <- function(S = NULL) {
  set.seed(1501)
  X <- matrix(rnorm(144), nrow = 24, ncol = 6)
  raw <- fit_mpcurve(X, intrinsic_dim = 3, num_bins = 5, max_iter = 2,
    tol = 0, S = S,
    init_control = mpcurve_init_control(),
    control = mpcurve_control(lambda_init = 2, fix_lambda = TRUE))$fit
  raw$dimension_initialization <- list(automatic = TRUE, selected_M = 3L,
    requested_max_intrinsic_dim = 3L, min_cluster_size = 2L)
  raw
}

test_that("automatic finalization removes numerical residue without refitting", {
  finalize <- getFromNamespace(".mpcurve_finalize_automatic_dimension", "MPCurver")
  for (S in list(NULL, rep(0.8, 6),
                 matrix(seq(0.6, 1, length.out = 144), 24, 6))) {
    raw <- automatic_dimension_fixture(S)
    # Leave a small positive mass in A; B and C retain appreciable support.
    raw$pi_weights[,] <- rep(c(1e-10, 0.4, 0.6 - 1e-10), each = raw$d)
    original <- raw
    reduced <- finalize(raw)
    fit <- getFromNamespace("as_mpcurve", "MPCurver")(reduced)
    expect_identical(raw, original)
    expect_identical(fit$intrinsic_dim, 2L)
    expect_identical(fit$model_intrinsic_dim, 2L)
    expect_identical(reduced$ordering_labels, c("B", "C"))
    expect_identical(reduced$gamma, raw$gamma[2:3])
    expect_identical(reduced$conditional_posterior$mean,
                     raw$conditional_posterior$mean[2:3])
    expect_identical(reduced$conditional_posterior$cov,
                     raw$conditional_posterior$cov[2:3])
    expect_identical(reduced$lambda_mat, raw$lambda_mat[, 2:3])
    expect_identical(reduced$params$sigma2, raw$params$sigma2)
    expect_identical(reduced$measurement_sd, S)
    expect_length(fit$fits, 2)
    expect_length(fit$locations, 2)
    expect_equal(rowSums(reduced$pi_weights), rep(1, raw$d))
    expect_equal(as.numeric(fitted_prior(fit, type = "partition")$omega),
                 unname(colMeans(reduced$pi_weights)))
    expect_true(all(reduced$assign %in% c("B", "C")))
    expect_identical(fit$iter, raw$iter)
    expect_length(fit$objective_history, 1)
    expect_length(reduced$weight_history, 1)
    expect_equal(dim(reduced$lambda_trace[[1]]), c(raw$d, 2L))
    # Independently reconstruct the reduced T=1 objective.
    cell <- sum(vapply(1:2, function(m) {
      R <- reduced$gamma[[m]]
      pi <- reduced$position_pi[[m]]
      sum(R * (rep(log(pi), each = nrow(R)) - log(R)))
    }, numeric(1)))
    w <- reduced$pi_weights
    omega <- colMeans(w)
    expected <- cell + sum(w * reduced$local_blocks) +
      sum(w * rep(log(omega), each = nrow(w))) - sum(w * log(w))
    expect_equal(tail(fit$objective_history, 1), expected, tolerance = 1e-10)
    record <- fit$dimension_estimation
    expect_identical(record$initial_intrinsic_dim, 3L)
    expect_identical(record$intrinsic_dim, 2L)
    expect_false(record$pruning[[1]]$refitted)
    expect_identical(record$pruning[[1]]$removed_orderings, "A")
    expect_identical(record$pruning[[1]]$fitting_history$objective,
                     raw$objective_history)
    expect_identical(summary(fit)$dimension_estimation, record)
    continued <- expect_no_warning(do_mpcurve(fit, max_iter = 1, tol = 0))
    expect_identical(continued$intrinsic_dim, 2L)
    expect_equal(continued$iter, fit$iter + 1)
    expect_equal(head(continued$objective_history, 1), fit$objective_history)
    expect_identical(continued$dimension_estimation$pruning, record$pruning)
  }
})

test_that("automatic finalization preserves the scaled prior-weight cutoff", {
  finalize <- getFromNamespace(".mpcurve_finalize_automatic_dimension", "MPCurver")
  raw <- automatic_dimension_fixture()
  raw$pi_weights[,] <- rep(c(1e-10, 0.4, 0.6 - 1e-10), each = raw$d)
  raw$control$effective_count_tol <- colSums(raw$pi_weights)[1]
  expect_identical(finalize(raw)$M, 2L)
  raw$control$effective_count_tol <- 0
  expect_identical(finalize(raw)$M, 3L)
  raw$pi_weights[, 1] <- 0
  raw$pi_weights <- raw$pi_weights / rowSums(raw$pi_weights)
  expect_identical(finalize(raw)$M, 2L)
  raw$control$partition_prior <- "fixed"
  expect_identical(finalize(raw)$M, 3L)
})

test_that("automatic reduction to one ordering returns the ordinary schema", {
  finalize <- getFromNamespace(".mpcurve_finalize_automatic_dimension", "MPCurver")
  for (S in list(NULL, rep(0.8, 6),
                 matrix(seq(0.6, 1, length.out = 144), 24, 6))) {
    raw <- automatic_dimension_fixture(S)
    raw$pi_weights[,] <- rep(c(1e-10, 1 - 2e-10, 1e-10), each = raw$d)
    reduced <- finalize(raw)
    fit <- getFromNamespace("as_mpcurve", "MPCurver")(reduced)
    expect_s3_class(reduced, "cavi")
    expect_identical(fit$intrinsic_dim, 1L)
    expect_null(fit$partition)
    expect_identical(fit$gamma, raw$gamma[[2]])
    expect_identical(fit$params$mu, raw$conditional_posterior$mean[[2]])
    expect_identical(reduced$posterior$cov, raw$conditional_posterior$cov[[2]])
    expect_identical(reduced$lambda_vec, raw$lambda_mat[, 2])
    expect_identical(fit$measurement_sd, S)
    expect_identical(fit$iter, raw$iter)
    expect_no_error(getFromNamespace(".mpcurve_dimension_single_score_check",
                                    "MPCurver")(fit, tail(fit$elbo_trace, 1)))
    continued <- expect_no_warning(do_mpcurve(fit, max_iter = 1, tol = 0))
    expect_identical(continued$intrinsic_dim, 1L)
    expect_equal(continued$iter, fit$iter + 1)
    expect_identical(continued$dimension_estimation, fit$dimension_estimation)
    pdf_file <- tempfile(fileext = ".pdf")
    grDevices::pdf(pdf_file)
    expect_no_error(plot(fit, plot_type = "mu", dims = 1))
    expect_no_error(plot(fit, plot_type = "scatterplot", dims = c(1, 2)))
    grDevices::dev.off()
    unlink(pdf_file)
  }
})

test_that("public auto fitting finalizes dimension while an explicit M is kept", {
  sim <- simulate_intrinsic_trajectories(n = 40, d_signal = c(4, 2),
    d_noise = 0, noise_sd = 0.01, seed = 404,
    trajectory_family = c("monotone", "monotone"))
  args <- list(X = sim$X, num_bins = 5, max_iter = 0,
    init_control = mpcurve_init_control(max_intrinsic_dim = 3),
    # A larger cutoff deliberately exercises public pruning at initialization.
    control = mpcurve_control(effective_count_tol = 2.5))
  expect_warning(
    automatic <- do.call(fit_mpcurve, c(args, list(intrinsic_dim = "auto"))),
    "effective_count_tol exceeds 1e-6", fixed = TRUE)
  explicit <- do.call(fit_mpcurve, c(args, list(intrinsic_dim = 2)))
  expect_identical(automatic$dimension_initialization$selected_M, 2L)
  expect_identical(automatic$intrinsic_dim, 1L)
  expect_true(is.matrix(automatic$gamma))
  expect_false(automatic$dimension_estimation$pruning[[1]]$refitted)
  expect_identical(explicit$intrinsic_dim, 2L)
  expect_length(explicit$conditional_posterior$mean, 2)
  expect_null(explicit$dimension_estimation)
  fixed <- do.call(fit_mpcurve, c(args, list(intrinsic_dim = "auto",
                                           partition_prior = "fixed")))
  expect_identical(fixed$intrinsic_dim, 2L)
  expect_length(fixed$dimension_estimation$pruning, 0)
})

test_that("automatic single-ordering initialization keeps its estimation record", {
  set.seed(1502)
  X <- matrix(rnorm(144), 24, 6)
  fit <- fit_mpcurve(X, intrinsic_dim = "auto", num_bins = 5, max_iter = 1,
    tol = 0, init_control = mpcurve_init_control(max_intrinsic_dim = 1))
  expect_identical(fit$intrinsic_dim, 1L)
  expect_identical(fit$dimension_estimation$initial_intrinsic_dim, 1L)
  expect_length(fit$dimension_estimation$pruning, 0)
  continued <- expect_no_warning(do_mpcurve(fit, max_iter = 1, tol = 0))
  expect_identical(continued$dimension_estimation, fit$dimension_estimation)
  expect_equal(continued$iter, fit$iter + 1)
})

test_that("later removal preserves earlier dimension-estimation history", {
  finalize <- getFromNamespace(".mpcurve_finalize_automatic_dimension", "MPCurver")
  raw <- automatic_dimension_fixture()
  raw$pi_weights[,] <- rep(c(1e-10, 0.4, 0.6 - 1e-10), each = raw$d)
  first <- finalize(raw)
  first$pi_weights[,] <- rep(c(1e-10, 1 - 1e-10), each = raw$d)
  second <- finalize(first)
  expect_s3_class(second, "cavi")
  expect_identical(second$iter, raw$iter)
  expect_length(second$dimension_estimation$pruning, 2)
  expect_identical(second$dimension_estimation$pruning[[1]],
                   first$dimension_estimation$pruning[[1]])
  expect_identical(second$dimension_estimation$pruning[[2]]$removed_orderings, "B")
  expect_identical(second$dimension_estimation$initial_intrinsic_dim, 3L)
  expect_identical(second$dimension_estimation$intrinsic_dim, 1L)
  expect_identical(second$gamma, raw$gamma[[3]])
})

test_that("large numerical-zero tolerances warn only when adaptive pruning applies", {
  finalize <- getFromNamespace(".mpcurve_finalize_automatic_dimension", "MPCurver")
  raw <- automatic_dimension_fixture()
  for (tol in c(0, 1e-8, 1e-6)) {
    raw$control$effective_count_tol <- tol
    expect_no_warning(finalize(raw))
  }
  raw$control$effective_count_tol <- 1.001e-6
  expect_warning(finalize(raw), "effective_count_tol exceeds 1e-6", fixed = TRUE)
  raw$control$partition_prior <- "fixed"
  expect_no_warning(finalize(raw))
  raw$control$partition_prior <- "adaptive"
  raw$dimension_initialization$automatic <- FALSE
  expect_no_warning(finalize(raw))

  X <- simulate_intrinsic_trajectories(n = 40, d_signal = c(4, 2),
    d_noise = 0, noise_sd = 0.01, seed = 404,
    trajectory_family = c("monotone", "monotone"))$X
  for (control in list(list(effective_count_tol = 1e-5),
                      mpcurve_control(effective_count_tol = 1e-5))) {
    warnings <- character()
    fit <- withCallingHandlers(
      fit_mpcurve(X, intrinsic_dim = "auto", num_bins = 5, max_iter = 1,
        init_control = mpcurve_init_control(max_intrinsic_dim = 2),
        control = control),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    expect_length(warnings, 1)
    expect_match(warnings, "estimated prior weight <= effective_count_tol / ncol(X)",
                 fixed = TRUE)
    expect_true(is.finite(tail(fit$elbo_trace, 1)))
    expect_warning(do_mpcurve(fit, max_iter = 1),
                   "effective_count_tol exceeds 1e-6", fixed = TRUE)
    expect_no_warning(fit_mpcurve(X, intrinsic_dim = 2, num_bins = 5,
                                  max_iter = 1, control = control))
    expect_no_warning(fit_mpcurve(X, intrinsic_dim = "auto", num_bins = 5,
      max_iter = 1, partition_prior = "fixed", control = control,
      init_control = mpcurve_init_control(max_intrinsic_dim = 2)))
  }
})

test_that("pruning rejects zero retained mass for soft assignments", {
  sim <- simulate_dual_trajectory(
    n = 40, d1 = 6, d2 = 4, d_noise = 0,
    noise_sd = 0.01,
    trajectory_family = c("quadratic", "monotone"), seed = 1)
  expect_warning(
    expect_error(
      fit_mpcurve(sim$X, intrinsic_dim = "auto",
        num_bins = 10, max_iter = 20,
        tol = 0, init_control = mpcurve_init_control(max_intrinsic_dim = 2),
        control = mpcurve_control(effective_count_tol = 4.5)),
      "no positive assignment probability", fixed = TRUE),
    "effective_count_tol exceeds 1e-6", fixed = TRUE)
})
