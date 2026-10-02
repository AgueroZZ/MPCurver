mpcurve_active_dim <- function(fit) {
  if (!inherits(fit$fit, "soft_partition_cavi")) {
    return(1L)
  }
  active <- fit$fit$active_orderings
  if (is.null(active)) {
    M <- fit$fit$M
    if (is.null(M)) {
      M <- length(fit$fit$fits)
    }
    return(as.integer(M))
  }
  as.integer(sum(active))
}

test_that("fit_mpcurve defaults to cavi and returns a valid mpcurve object", {
  sim <- simulate_mpcurve(
    n = 100,
    d = 10,
    num_bins = 6,
    seed = 11,
    control = list(rw_order = 2)
  )

  fit <- fit_mpcurve(
    sim$X,
    num_bins = 6,
    max_iter = 40
  )

  expect_s3_class(fit, "mpcurve")
  expect_equal(fit$algorithm, "cavi")
  expect_equal(dim(fit$params$mu), c(10, 6))
  expect_equal(dim(fit$gamma), c(100, 6))
  expect_true(!is.null(fit$locations))
  expect_length(fit$locations$mean$index, 100)
  expect_length(fit$locations$map$index, 100)
  expect_true(all(fit$locations$mean$pseudotime >= 0 & fit$locations$mean$pseudotime <= 1))
  expect_true(all(fit$locations$map$pseudotime >= 0 & fit$locations$map$pseudotime <= 1))
  expect_equal(length(fit$elbo_trace), fit$iter + 1L)
  expect_gte(min(diff(fit$elbo_trace)), -1e-8)

  s <- summary(fit)
  expect_s3_class(s, "summary.mpcurve")
  expect_s3_class(s$underlying, "summary.cavi")
})

test_that("single-ordering fits expose top-level data and fitted_prior access", {
  sim <- simulate_mpcurve(
    n = 60,
    d = 6,
    num_bins = 5,
    seed = 14,
    control = list(rw_order = 2)
  )

  raw_fit <- cavi(
    sim$X,
    K = 5,
    method = "PCA",
    max_iter = 3,
    verbose = FALSE
  )
  fit <- as_mpcurve(raw_fit)

  expect_true(is.matrix(raw_fit$data))
  expect_true(is.list(raw_fit$priors))
  expect_true(is.matrix(fit$data))
  expect_true(is.list(fit$priors))
  expect_equal(fitted_prior(raw_fit, type = "position"), raw_fit$params$pi)
  expect_equal(fitted_prior(fit, type = "position"), fit$params$pi)
})

test_that("do_cavi and do_mpcurve extend cavi traces", {
  sim <- simulate_mpcurve(
    n = 90,
    d = 8,
    num_bins = 5,
    seed = 12,
    control = list(rw_order = 2)
  )

  raw_fit <- cavi(
    sim$X,
    K = 5,
    method = "PCA",
    rw_q = 2,
    max_iter = 4,
    tol = 0,
    verbose = FALSE
  )
  raw_more <- do_cavi(raw_fit, iter = 3, tol = 0, verbose = FALSE)

  expect_s3_class(raw_more, "cavi")
  expect_gt(length(raw_more$elbo_trace), length(raw_fit$elbo_trace))
  expect_equal(length(raw_more$elbo_trace), raw_more$iter + 1L)
  expect_gte(min(diff(raw_more$elbo_trace)), -1e-8)

  mp_fit <- as_mpcurve(raw_fit)
  mp_more <- do_mpcurve(mp_fit, max_iter = 3, verbose = FALSE)

  expect_s3_class(mp_more, "mpcurve")
  expect_equal(mp_more$algorithm, "cavi")
  expect_gt(length(mp_more$elbo_trace), length(mp_fit$elbo_trace))
  expect_equal(length(mp_more$elbo_trace), mp_more$iter + 1L)
  expect_equal(length(mp_more$locations$mean$index), 90)
  expect_gte(min(diff(mp_more$elbo_trace)), -1e-8)
})

test_that("fit_mpcurve and do_mpcurve expose a cavi-only public wrapper", {
  fit_formals <- names(formals(fit_mpcurve))
  do_formals <- names(formals(do_mpcurve))

  expect_false("greedy" %in% fit_formals)
  expect_true("S" %in% fit_formals)
  expect_false("algorithm" %in% fit_formals)
  expect_true("position_prior" %in% fit_formals)
  expect_true("control" %in% fit_formals)
  expect_true("partition_prior" %in% fit_formals)
  expect_true("init_control" %in% fit_formals)
  expect_false("relative_lambda" %in% fit_formals)
  expect_false("adaptive" %in% fit_formals)
  expect_false("sigma_update" %in% fit_formals)
  expect_false("check_decrease" %in% fit_formals)
  expect_false("tol_decrease" %in% fit_formals)

  expect_false("adaptive" %in% do_formals)
  expect_false("S" %in% do_formals)
  expect_false("sigma_update" %in% do_formals)
  expect_false("check_decrease" %in% do_formals)
  expect_false("tol_decrease" %in% do_formals)
})

test_that("fit_mpcurve forwards fixed position_prior to the cavi backend", {
  sim <- simulate_mpcurve(
    n = 70,
    d = 7,
    num_bins = 5,
    seed = 45,
    control = list(rw_order = 2)
  )
  pi_fixed <- c(0.4, 0.25, 0.15, 0.1, 0.1)

  fit <- fit_mpcurve(
    sim$X,
    num_bins = 5,
    position_prior = "fixed",
    max_iter = 3,
    tol = 0,
    verbose = FALSE,
    control = mpcurve_control(position_prior_weights = pi_fixed)
  )

  expect_equal(fit$params$pi, pi_fixed, tolerance = 1e-10)
  expect_equal(fit$fit$control$position_prior, "fixed")
})

test_that("fit_mpcurve and do_mpcurve preserve known measurement sd", {
  sim <- simulate_mpcurve(
    n = 70,
    d = 6,
    num_bins = 5,
    seed = 44,
    control = list(rw_order = 2)
  )
  S <- seq(0.09, 0.14, length.out = ncol(sim$X))

  fit <- fit_mpcurve(
    sim$X,
    num_bins = 5,
    S = S,
    max_iter = 4,
    tol = 0,
    verbose = FALSE
  )
  fit_more <- do_mpcurve(fit, max_iter = 2, tol = 0, verbose = FALSE)

  expect_s3_class(fit, "mpcurve")
  expect_equal(fit$measurement_sd, S, tolerance = 1e-12)
  expect_equal(fit$fit$control$noise_model, "known_feature_sd")
  expect_null(fit$params$sigma2)
  expect_equal(fit_more$measurement_sd, S, tolerance = 1e-12)
  expect_equal(fit_more$fit$control$noise_model, "known_feature_sd")
  expect_null(fit_more$params$sigma2)
  expect_gte(min(diff(fit_more$elbo_trace)), -1e-8)
})









test_that("legacy hard CAVI partition selectors are no longer exported", {
  ns_exports <- getNamespaceExports("MPCurver")
  expect_false("forward_two_ordering_partition_cavi" %in% ns_exports)
  expect_false("backward_two_ordering_partition_cavi" %in% ns_exports)
})

test_that("legacy wrapper controls are rejected cleanly", {
  sim <- simulate_mpcurve(
    n = 60,
    d = 6,
    num_bins = 4,
    seed = 13,
    control = list(rw_order = 2)
  )

  expect_s3_class(
    fit_mpcurve(
    sim$X,
    initial_method = "PCA",
    num_bins = 4,
    max_iter = 2,
    verbose = FALSE
  ),
    "mpcurve"
  )

  expect_error(
    fit_mpcurve(
    sim$X,
    algorithm = "csmooth_em",
    initial_method = "PCA",
    num_bins = 4,
    max_iter = 2
  ),
    "unused argument"
  )

  expect_error(
    fit_mpcurve(
    sim$X,
    num_bins = 4,
    max_iter = 2,
    adaptive = "ml"
  ),
    "unused argument"
  )

  legacy_raw <- initialize_csmoothEM(sim$X, method = "PCA", K = 4, num_iter = 1)
  legacy_mp <- as_mpcurve(legacy_raw)
  expect_error(
    do_mpcurve(legacy_mp, max_iter = 1),
    "Only CAVI fits"
  )
})
