similarity_dimension_toy_distance <- function() {
  distance <- matrix(0.95, nrow = 7L, ncol = 7L)
  diag(distance) <- 0

  for (pair in list(c(1L, 2L), c(3L, 4L), c(5L, 6L))) {
    distance[pair, pair] <- 0.05
  }
  diag(distance) <- 0

  distance[7L, c(1L, 2L)] <- 0.80
  distance[c(1L, 2L), 7L] <- 0.80
  distance[c(1L, 2L), c(3L, 4L)] <- 0.90
  distance[c(3L, 4L), c(1L, 2L)] <- 0.90
  distance
}

test_that("similarity dimension selection maximizes silhouette over feasible cuts", {
  select_dimension <- ns_fn(".cavi_select_similarity_dimension")
  distance <- similarity_dimension_toy_distance()

  unrestricted <- select_dimension(
    distance = distance,
    cluster_linkage = "single",
    max_intrinsic_dim = 5L,
    min_cluster_size = 1L
  )
  constrained <- select_dimension(
    distance = distance,
    cluster_linkage = "single",
    max_intrinsic_dim = 5L,
    min_cluster_size = 2L
  )

  expect_identical(unrestricted$selected_M, 4L)
  expect_equal(sort(unname(unrestricted$cluster_sizes)), c(1L, 2L, 2L, 2L))

  expect_identical(constrained$selected_M, 3L)
  expect_equal(sort(unname(constrained$cluster_sizes)), c(2L, 2L, 3L))
  expect_true(all(constrained$cluster_sizes >= 2L))
  expect_equal(length(constrained$feature_cluster), 7L)
  expect_setequal(unique(unname(constrained$feature_cluster)), seq_len(3L))

  diagnostics <- constrained$diagnostics
  expect_s3_class(diagnostics, "data.frame")
  expect_true(all(c(
    "M", "mean_silhouette", "minimum_cluster_size", "feasible"
  ) %in% names(diagnostics)))

  row_M2 <- diagnostics[diagnostics$M == 2L, , drop = FALSE]
  row_M3 <- diagnostics[diagnostics$M == 3L, , drop = FALSE]
  row_M4 <- diagnostics[diagnostics$M == 4L, , drop = FALSE]
  expect_equal(row_M2$mean_silhouette, 0.4436090, tolerance = 1e-6)
  expect_equal(row_M3$mean_silhouette, 0.7141566, tolerance = 1e-6)
  expect_equal(row_M4$mean_silhouette, 0.8083751, tolerance = 1e-6)
  expect_true(row_M3$feasible)
  expect_false(row_M4$feasible)
  expect_identical(row_M4$minimum_cluster_size, 1L)
})

test_that("similarity dimension selection validates controls and falls back to one", {
  select_dimension <- ns_fn(".cavi_select_similarity_dimension")
  distance <- similarity_dimension_toy_distance()

  fallback <- select_dimension(
    distance = distance,
    cluster_linkage = "single",
    max_intrinsic_dim = 5L,
    min_cluster_size = 4L
  )
  expect_identical(fallback$selected_M, 1L)
  expect_identical(unname(fallback$cluster_sizes), 7L)
  expect_identical(unname(fallback$feature_cluster), rep(1L, 7L))
  expect_false(any(fallback$diagnostics$feasible))

  for (bad_minimum in list(0, -1, 1.5, NA_real_, Inf, c(2, 3))) {
    expect_error(
      select_dimension(distance, "single", 5L, bad_minimum),
      "min_cluster_size"
    )
  }
  for (bad_maximum in list(0, -1, 2.5, NA_real_, Inf, c(3, 4))) {
    expect_error(
      select_dimension(distance, "single", bad_maximum, 2L),
      "max_intrinsic_dim"
    )
  }
  expect_error(
    select_dimension(distance, "not-a-linkage", 5L, 2L),
    "cluster_linkage"
  )
})

test_that("spline-r2 similarity is symmetric, bounded, and uses fixed-df fits", {
  compute_similarity <- ns_fn(".compute_same_ordering_similarity")
  n <- 40L
  position <- (seq_len(n) - 0.5) / n
  set.seed(901)
  X <- cbind(
    anchor = position,
    smooth = 2 * (position - 0.35)^2 - 0.3 * position,
    unrelated = stats::rnorm(n),
    constant = rep(1, n)
  )

  similarity <- compute_similarity(
    X = X,
    metric = "spline_r2",
    spline_r2_df = 5L,
    min_feature_sd = 1e-8
  )

  S <- similarity$S
  expect_equal(dim(S), c(4L, 4L))
  expect_equal(S, t(S), tolerance = 1e-12)
  expect_true(all(is.finite(S)))
  expect_true(all(S >= 0 & S <= 1))
  expect_equal(unname(diag(S)), rep(1, 4L))
  expect_equal(similarity$distance, 1 - S, tolerance = 1e-12)
  expect_gt(S["anchor", "smooth"], 0.95)
  expect_equal(unname(S["constant", -4L]), rep(0, 3L))
  expect_true(similarity$feature_info$low_variance[4L])

  spline_basis <- cbind(
    "(Intercept)" = 1,
    splines::ns(position, df = 5L, intercept = FALSE)
  )
  directional_r2 <- function(predictor, response) {
    ordered_response <- response[order(predictor, method = "radix")]
    fit <- stats::lm.fit(spline_basis, ordered_response)
    total_sum_squares <- sum((response - mean(response))^2)
    1 - sum(fit$residuals^2) / total_sum_squares
  }
  expected_pair_similarity <- max(
    directional_r2(X[, "anchor"], X[, "smooth"]),
    directional_r2(X[, "smooth"], X[, "anchor"])
  )
  expected_pair_similarity <- min(max(expected_pair_similarity, 0), 1)
  expect_equal(
    unname(S["anchor", "smooth"]),
    expected_pair_similarity,
    tolerance = 1e-10
  )

  expect_error(
    compute_similarity(X, metric = "spline_r2", spline_r2_df = 0),
    "spline_r2_df"
  )
  expect_error(
    compute_similarity(X, metric = "spline_r2", spline_r2_df = 2.5),
    "spline_r2_df"
  )
})

test_that("minimum cluster size does not change an explicit numeric M", {
  set.seed(92)
  X <- cbind(
    anchor = seq_len(12L) / 12,
    noise_1 = stats::rnorm(12L),
    noise_2 = stats::rnorm(12L)
  )

  fit <- suppressWarnings(fit_mpcurve(
    X,
    intrinsic_dim = 3L,
    similarity_min_cluster_size = 2L,
    partition_init = "similarity",
    similarity_metric = "spline_r2",
    spline_r2_df = 3L,
    K = 3L,
    T_start = 1,
    T_end = 1,
    n_outer = 1L,
    inner_iter = 1L,
    max_converge_iter = 0L,
    verbose = FALSE
  ))

  expect_identical(fit$intrinsic_dim, 3L)
  expect_null(fit$dimension_initialization)
  expect_equal(
    sort(unname(fit$similarity_init$cluster_sizes)),
    rep(1L, 3L)
  )
})

test_that("fit_mpcurve auto dimension uses cluster-size adaptive initialization", {
  sim <- simulate_intrinsic_trajectories(
    n = 40,
    d_signal = c(4, 2),
    d_noise = 0,
    sigma = 0.01,
    seed = 404,
    trajectory_family = c("monotone", "monotone")
  )

  fit <- suppressWarnings(fit_mpcurve(
    sim$X,
    intrinsic_dim = "auto",
    max_intrinsic_dim = 3L,
    similarity_min_cluster_size = 2L,
    partition_prior = "adaptive",
    partition_init = "similarity",
    spline_r2_df = 5L,
    cluster_linkage = "single",
    K = 5L,
    T_start = 1,
    T_end = 1,
    n_outer = 1L,
    inner_iter = 1L,
    max_converge_iter = 0L,
    verbose = FALSE
  ))

  expect_s3_class(fit, "mpcurve")
  expect_identical(fit$intrinsic_dim, 2L)
  expect_true(is.list(fit$dimension_initialization))

  dimension_init <- fit$dimension_initialization
  expect_identical(dimension_init$selected_M, 2L)
  expect_equal(sort(unname(dimension_init$cluster_sizes)), c(2L, 4L))
  expect_true(all(dimension_init$cluster_sizes >= 2L))
  expect_equal(
    sort(unname(dimension_init$initial_partition_probabilities)),
    c(1 / 3, 2 / 3),
    tolerance = 1e-12
  )

  initial_weights <- fit$fit$weight_history[[1L]]
  expect_equal(
    sort(unname(colMeans(initial_weights))),
    c(1 / 3, 2 / 3),
    tolerance = 1e-12
  )
  expect_true(all(abs(rowSums(initial_weights) - 1) < 1e-12))
  expect_equal(
    initial_weights,
    matrix(
      rep(dimension_init$initial_partition_probabilities, each = ncol(sim$X)),
      nrow = ncol(sim$X),
      ncol = 2L,
      dimnames = dimnames(initial_weights)
    ),
    tolerance = 1e-12
  )

  expect_identical(fit$fit$control$similarity_metric, "spline_r2")
  expect_identical(fit$fit$control$spline_r2_df, 5L)
  expect_identical(fit$fit$control$similarity_min_cluster_size, 2L)
})
