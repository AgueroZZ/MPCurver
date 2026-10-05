# Independent exhaustive reference using all pairwise Euclidean distances.
isomap_reference_graph <- function(X, k) {
  distances <- as.matrix(stats::dist(X))
  diag(distances) <- Inf
  neighbors <- t(vapply(seq_len(nrow(X)), function(row) {
    order(distances[row, ], seq_len(nrow(X)))[seq_len(k)]
  }, integer(k)))
  edges <- cbind(rep(seq_len(nrow(X)), each = k), as.vector(t(neighbors)))
  igraph::graph_from_edgelist(edges, directed = FALSE)
}

test_that("automatic Isomap finds the exact smallest connected neighborhood", {
  set.seed(2602)
  fixtures <- list(
    matrix(c(0, 1, 3, 7, 15), ncol = 1), # Already connected at k=1.
    matrix(c(0, 1, 2, 20, 21, 22, 40, 41, 42), ncol = 1),
    matrix(rnorm(80), 40, 2),
    matrix(c(0, 1, 10), ncol = 1) # Smallest supported sample count.
  )
  for (X in fixtures) {
    connected <- vapply(seq_len(nrow(X) - 1L), function(k) {
      igraph::is_connected(isomap_reference_graph(X, k))
    }, logical(1))
    expected <- which(connected)[1L]
    automatic <- expect_no_warning(isomap_ordering(X, seed = 19))
    expect_identical(automatic$k_used, expected)
    expect_equal(automatic$n_components, 1L)
    expect_identical(automatic$keep_idx, seq_len(nrow(X)))
    expect_true(all(is.finite(automatic$t)))
    expect_equal(automatic, isomap_ordering(X, num_neighbors = expected, seed = 19))
  }
})

test_that("explicit Isomap neighborhoods stay fixed when disconnected", {
  X <- matrix(c(0, 1, 2, 20, 21, 22, 40, 41, 42), ncol = 1)
  expect_warning(fixed <- isomap_ordering(X, num_neighbors = 2, seed = 19),
                 "connected components")
  expect_identical(fixed$k_used, 2L)
  expect_equal(fixed$n_components, 3L)
  expect_length(fixed$keep_idx, 3L)
  expect_identical(sum(is.na(fixed$t)), 6L)
  expect_gt(isomap_ordering(X, seed = 19)$k_used, 2L)
})

test_that("ties and duplicate observations use nested graphs without self neighbors", {
  fixtures <- list(
    matrix(c(0, 0, 0, 1, 1, 2, 2, 5), ncol = 1),
    cbind(rep(0:3, 4), rep(0:3, each = 4)),
    matrix(rep(0, 8), ncol = 1)
  )
  for (X in fixtures) {
    minimum <- NA_integer_
    for (k in seq_len(nrow(X) - 1L)) {
      neighbors <- MPCurver:::.isomap_neighbors(X, k)
      expect_true(all(neighbors$idx != row(neighbors$idx)))
      expect_true(all(apply(neighbors$idx, 1, function(ids) !anyDuplicated(ids))))
      graph <- MPCurver:::.isomap_graph(X, k)
      reference <- isomap_reference_graph(X, k)
      expect_equal(igraph::as_adjacency_matrix(graph$graph, sparse = FALSE),
                   igraph::as_adjacency_matrix(reference, sparse = FALSE),
                   ignore_attr = TRUE)
      if (k > 1L) expect_identical(neighbors$idx[, seq_len(k - 1L), drop = FALSE], previous)
      previous <- neighbors$idx
      if (is.na(minimum) && igraph::is_connected(reference)) minimum <- k
    }
    expect_identical(MPCurver:::.isomap_connected_graph(X)$k, minimum)
  }
})

test_that("single and grouped fits preserve their own realized Isomap counts", {
  set.seed(720)
  X <- matrix(rnorm(160), 40, 4)
  for (dimension in list(1L, 2L, "auto")) {
    fit <- expect_no_warning(fit_mpcurve(X, num_bins = 5,
      intrinsic_dim = dimension, initial_method = "isomap", max_iter = 0,
      init_control = mpcurve_init_control(max_intrinsic_dim = 2,
        on_failure = "error", method_args = list(num_neighbors = NULL, seed = 19))))
    info <- fit$fit$init_info
    if (fit$intrinsic_dim == 1L) {
      expect_identical(info$k_used, isomap_ordering(X, seed = 19)$k_used)
    } else {
      for (entry in info) {
        expect_identical(entry$k_used,
          isomap_ordering(X[, entry$feature_idx, drop = FALSE], seed = 19)$k_used)
      }
    }
    expect_identical(do_mpcurve(fit, max_iter = 1)$fit$init_info, info)
  }
  fixed <- fit_mpcurve(X, num_bins = 5, initial_method = "isomap", max_iter = 0,
    init_control = mpcurve_init_control(method_args = list(num_neighbors = 15)))
  expect_identical(fixed$fit$init_info$k_used, 15L)
})
