test_that("fiedler_ordering auto-increases the default k until the graph is connected", {
  centers <- rbind(
    c(-8, 0),
    c(0, 0),
    c(8, 0)
  )
  X <- do.call(
    rbind,
    lapply(seq_len(nrow(centers)), function(g) {
      matrix(
        stats::rnorm(40, mean = rep(centers[g, ], each = 20), sd = 0.15),
        ncol = 2,
        byrow = FALSE
      )
    })
  )

  fit_default <- expect_no_warning(fiedler_ordering(X))
  expect_equal(fit_default$n_components, 1L)
  expect_equal(length(fit_default$keep_idx), nrow(X))
  expect_equal(sum(!is.na(fit_default$t)), nrow(X))
  expect_gt(fit_default$k_used, 15)

  expect_warning(
    fit_explicit <- fiedler_ordering(X, num_neighbors = 15),
    "connected components"
  )
  expect_gt(fit_explicit$n_components, 1L)
  expect_lt(length(fit_explicit$keep_idx), nrow(X))
  expect_equal(fit_explicit$k_used, 15)
})

test_that("random starts bin a random sample permutation like any ordering score", {
  initialize_ordering <- ns_fn(".cavi_initialize_ordering")
  bin_ordering <- ns_fn(".cavi_bin_ordering")
  scale01 <- ns_fn(".scale01")
  set.seed(5)
  X <- matrix(rnorm(43 * 3), 43, 3)

  set.seed(1)
  init <- initialize_ordering(X, K = 6, method = "random", discretization = "quantile")
  set.seed(1)
  reference <- bin_ordering(X, scale01(sample.int(nrow(X))), K = 6,
                            discretization = "quantile")
  expect_identical(init$cluster_rank, reference$cluster_rank)
  expect_identical(init$params, reference[c("pi", "mu", "sigma2")])
  # Quantile bins of a permutation hold floor(43/6) or ceiling(43/6) samples.
  expect_true(all(tabulate(init$cluster_rank, 6) %in% c(7L, 8L)))
  info <- attr(init$ordering, "init_info")
  expect_identical(info$method_used, "random")
  expect_false(info$fallback)

  set.seed(1)
  expect_identical(initialize_ordering(X, K = 6, method = "random")$cluster_rank,
                   init$cluster_rank)
  set.seed(2)
  expect_false(identical(initialize_ordering(X, K = 6, method = "random")$cluster_rank,
                         init$cluster_rank))
})

test_that("random starts fit single- and multi-ordering models", {
  sim <- simulate_mpcurve(n = 60, d = 5, num_bins = 6, seed = 1)
  set.seed(3)
  start <- fit_mpcurve(sim$X, initial_method = "random", num_bins = 6, max_iter = 0)
  expect_equal(colSums(start$gamma), rep(10, 6))
  expect_identical(start$fit$init_info$method_used, "random")

  set.seed(3)
  fit <- fit_mpcurve(sim$X, initial_method = "random", num_bins = 6, max_iter = 20)
  expect_true(all(is.finite(fit$fit$elbo_trace)))
  expect_true(all(diff(fit$fit$elbo_trace) >= -1e-8))

  dual <- simulate_dual_trajectory(n = 60, d1 = 3, d2 = 3, d_noise = 0,
                                   noise_sd = 0.1, seed = 2)
  set.seed(4)
  fit2 <- fit_mpcurve(dual$X, intrinsic_dim = 2, initial_method = "random",
                      num_bins = 6, max_iter = 3)
  expect_identical(fit2$intrinsic_dim, 2L)
  expect_identical(unname(vapply(fit2$init_info, `[[`, character(1), "method_used")),
                   c("random", "random"))
  expect_true(all(is.finite(fit2$objective_history)))
  expect_true(all(diff(fit2$objective_history) >= -1e-8))
})
