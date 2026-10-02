test_that("a scalar sinusoid frequency specifies the actual simulated frequency", {
  t <- seq(0, 1, length.out = 80)
  sim <- simulate_intrinsic_trajectories(n = 80, d_signal = 8, d_noise = 0,
    noise_sd = 0, seed = 1, control = list(sinusoid_freq = 4,
      signal_range = c(1, 1), latent_positions = matrix(t, ncol = 1)))
  basis <- qr(cbind(sin(4 * pi * t), cos(4 * pi * t)))
  expect_lt(max(abs(qr.resid(basis, sim$X))), 1e-12)
})

test_that("singleton GP blocks retain sample by feature dimensions", {
  for (permute in c(FALSE, TRUE)) {
    sim <- simulate_two_order_gp_dataset(n = 10, d = 2, seed = 1,
      control = list(permute_rows_block2 = permute))
    expect_identical(dim(sim$X), c(10L, 2L))
    expect_true(all(is.finite(sim$X)))
    expect_setequal(sim$true_group, 1:2)
  }
})

test_that("automatic grids respect the RW order across dimensions", {
  set.seed(10)
  X <- matrix(rnorm(40), 10, 4)
  for (q in 1:3) for (M in list(1L, 2L, "auto")) {
    fit <- fit_mpcurve(X, intrinsic_dim = M, max_iter = 0,
      control = mpcurve_control(rw_order = q),
      init_control = mpcurve_init_control(max_intrinsic_dim = 2))
    expect_gte(fit$K, q + 1L)
  }
  expect_error(fit_mpcurve(X[1:2, ], max_iter = 0), "at least rw_order")
  expect_error(fit_mpcurve(X, num_bins = 2, max_iter = 0), "at least rw_order")
})

test_that("plot labels override defaults and narrow pseudotimes retain absolute colors", {
  set.seed(44)
  X <- matrix(rnorm(120), 40, 3)
  single <- fit_mpcurve(X, S = rep(100, 3), num_bins = 5, max_iter = 2,
    position_prior = "fixed")
  multi <- fit_mpcurve(X, intrinsic_dim = 2, num_bins = 5, max_iter = 2)
  output <- tempfile(fileext = ".pdf")
  grDevices::pdf(output)
  on.exit({grDevices::dev.off(); unlink(output)}, add = TRUE)
  calls <- list()
  original <- MPCurver:::.mpcurve_plot_call
  local_mocked_bindings(.mpcurve_plot_call = function(fun, defaults, ...) {
    args <- defaults
    args[names(list(...))] <- list(...)
    calls[[length(calls) + 1L]] <<- args
    original(fun, defaults, ...)
  }, .package = "MPCurver")
  for (fit in list(single, multi)) for (kind in c("scatterplot", "mu", "elbo")) {
    for (dims in list(1L, 1:2)) {
      calls <- list()
      expect_no_warning(plot(fit, plot_type = kind, dims = dims,
        main = "Custom title", xlab = "Custom x", ylab = "Custom y"))
      expect_true(length(calls) > 0L)
      expect_true(all(vapply(calls, function(a) identical(a$main, "Custom title"), logical(1))))
      expect_true(all(vapply(calls, function(a) identical(a$xlab, "Custom x"), logical(1))))
    }
  }
  pal <- grDevices::rainbow(101)
  for (dims in list(1L, 1:2)) {
    calls <- list()
    plot(single, dims = dims, pal = pal, add_legend = FALSE)
    colors <- calls[[1]]$col
    positions <- single$gamma %*% seq(0, 1, length.out = single$K)
    expect_identical(colors, pal[1L + floor(as.numeric(positions) * 100)])
    expect_lte(length(unique(colors)), 2L)
  }
})

test_that("partition fitting retains soft weights and counts only CAVI sweeps", {
  expect_error(mpcurve_control(hard_assign_final = TRUE), "unused argument")
  set.seed(9)
  fit <- fit_mpcurve(matrix(rnorm(160), 40, 4), intrinsic_dim = 2,
    num_bins = 5, max_iter = 0)
  expect_equal(fit$iter, 0L)
  continued <- do_mpcurve(fit, max_iter = 2, tol = 0)
  expect_equal(continued$iter, 2L)
  weights <- continued$partition$pi_weights
  expect_true(any(weights > 0 & weights < 1))
  expect_equal(rowSums(weights), rep(1, nrow(weights)), ignore_attr = TRUE)
})
