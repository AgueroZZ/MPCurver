test_that("simulate_dual_trajectory supports multiple trajectory families", {
  families <- c("linear", "monotone", "quadratic", "sinusoidal")

  for (fam in families) {
    sim <- simulate_dual_trajectory(
      n = 40,
      d1 = 3,
      d2 = 2,
      d_noise = 1,
      noise_sd = 0.05,
      seed = 10,
      trajectory_family = fam
    )

    expect_equal(dim(sim$X), c(40, 6))
    expect_equal(length(sim$true_assign), 6)
    expect_equal(sim$trajectory_family, c(fam, fam))
  }
})

test_that("simulate_intrinsic_trajectories supports multiple intrinsic dimensions", {
  latent <- cbind(
    seq(0, 1, length.out = 30),
    seq(1, 0, length.out = 30),
    rep(0.5, 30)
  )
  sim <- simulate_intrinsic_trajectories(
    n = 30,
    d_signal = c(2, 3, 1),
    d_noise = 2,
    noise_sd = 0.05,
    seed = 11,
    trajectory_family = c("linear", "quadratic", "monotone"),
    control = list(latent_positions = latent)
  )

  expect_equal(dim(sim$X), c(30, 8))
  expect_equal(dim(sim$latent_positions), c(30, 3))
  expect_equal(sim$ordering_labels, c("A", "B", "C"))
  expect_equal(as.integer(sort(table(sim$true_assign)[c("A", "B", "C", "noise")])),
               c(1, 2, 2, 3))
  expect_equal(sim$trajectory_family, c("linear", "quadratic", "monotone"))
  expect_equal(sim$latent_positions, latent)
})
