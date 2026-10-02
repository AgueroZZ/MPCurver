test_that("public namespace exports the curated mpcurve-first API", {
  expect_setequal(
    getNamespaceExports("MPCurver"),
    c(
      "PCA_ordering",
      "do_mpcurve",
      "fiedler_ordering",
      "fit_mpcurve",
      "mpcurve_control",
      "mpcurve_continue_control",
      "mpcurve_init_control",
      "select_mpcurve_dimension",
      "fitted_prior",
      "fitted_positions",
      "fitted_trajectories",
      "fitted_assignments",
      "isomap_ordering",
      "pcurve_ordering",
      "simulate_mpcurve",
      "simulate_dual_trajectory",
      "simulate_intrinsic_trajectories",
      "simulate_spiral2d",
      "simulate_swiss_roll_1d_2d",
      "simulate_two_order_gp_dataset",
      "tSNE_ordering"
    )
  )

  namespace_path <- testthat::test_path("..", "..", "NAMESPACE")
  if (!file.exists(namespace_path)) {
    namespace_path <- system.file("NAMESPACE", package = "MPCurver")
  }
  namespace_text <- paste(readLines(namespace_path), collapse = "\n")

  expect_false(grepl("S3method\\(print,cavi\\)", namespace_text, fixed = FALSE))
  expect_false(grepl("S3method\\(summary,cavi\\)", namespace_text, fixed = FALSE))
  expect_false(grepl("S3method\\(plot,cavi\\)", namespace_text, fixed = FALSE))
  expect_true(is.function(utils::getS3method("print", "mpcurve")))
  expect_true(is.function(utils::getS3method("summary", "mpcurve")))
  expect_true(is.function(utils::getS3method("plot", "mpcurve")))
})

test_that("fit_mpcurve has one initialization method and twelve public arguments", {
  expect_identical(names(formals(fit_mpcurve)), c(
    "X", "S", "num_bins", "intrinsic_dim", "initial_method", "position_prior",
    "partition_prior", "max_iter", "tol", "verbose", "init_control", "control"
  ))
  X <- matrix(0, 10, 4)
  for (M in list(1L, 2L, "auto")) {
    expect_error(fit_mpcurve(X, intrinsic_dim = M,
                            initial_method = c("PCA", "isomap")), "exactly one")
    for (argument in c("algorithm", "greedy", "num_cores", "tol_outer",
                       "assignment_prior", "ordering_alpha", "freeze_feature",
                       "freeze_unused_ordering", "drop_unused_ordering")) {
      expect_error(do.call(fit_mpcurve, c(list(X = X, intrinsic_dim = M),
                                         setNames(list(NULL), argument))), "unused argument")
    }
  }
})

test_that("feature grouping is the only automatic partition initialization path", {
  expect_false("partition_init" %in% names(formals(mpcurve_init_control)))
  expect_false("partition_init" %in% names(mpcurve_init_control()))
  expect_error(mpcurve_init_control(partition_init = "similarity"), "unused argument")
  expect_error(mpcurve_init_control(partition_init = "ordering_methods"), "unused argument")
  X <- matrix(0, 10, 4)
  for (M in list(1L, 2L, "auto")) {
    expect_error(fit_mpcurve(X, intrinsic_dim = M,
      init_control = list(partition_init = "ordering_methods")), "Unknown init_control")
  }
  expect_error(select_mpcurve_dimension(X, max_intrinsic_dim = 2,
    init_control = list(partition_init = "ordering_methods")), "Unknown init_control")
})

test_that("print.mpcurve is shorter than summary for single and partition fits", {
  sim_single <- simulate_mpcurve(
    n = 60,
    d = 8,
    num_bins = 5,
    seed = 101,
    control = list(rw_order = 2)
  )
  fit_single <- fit_mpcurve(
    sim_single$X,
    initial_method = "PCA",
    num_bins = 5,
    max_iter = 4,
    tol = 0,
    verbose = FALSE
  )

  single_print <- paste(capture.output(print(fit_single)), collapse = "\n")
  single_summary <- paste(capture.output(summary(fit_single)), collapse = "\n")
  single_summary_obj <- summary(fit_single)
  expect_lt(length(strsplit(single_print, "\n", fixed = TRUE)[[1]]),
            length(strsplit(single_summary, "\n", fixed = TRUE)[[1]]))
  expect_true(all(c("priors", "converged", "underlying") %in% names(single_summary_obj)))
  expect_snapshot_output(print(fit_single))

  sim_partition <- simulate_dual_trajectory(
    n = 60,
    d1 = 4,
    d2 = 4,
    d_noise = 0,
    noise_sd = 0.15,
    seed = 102
  )
  fit_partition <- suppressWarnings(fit_mpcurve(
    sim_partition$X,
    intrinsic_dim = 2,
    initial_method = "PCA",
    num_bins = 6,
    max_iter = 2L,
    verbose = FALSE,
    control = mpcurve_control(anneal_steps = 2L, anneal_sweeps = 1L)
  ))

  partition_print <- paste(capture.output(print(fit_partition)), collapse = "\n")
  partition_summary <- paste(capture.output(summary(fit_partition)), collapse = "\n")
  partition_summary_obj <- summary(fit_partition)
  expect_lt(length(strsplit(partition_print, "\n", fixed = TRUE)[[1]]),
            length(strsplit(partition_summary, "\n", fixed = TRUE)[[1]]))
  expect_true(all(c(
    "priors",
    "partition",
    "objective_history",
    "convergence_info",
    "sigma2",
    "measurement_sd",
    "variational_family",
    "greedy_selection"
  ) %in% names(partition_summary_obj)))
  expect_false("underlying" %in% names(partition_summary_obj))
  expect_identical(partition_summary_obj$variational_family, "structured")
  expect_equal(partition_summary_obj$sigma2, fit_partition$params$sigma2)
  expect_equal(partition_summary_obj$active_intrinsic_dim,
               fit_partition$model_intrinsic_dim)
  expect_equal(partition_summary_obj$displayed_intrinsic_dim,
               fit_partition$model_intrinsic_dim)
  expect_match(partition_summary, "structural VI", fixed = TRUE)
  expect_match(partition_summary, "Shared sigma2 range", fixed = TRUE)
  expect_match(partition_summary, "fixed-M structural ELBO", fixed = TRUE)

  position_prior <- fitted_prior(fit_partition, type = "position")
  expect_named(position_prior, names(fit_partition$gamma))
  expect_true(all(vapply(position_prior, length, integer(1)) == fit_partition$K))
  expect_equal(
    fitted_prior(fit_partition, type = "position", ordering = names(position_prior)[1]),
    position_prior[[1]]
  )
  partition_prior <- fitted_prior(fit_partition, type = "partition")
  expect_named(partition_prior$omega, names(fit_partition$gamma))
  expect_equal(sum(partition_prior$omega), 1, tolerance = 1e-12)
  expect_equal(
    fitted_prior(fit_partition, type = "partition", ordering = 1L),
    unname(partition_prior$omega[1])
  )
  expect_snapshot_output(print(fit_partition))
})
