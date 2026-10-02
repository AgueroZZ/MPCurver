test_that("uniform-prior dimension selection records both greedy directions", {
  set.seed(1401)
  X <- matrix(rnorm(120), nrow = 30, ncol = 4)
  common <- list(
    X = X,
    num_bins = 5,
    max_iter = 1,
    verbose = FALSE,
    control = mpcurve_control(anneal_start = 1, anneal_steps = 1, anneal_sweeps = 1)
  )

  forward <- do.call(select_mpcurve_dimension, c(
    common, list(max_intrinsic_dim = 2, direction = "forward")
  ))
  backward <- do.call(select_mpcurve_dimension, c(
    common, list(max_intrinsic_dim = 2, direction = "backward")
  ))

  for (fit in list(forward, backward)) {
    selection <- fit$dimension_selection
    expect_s3_class(fit, "mpcurve")
    expect_equal(selection$max_intrinsic_dim, 2L)
    expect_equal(selection$selected_M, fit$intrinsic_dim)
    expect_identical(selection$criterion, "final_soft_T1_ELBO")
    expect_false(selection$continued_after_selection)
    expect_equal(nrow(selection$candidates), 2L)
    expect_equal(nrow(selection$history), 1L)
    expect_equal(sort(selection$candidates$M), c(1L, 2L))
    expect_equal(sum(selection$candidates$selected_for_M), 2L)
    expect_true(all(selection$candidates$status == "success"))
    expect_true(all(selection$candidates$final_temperature == 1))
    expect_true(all(selection$candidates$K == 5L))
    expect_equal(selection$history$delta,
                 selection$history$candidate_score - selection$history$current_score)
    expect_equal(selection$history$accepted, selection$history$delta > 0)
    expect_equal(selection$selected_M,
                 if (selection$history$accepted) selection$history$candidate_M else
                   selection$history$current_M)
  }

  expect_equal(forward$dimension_selection$history$current_M, 1L)
  expect_equal(forward$dimension_selection$history$candidate_M, 2L)
  expect_equal(backward$dimension_selection$history$current_M, 2L)
  expect_equal(backward$dimension_selection$history$candidate_M, 1L)

  args <- common
  args$partition_prior <- "fixed"
  args$intrinsic_dim <- 1
  direct1 <- do.call(fit_mpcurve, args)
  args$intrinsic_dim <- 2
  args$init_control <- mpcurve_init_control()
  direct2_similarity <- do.call(fit_mpcurve, args)
  expected_scores <- c(
    tail(direct1$elbo_trace, 1L),
    tail(direct2_similarity$objective_history, 1L)
  )
  forward_rows <- forward$dimension_selection$candidates
  backward_rows <- backward$dimension_selection$candidates
  forward_rows <- forward_rows[order(forward_rows$M, forward_rows$initialization), ]
  backward_rows <- backward_rows[order(backward_rows$M, backward_rows$initialization), ]
  expected_rows <- data.frame(
    M = c(1L, 2L),
    initialization = c("single", "similarity"),
    score = expected_scores
  )
  expect_equal(forward_rows$M, expected_rows$M)
  expect_equal(forward_rows$initialization, expected_rows$initialization)
  expect_equal(forward_rows$score, expected_rows$score, tolerance = 1e-8)
  expect_equal(backward_rows$score, expected_rows$score, tolerance = 1e-8)

  if (forward$intrinsic_dim == 2L) {
    expect_identical(fitted_prior(forward, type = "partition")$mode, "fixed")
    expect_equal(as.numeric(fitted_prior(forward, type = "partition")$omega),
                 rep(1 / 2, 2))
  }
  if (backward$intrinsic_dim == 2L) {
    expect_identical(fitted_prior(backward, type = "partition")$mode, "fixed")
    expect_equal(as.numeric(fitted_prior(backward, type = "partition")$omega),
                 rep(1 / 2, 2))
  }

  expect_equal(summary(forward)$dimension_selection$selected_M,
               forward$intrinsic_dim)
  expect_match(paste(capture.output(print(forward)), collapse = "\n"),
               "M selection", fixed = TRUE)
  expect_match(paste(capture.output(print(summary(forward))), collapse = "\n"),
               "M selection", fixed = TRUE)

  continued <- do_mpcurve(forward, max_iter = 1L)
  expect_equal(continued$dimension_selection$selected_M,
               forward$dimension_selection$selected_M)
  expect_true(continued$dimension_selection$continued_after_selection)
  expect_equal(continued$model_intrinsic_dim,
               continued$intrinsic_dim)

  direct2_similarity$dimension_selection <- list(
    direction = "backward",
    max_intrinsic_dim = 2L,
    selected_M = 2L,
    continued_after_selection = FALSE
  )
  continued_partition <- do_mpcurve(direct2_similarity, max_iter = 1L)
  expect_true(continued_partition$dimension_selection$continued_after_selection)
  expect_equal(continued_partition$dimension_selection$selected_M, 2L)
})

test_that("dimension selection uses similarity groups and the chosen method", {
  set.seed(1403)
  X <- matrix(rnorm(120), nrow = 30, ncol = 4)
  for (method in c("PCA", "fiedler")) {
    common <- list(
      X = X, S = rep(0.8, ncol(X)), num_bins = 5, max_iter = 1,
      initial_method = method,
      control = mpcurve_control(lambda_init = 2, fix_lambda = TRUE)
    )
    direct <- lapply(1:2, function(M) {
      do.call(fit_mpcurve, c(common,
        list(intrinsic_dim = M, partition_prior = "fixed")))
    })
    expected_scores <- c(tail(direct[[1]]$elbo_trace, 1L),
                         tail(direct[[2]]$objective_history, 1L))
    for (direction in c("forward", "backward")) {
      fit <- do.call(select_mpcurve_dimension,
        c(common, list(max_intrinsic_dim = 2, direction = direction)))
      rows <- fit$dimension_selection$candidates
      rows <- rows[order(rows$M), ]
      expect_equal(rows$M, 1:2)
      expect_equal(rows$initialization, c("single", "similarity"))
      expect_equal(rows$score, expected_scores, tolerance = 1e-8)
      expect_true(all(rows$status == "success"))
      expect_true(all(rows$selected_for_M))
    }
  }
})

test_that("dimension selection accepts a one-dimensional bound", {
  set.seed(1402)
  X <- matrix(rnorm(60), nrow = 20, ncol = 3)
  for (direction in c("forward", "backward")) {
    fit <- select_mpcurve_dimension(
  X,
  max_intrinsic_dim = 1,
    direction = direction,
    num_bins = 4,
    max_iter = 1,
    position_prior = "fixed",
    verbose = FALSE,
    control = mpcurve_control(partition_prior_weights = NULL, position_prior_weights = rep(1 / 4, 4))
  )
    expect_equal(fit$intrinsic_dim, 1L)
    expect_equal(fit$dimension_selection$selected_M, 1L)
    expect_equal(nrow(fit$dimension_selection$candidates), 1L)
    expect_equal(nrow(fit$dimension_selection$history), 0L)
  }
})

test_that("dimension selection rejects controls that change the comparison", {
  X <- matrix(0, 10, 3)
  expect_error(select_mpcurve_dimension(X, 0), "integer")
  expect_error(select_mpcurve_dimension(X, 2.5), "integer")
  expect_error(select_mpcurve_dimension(X, 2, intrinsic_dim = 2), "unused argument")
  expect_error(select_mpcurve_dimension(X, 2, partition_prior = "adaptive"), "unused argument")
  expect_error(select_mpcurve_dimension(X, 2,
    control = mpcurve_control(partition_prior_weights = c(.5, .5))), "uniform partition prior")
  expect_error(select_mpcurve_dimension(X, 2,
    control = mpcurve_control(hard_assign_final = TRUE)), "hard_assign_final")
  expect_error(select_mpcurve_dimension(X, 2, initial_method = c("PCA", "fiedler")), "exactly one")
  expect_error(select_mpcurve_dimension(X, 2, greedy = "forward"), "unused argument")
})
