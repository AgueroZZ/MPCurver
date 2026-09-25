test_that("uniform-prior dimension selection records both greedy directions", {
  set.seed(1401)
  X <- matrix(rnorm(120), nrow = 30, ncol = 4)
  common <- list(
    X = X,
    K = 5,
    iter = 2,
    T_start = 1,
    T_end = 1,
    n_outer = 1,
    inner_iter = 1,
    max_converge_iter = 1,
    verbose = FALSE
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
    expect_equal(nrow(selection$candidates), 3L)
    expect_equal(nrow(selection$history), 1L)
    expect_equal(sort(selection$candidates$M), c(1L, 2L, 2L))
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

  direct1 <- do.call(fit_mpcurve, c(common, list(
    intrinsic_dim = 1,
    partition_prior = "fixed"
  )))
  direct2_similarity <- do.call(fit_mpcurve, c(common, list(
    intrinsic_dim = 2,
    partition_init = "similarity",
    partition_prior = "fixed",
    hard_assign_final = FALSE
  )))
  direct2_ordering <- do.call(fit_mpcurve, c(common, list(
    intrinsic_dim = 2,
    partition_init = "ordering_methods",
    partition_prior = "fixed",
    hard_assign_final = FALSE
  )))
  expected_scores <- c(
    tail(direct1$elbo_trace, 1L),
    tail(direct2_similarity$objective_history, 1L),
    tail(direct2_ordering$objective_history, 1L)
  )
  forward_rows <- forward$dimension_selection$candidates
  backward_rows <- backward$dimension_selection$candidates
  forward_rows <- forward_rows[order(forward_rows$M, forward_rows$initialization), ]
  backward_rows <- backward_rows[order(backward_rows$M, backward_rows$initialization), ]
  expected_rows <- data.frame(
    M = c(1L, 2L, 2L),
    initialization = c("single", "ordering_methods", "similarity"),
    score = expected_scores[c(1L, 3L, 2L)]
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

  continued <- do_mpcurve(forward, iter = 1L)
  expect_equal(continued$dimension_selection$selected_M,
               forward$dimension_selection$selected_M)
  expect_true(continued$dimension_selection$continued_after_selection)
  expect_equal(continued$effective_intrinsic_dim,
               continued$intrinsic_dim)

  direct2_similarity$dimension_selection <- list(
    direction = "backward",
    max_intrinsic_dim = 2L,
    selected_M = 2L,
    continued_after_selection = FALSE
  )
  continued_partition <- do_mpcurve(direct2_similarity, iter = 1L)
  expect_true(continued_partition$dimension_selection$continued_after_selection)
  expect_equal(continued_partition$dimension_selection$selected_M, 2L)
})

test_that("dimension selection accepts a one-dimensional bound", {
  set.seed(1402)
  X <- matrix(rnorm(60), nrow = 20, ncol = 3)
  for (direction in c("forward", "backward")) {
    fit <- select_mpcurve_dimension(
      X,
      max_intrinsic_dim = 1,
      direction = direction,
      K = 4,
      iter = 1,
      hard_assign_final = FALSE,
      partition_prior = "fixed",
      partition_prior_init = NULL,
      position_prior = "fixed",
      position_prior_init = rep(1 / 4, 4),
      greedy = "none",
      verbose = FALSE
    )
    expect_equal(fit$intrinsic_dim, 1L)
    expect_equal(fit$dimension_selection$selected_M, 1L)
    expect_equal(nrow(fit$dimension_selection$candidates), 1L)
    expect_equal(nrow(fit$dimension_selection$history), 0L)
  }
})

test_that("dimension selection rejects controls that change the comparison", {
  X <- matrix(0, nrow = 10, ncol = 3)
  expect_error(select_mpcurve_dimension(X, 0), "positive integer")
  expect_error(select_mpcurve_dimension(X, 2.5), "positive integer")
  expect_error(select_mpcurve_dimension(X, 2, intrinsic_dim = 2),
               "selector controls")
  expect_error(select_mpcurve_dimension(X, 2, partition_init = "similarity"),
               "selector controls")
  expect_error(select_mpcurve_dimension(X, 2, partition_prior = "adaptive"),
               "partition_prior")
  expect_error(select_mpcurve_dimension(X, 2, partition_prior_init = c(0.5, 0.5)),
               "uniform partition prior")
  expect_error(select_mpcurve_dimension(X, 2, greedy = "forward"),
               "greedy")
  expect_error(select_mpcurve_dimension(X, 2, hard_assign_final = TRUE),
               "hard_assign_final")
  expect_error(select_mpcurve_dimension(X, 2, method = c("PCA", "fiedler")),
               "method must have length one")
})
