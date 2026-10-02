test_that("effective ordering count ignores numerical prior dust", {
  count_effective <- getFromNamespace(
    ".mpcurve_effective_intrinsic_dim", "MPCurver"
  )
  priors <- list(partition = list(
    mode = "adaptive",
    omega = c(A = 0.5, B = 0.5, C = 8.89318162514244e-323, D = 0)
  ))

  expect_identical(count_effective(priors, 4L, 1e-8, 100 * priors$partition$omega), 2L)
  expect_identical(count_effective(priors, 4L, 0, 100 * priors$partition$omega), 3L)
  priors$partition$mode <- "fixed"
  expect_identical(count_effective(priors, 4L, 1e-8, 100 * priors$partition$omega), 4L)
})

test_that("explicit fitted and continued models retain the specified dimension", {
  set.seed(117)
  X <- matrix(rnorm(144), nrow = 24, ncol = 6)
  common <- list(X = X, intrinsic_dim = 2L, num_bins = 5L,
                 initial_method = "PCA", max_iter = 1L,
                 init_control = mpcurve_init_control(),
                 control = mpcurve_control(anneal_start = 1, anneal_steps = 1L,
                                          effective_count_tol = 0.05), verbose = FALSE)

  adaptive <- suppressWarnings(do.call(fit_mpcurve, common))
  expect_identical(adaptive$intrinsic_dim, 2L)
  expect_identical(adaptive$effective_count_tol, 0.05)
  expect_identical(adaptive$fit$control$effective_count_tol, 0.05)
  expect_identical(adaptive$active_intrinsic_dim, 2L)
  expect_identical(summary(adaptive)$intrinsic_dim,
                   adaptive$intrinsic_dim)
  expect_match(paste(capture.output(print(adaptive)), collapse = "\n"),
               "Intrinsic dim", fixed = TRUE)

  continued <- suppressWarnings(do_mpcurve(adaptive, max_iter = 1L))
  expect_identical(continued$effective_count_tol, 0.05)
  expect_identical(continued$intrinsic_dim, 2L)

  fixed <- suppressWarnings(do.call(
    fit_mpcurve, c(common, list(partition_prior = "fixed"))
  ))
  expect_identical(fixed$intrinsic_dim, 2L)
  expect_identical(summary(fixed)$intrinsic_dim, 2L)

  single <- suppressWarnings(fit_mpcurve(
    X,
    num_bins = 5L,
    max_iter = 1L
  ))
  expect_identical(single$intrinsic_dim, 1L)
  expect_identical(summary(single)$intrinsic_dim, 1L)
})

test_that("effective ordering threshold is validated before fitting", {
  X <- matrix(rnorm(24), nrow = 6, ncol = 4)
  expect_error(
    fit_mpcurve(
    X,
    intrinsic_dim = 2L,
    control = mpcurve_control(effective_count_tol = -1)
  ),
    "effective_count_tol"
  )
  expect_error(
    fit_mpcurve(
    X,
    intrinsic_dim = 2L,
    control = mpcurve_control(effective_count_tol = 2)
  ),
    "effective_count_tol"
  )
  expect_error(
    fit_mpcurve(
    X,
    intrinsic_dim = 2L,
    control = mpcurve_control(effective_count_tol = NA_real_)
  ),
    "effective_count_tol"
  )
})

test_that("legacy reported dimensions remain independent of retained slots", {
  set.seed(427)
  X <- matrix(rnorm(144), nrow = 24, ncol = 6)
  fit <- suppressWarnings(fit_mpcurve(
    X, intrinsic_dim = 2L, num_bins = 5L, max_iter = 1L,
    init_control = mpcurve_init_control(),
    control = mpcurve_control(anneal_steps = 1L, anneal_start = 1,
                             effective_count_tol = 1.2)
  ))
  # Exercise conversion with one supported ordering while retaining two slots.
  # Only prior metadata is changed; all conditional trajectories remain present.
  raw <- fit$fit
  raw$control$intrinsic_dim_semantics <- NULL
  raw$priors <- fit$priors
  raw$pi_weights[,] <- rep(c(0.9, 0.1), each = nrow(raw$pi_weights))
  raw$priors$partition$omega[] <- c(0.9, 0.1)
  reduced <- getFromNamespace("as_mpcurve", "MPCurver")(raw)
  expect_identical(reduced$intrinsic_dim, 1L)
  expect_identical(reduced$model_intrinsic_dim, 2L)
  expect_null(reduced$effective_intrinsic_dim)
  expect_length(reduced$conditional_posterior$mean, 2L)
  expect_equal(dim(reduced$lambda_mat), c(ncol(X), 2L))
  expect_identical(summary(reduced)$intrinsic_dim, 1L)
  expect_identical(summary(reduced)$model_intrinsic_dim, 2L)
  expect_null(summary(reduced)$effective_intrinsic_dim)

  pdf_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  on.exit({ grDevices::dev.off(); unlink(pdf_file) }, add = TRUE)
  expect_no_error(plot(reduced, plot_type = "mu", dims = 1L))
  expect_no_error(plot(reduced, plot_type = "scatterplot", dims = c(1L, 2L)))

  # A saved object from the old schema reports the same estimate and accesses
  # all trajectories. Continuation must produce the same numerical state.
  legacy <- reduced
  legacy$model_intrinsic_dim <- NULL
  legacy$effective_intrinsic_dim <- reduced$intrinsic_dim
  legacy$intrinsic_dim <- 2L
  expect_identical(summary(legacy)$intrinsic_dim, 1L)
  expect_identical(summary(legacy)$model_intrinsic_dim, 2L)
  expect_no_error(plot(legacy, plot_type = "mu", dims = 1L))
  new_continued <- suppressWarnings(do_mpcurve(reduced, max_iter = 1L))
  old_continued <- suppressWarnings(do_mpcurve(legacy, max_iter = 1L))
  expect_identical(new_continued$fit, old_continued$fit)
  expect_identical(new_continued$model_intrinsic_dim, 2L)
  expect_length(new_continued$conditional_posterior$mean, 2L)
  expect_identical(new_continued$intrinsic_dim, as.integer(sum(colSums(new_continued$partition$pi_weights) > 1.2)))
  expect_null(old_continued$effective_intrinsic_dim)
})

test_that("small ordering proportions count according to total feature mass", {
  count_effective <- getFromNamespace(".mpcurve_effective_intrinsic_dim", "MPCurver")
  omega <- c(A = 1 - 1e-13, B = 1e-13)
  priors <- list(partition = list(mode = "adaptive", omega = omega))
  # Same proportion, different feature counts; no huge matrix is needed.
  expect_identical(count_effective(priors, 2L, 1e-8, 100 * omega), 1L)
  expect_identical(count_effective(priors, 2L, 1e-8, 1e8 * omega), 2L)
  expect_identical(count_effective(priors, 2L, 1e-8, c(10, 1e-8)), 1L)
  expect_identical(count_effective(priors, 2L, 1e-8, c(10, 1.001e-8)), 2L)
  expect_error(count_effective(priors, 2L, 1e-8, c(1, NA)), "expected feature counts")
  expect_error(count_effective(priors, 2L, 1e-8, c(1, -1)), "expected feature counts")
})

test_that("adaptive EB weights agree with total soft assignment counts", {
  update <- getFromNamespace(".structural_partition_assignment_info", "MPCurver")
  weights <- rbind(c(0.9, 0.1), c(0.25, 0.75), c(0.8, 0.2))
  info <- update(weights, 1, list(partition_prior = "adaptive"))
  expect_equal(nrow(weights) * info$omega, colSums(weights))
  expect_identical(mpcurve_control()$effective_count_tol, 1e-8)
  expect_error(mpcurve_control(effective_weight_tol = 1e-12), "unused argument")
})

test_that("saved fraction thresholds retain their meaning on continuation", {
  count_tol <- getFromNamespace(".mpcurve_effective_count_tol", "MPCurver")
  old <- list(control = list(effective_weight_tol = 0.05))
  expect_identical(count_tol(old, 2L, 100L), 5)
  old$control$effective_count_tol <- 1e-8
  expect_identical(count_tol(old, 2L, 100L), 1e-8)
})


test_that("draft ordering-count schemas retain their model dimension", {
  dimension <- getFromNamespace(".mpcurve_model_intrinsic_dim", "MPCurver")
  draft <- structure(list(num_orderings = 3L, intrinsic_dim = 1L),
                     class = "summary.mpcurve")
  expect_identical(dimension(draft), 3L)
  draft$model_intrinsic_dim <- 4L
  expect_identical(dimension(draft), 4L)
})
