test_that("invalid initializer settings error before computation in every dimension", {
  set.seed(614)
  X <- matrix(rnorm(200), 40, 5)
  cases <- list(
    list("fiedler", list(num_neighbors = 1), "num_neighbors"),
    list("fiedler", list(num_neighbors = 40), "num_neighbors"),
    list("isomap", list(num_neighbors = NULL), "num_neighbors"),
    list("isomap", list(component = 2), "embedding_dims"),
    list("tSNE", list(perplexity = 20), "perplexity"),
    list("tSNE", list(max_iter = 1.5), "max_iter"),
    list("tSNE", list(seed = c(1, 2)), "seed"),
    list("pcurve", list(smoother = "invalid"), "arg"),
    list("pcurve", list(tol = 0), "tol"),
    list("PCA", list(scale = NA), "scale")
  )
  for (M in list(1L, 2L, "auto")) {
    for (case in cases) {
      expect_error(fit_mpcurve(X, num_bins = 5, intrinsic_dim = M,
        initial_method = case[[1]], max_iter = 0,
        init_control = list(method_args = case[[2]])), case[[3]])
    }
    expect_error(fit_mpcurve(X, intrinsic_dim = M,
      init_control = list(pca_components = 1.9)), "pca_components")
  }
  expect_error(mpcurve_init_control(on_failure = "silent"), "arg")
  expect_error(MPCurver:::.cavi_resolve_similarity_methods(
    pca_components = 1.9, M = 2L), "pca_components")
  expect_error(MPCurver:::.cavi_similarity_subset_fit(
    X[, 1, drop = FALSE], K = 5, method = "PCA", pca_component = 2),
    "available PCs")
})

test_that("ordering failures warn and preserve PCA fallback provenance", {
  set.seed(614)
  X <- matrix(rnorm(200), 40, 5)
  local_mocked_bindings(fiedler_ordering = function(X, num_neighbors = NULL,
                                                  control = NULL) {
    stop("eigensolver failed to converge")
  }, .package = "MPCurver")
  for (M in list(1L, 2L, "auto")) {
    messages <- character()
    fit <- withCallingHandlers(
      fit_mpcurve(X, num_bins = 5, intrinsic_dim = M, max_iter = 0,
        initial_method = "fiedler", init_control = list(max_intrinsic_dim = 2)),
      warning = function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    expect_true(length(messages) >= 1L)
    expect_true(all(grepl("eigensolver failed to converge.*Falling back to PCA", messages)))
    expect_true(all(grepl("on_failure = 'error'", messages, fixed = TRUE)))
    info <- fit$fit$init_info
    if (fit$intrinsic_dim == 1L) info <- list(info)
    for (entry in info) {
      expect_identical(entry$method_requested, "fiedler")
      expect_identical(entry$method_used, "PCA")
      expect_true(entry$fallback)
      expect_identical(entry$fallback_reason, "eigensolver failed to converge")
    }
    expect_no_warning(continued <- do_mpcurve(fit, max_iter = 1))
    expect_identical(continued$fit$init_info, fit$fit$init_info)
    expect_error(fit_mpcurve(X, num_bins = 5, intrinsic_dim = M, max_iter = 0,
      initial_method = "fiedler", init_control = list(on_failure = "error")),
      "fiedler initialization failed: eigensolver failed to converge")
  }
})

test_that("PCA failure stops and downstream fitting errors do not trigger fallback", {
  X <- matrix(seq_len(120), 30, 4)
  local_mocked_bindings(
    fiedler_ordering = function(X, num_neighbors = NULL, control = NULL) stop("graph failure"),
    PCA_ordering = function(X, component = 1L, center = TRUE, scale = FALSE,
                            scale01 = TRUE) stop("SVD failure"),
    .package = "MPCurver")
  expect_warning(expect_error(fit_mpcurve(X, num_bins = 5,
    initial_method = "fiedler", max_iter = 0), "PCA fallback failed: SVD failure"),
    "Falling back to PCA")
  expect_error(fit_mpcurve(X, num_bins = 5, max_iter = 0),
    "PCA initialization failed: SVD failure")
})

test_that("model fitting errors propagate without another initializer", {
  set.seed(6)
  X <- matrix(rnorm(120), 30, 4)
  local_mocked_bindings(.cavi_build_from_ordering = function(...) {
    stop("model fitting failure")
  }, .package = "MPCurver")
  expect_error(MPCurver:::.cavi_similarity_subset_fit(X, K = 5, method = "PCA"),
    "model fitting failure")
})
