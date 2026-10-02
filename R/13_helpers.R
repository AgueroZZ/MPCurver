# Merge advanced settings without accepting misspellings or unnamed options.
.mpcurve_named_settings <- function(control, defaults, name = "control") {
  if (is.null(control)) return(defaults)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || any(!nzchar(names(control))) || anyDuplicated(names(control))))) {
    stop(name, " must be a list with unique named settings.", call. = FALSE)
  }
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) stop("Unknown ", name, " option(s): ",
                           paste(unknown, collapse = ", "), ".", call. = FALSE)
  defaults[names(control)] <- control
  defaults
}

.mpcurve_ordering_options <- function(method, control) {
  defaults <- switch(method,
    fiedler = list(weight = "binary", sigma = NULL, keep = "giant", return_full = TRUE),
    tSNE = list(embedding_dims = 1L, scale01 = TRUE, orient_by_pc1 = TRUE),
    pcurve = list(stretch = 2, approx_points = FALSE, scale01 = TRUE, orient_by_pc1 = TRUE),
    isomap = list(embedding_dims = 1L, num_landmarks = NULL,
      landmark_method = "random", keep = "giant", return_full = TRUE,
      orient_by_pc1 = TRUE, scale01 = TRUE))
  opts <- .mpcurve_named_settings(control, defaults, paste0(method, " control"))
  for (nm in intersect(names(opts), c("scale01", "orient_by_pc1", "return_full"))) {
    opts[[nm]] <- .mpcurve_logical_setting(opts[[nm]], nm)
  }
  if ("keep" %in% names(opts)) opts$keep <- match.arg(opts$keep, c("giant", "all"))
  if ("embedding_dims" %in% names(opts)) {
    opts$embedding_dims <- .mpcurve_integer_setting(opts$embedding_dims, "embedding_dims")
    if (method == "tSNE" && opts$embedding_dims > 3L) {
      stop("tSNE embedding_dims must be at most 3.", call. = FALSE)
    }
  }
  if (method == "fiedler") {
    opts$weight <- match.arg(opts$weight, c("binary", "rbf", "inv"))
    if (!is.null(opts$sigma)) {
      .mpcurve_positive_setting(opts$sigma, "sigma")
      if (length(opts$sigma) != 1L) stop("sigma must be a scalar.", call. = FALSE)
    }
  }
  if (method == "isomap") {
    if (!is.null(opts$num_landmarks)) opts$num_landmarks <-
      .mpcurve_integer_setting(opts$num_landmarks, "num_landmarks", 2L)
    opts$landmark_method <- match.arg(opts$landmark_method, c("random", "kmeans"))
  }
  if (method == "pcurve") {
    .mpcurve_positive_setting(opts$stretch, "stretch", allow_zero = TRUE)
    if (length(opts$stretch) != 1L) stop("stretch must be a scalar.", call. = FALSE)
    if (!identical(opts$approx_points, FALSE)) opts$approx_points <-
      .mpcurve_integer_setting(opts$approx_points, "approx_points", 2L)
  }
  opts
}

# Check deterministic argument errors before any initialization fallback.
.mpcurve_validate_ordering_args <- function(X, method, method_args = list()) {
  helper <- switch(method, PCA = PCA_ordering, fiedler = fiedler_ordering,
    pcurve = pcurve_ordering, tSNE = tSNE_ordering, isomap = isomap_ordering,
    random = function() NULL)
  if (!is.list(method_args) || (length(method_args) &&
      (is.null(names(method_args)) || any(!nzchar(names(method_args))) ||
       anyDuplicated(names(method_args))))) {
    stop("method_args must contain uniquely named ordering-helper arguments.", call. = FALSE)
  }
  unknown <- setdiff(names(method_args), setdiff(names(formals(helper)), "X"))
  if (length(unknown)) stop("Unknown initialization argument(s) for ", method, ": ",
    paste(unknown, collapse = ", "), ".", call. = FALSE)
  if (method == "random") return(invisible(NULL))
  defaults <- switch(method,
    PCA = list(component = 1L, center = TRUE, scale = FALSE, scale01 = TRUE),
    fiedler = list(num_neighbors = NULL),
    tSNE = list(component = 1L, perplexity = 10, max_iter = 500L, seed = NULL),
    pcurve = list(smoother = "smooth_spline", max_iter = 10L, tol = 0.001),
    isomap = list(num_neighbors = 15L, component = 1L, seed = NULL))
  args <- utils::modifyList(defaults, method_args, keep.null = TRUE)
  opts <- if (method != "PCA") .mpcurve_ordering_options(method, args$control) else list()
  for (nm in intersect(names(args), c("component", "max_iter"))) {
    .mpcurve_integer_setting(args[[nm]], nm)
  }
  if (!is.null(args$seed)) .mpcurve_integer_setting(args$seed, "seed", -.Machine$integer.max)
  if (method == "PCA") {
    for (nm in c("center", "scale", "scale01")) .mpcurve_logical_setting(args[[nm]], nm)
    if (args$component > min(dim(X))) stop("component exceeds available PCs (",
      min(dim(X)), ").", call. = FALSE)
  }
  if (method %in% c("fiedler", "isomap")) {
    if (nrow(X) < 3L) stop("Ordering initialization requires at least 3 samples.", call. = FALSE)
    if (method == "isomap" || !is.null(args$num_neighbors)) {
      .mpcurve_integer_setting(args$num_neighbors, "num_neighbors", 2L)
      if (args$num_neighbors >= nrow(X)) stop("num_neighbors must be less than nrow(X).", call. = FALSE)
    }
  }
  if (method %in% c("tSNE", "isomap") && args$component > opts$embedding_dims) {
    stop("component must not exceed embedding_dims.", call. = FALSE)
  }
  if (method == "tSNE") {
    .mpcurve_positive_setting(args$perplexity, "perplexity")
    if (length(args$perplexity) != 1L || args$perplexity >= (nrow(X) - 1) / 3) {
      stop("perplexity must be a scalar less than (nrow(X)-1)/3.", call. = FALSE)
    }
  }
  if (method == "pcurve") {
    match.arg(args$smoother, c("smooth_spline", "lowess", "periodic_lowess"))
    .mpcurve_positive_setting(args$tol, "tol")
    if (length(args$tol) != 1L) stop("tol must be a scalar.", call. = FALSE)
  }
  invisible(NULL)
}

#' Fiedler ordering from a nearest-neighbor graph
#'
#' Order samples by the second graph-Laplacian eigenvector.
#' @param X Numeric sample-by-feature matrix.
#' @param num_neighbors Number of nearest neighbors, between 2 and `nrow(X)-1`.
#'   `NULL` starts at `min(15, nrow(X)-1)` and increases the count until the
#'   graph is connected. An explicit count is held fixed.
#' @param control Named advanced settings list. Available settings:
#'   `weight = "binary"` (`"rbf"` or `"inv"` also supported);
#'   `sigma = NULL` (RBF distance bandwidth; defaults to median edge distance);
#'   `keep = "giant"` (`"all"` also supported); `return_full = TRUE`.
#'   Giant-component handling retains the largest connected component;
#'   full output has `NA` for excluded samples. Unknown settings are rejected.
#' @return A list with ordering scores `t`, original-row indices `keep_idx`,
#'   component information, and the realized `k_used`.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(80), 40, 2)
#' fiedler_ordering(X, num_neighbors = 10)
#' @md
#' @export
fiedler_ordering <- function(X, num_neighbors = NULL, control = NULL) {
  opts <- .mpcurve_ordering_options("fiedler", control)
  if (!is.null(num_neighbors)) {
    opts$k <- .mpcurve_integer_setting(num_neighbors, "num_neighbors", 2L)
  }
  do.call(.fiedler_ordering, c(list(X = X), opts))
}

#' t-SNE ordering
#'
#' Use one coordinate of a t-SNE embedding as an initial sample ordering.
#' @inheritParams fiedler_ordering
#' @param component Positive embedding-coordinate index; default 1.
#' @param perplexity Positive t-SNE perplexity, less than `(nrow(X)-1)/3`.
#' @param max_iter Positive integer t-SNE iteration budget; default 500.
#' @param seed Optional random seed.
#' @param control Named advanced settings list: `embedding_dims = 1` (1 to 3,
#'   at least `component`), `scale01 = TRUE`, `orient_by_pc1 = TRUE`.
#'   Scaling maps scores to `[0,1]`; orientation aligns their sign with PC1.
#' @return A list with scores `t` and original-row indices `keep_idx`.
#' @examples
#' \dontrun{
#' tSNE_ordering(X, seed = 1, control = list(embedding_dims = 2))
#' }
#' @md
#' @export
tSNE_ordering <- function(X, component = 1L, perplexity = 10,
                          max_iter = 500L, seed = NULL, control = NULL) {
  opts <- .mpcurve_ordering_options("tSNE", control)
  component <- .mpcurve_integer_setting(component, "component")
  max_iter <- .mpcurve_integer_setting(max_iter, "max_iter")
  .mpcurve_positive_setting(perplexity, "perplexity")
  if (length(perplexity) != 1L) stop("perplexity must be a scalar.", call. = FALSE)
  opts$tSNE_dims <- opts$embedding_dims
  opts$embedding_dims <- NULL
  do.call(.tSNE_ordering, c(list(X = X, component = component,
                               perplexity = perplexity, max_iter = max_iter, seed = seed), opts))
}

#' Principal-curve ordering
#'
#' Order samples by arc length along a fitted principal curve.
#' @inheritParams fiedler_ordering
#' @param smoother One of `"smooth_spline"`, `"lowess"`, `"periodic_lowess"`.
#' @param max_iter Positive integer principal-curve iteration budget; default 10.
#' @param tol Positive principal-curve convergence threshold; default 0.001.
#' @param control Named advanced settings list: `stretch = 2` (nonnegative
#'   endpoint extension), `approx_points = FALSE` (or an integer number of
#'   approximation points), `scale01 = TRUE`, `orient_by_pc1 = TRUE`.
#'   Scaling maps scores to `[0,1]`; orientation aligns their sign with PC1.
#' @return A list with scores `t`, original-row indices `keep_idx`, and `fit`.
#' @examples
#' set.seed(1)
#' pcurve_ordering(matrix(rnorm(80), 40, 2), max_iter = 3)
#' @md
#' @export
pcurve_ordering <- function(X, smoother = c("smooth_spline", "lowess", "periodic_lowess"),
                            max_iter = 10L, tol = 0.001, control = NULL) {
  smoother <- match.arg(smoother)
  max_iter <- .mpcurve_integer_setting(max_iter, "max_iter")
  .mpcurve_positive_setting(tol, "tol")
  if (length(tol) != 1L) stop("tol must be a scalar.", call. = FALSE)
  do.call(.pcurve_ordering, c(list(X = X, smoother = smoother,
    maxit = max_iter, thresh = tol), .mpcurve_ordering_options("pcurve", control)))
}

#' Landmark Isomap ordering
#'
#' Build a nearest-neighbor graph, calculate geodesic distances to landmarks,
#' and use a landmark embedding coordinate to order samples.
#' @inheritParams fiedler_ordering
#' @param num_neighbors Integer nearest-neighbor count; default 15.
#' @param component Positive embedding-coordinate index; default 1.
#' @param seed Optional random seed for landmark selection.
#' @param control Named advanced settings list: `embedding_dims = 1` (at least
#'   `component`), `num_landmarks = NULL` (uses `min(1000, nrow(X))`),
#'   `landmark_method = "random"` (`"kmeans"` also supported), `keep = "giant"`,
#'   `return_full = TRUE`, `orient_by_pc1 = TRUE`, `scale01 = TRUE`.
#'   Giant-component handling retains the largest connected component;
#'   full output has `NA` for excluded samples. Scaling maps scores to `[0,1]`;
#'   orientation aligns their sign with PC1. Inverse-distance stabilization
#'   is an internal numerical rule.
#' @return A list with `t`, `keep_idx`, `n_components`, `embed`, `landmark_idx`,
#'   and `geodesic_to_landmark`.
#' @examples
#' set.seed(1)
#' isomap_ordering(matrix(rnorm(80), 40, 2), num_neighbors = 10, seed = 1)
#' @md
#' @export
isomap_ordering <- function(X, num_neighbors = 15L, component = 1L,
                            seed = NULL, control = NULL) {
  opts <- .mpcurve_ordering_options("isomap", control)
  opts$ndim <- opts$embedding_dims
  opts$landmark <- opts$num_landmarks
  opts$embedding_dims <- opts$num_landmarks <- NULL
  do.call(.isomap_ordering, c(list(X = X,
    k = .mpcurve_integer_setting(num_neighbors, "num_neighbors", 2L),
    component = .mpcurve_integer_setting(component, "component"), seed = seed), opts))
}

.mpcurve_trajectory_options <- function(control, allow_positions = TRUE) {
  defaults <- list(signal_range = c(0.5, 2), linear_slope_range = c(0.8, 2),
    monotone_power_range = c(0.7, 2.5), quadratic_curvature_range = c(0.8, 2),
    quadratic_center_range = c(0.25, 0.75), intercept_sd = 0, sinusoid_freq = 1:4)
  if (allow_positions) defaults["latent_positions"] <- list(NULL)
  .mpcurve_named_settings(control, defaults)
}

#' Simulate multiple latent trajectories
#'
#' Generate sample-by-feature data with one signal-feature block per ordering
#' and an optional noise-only block. Return truth alongside the observations.
#' @param n Number of samples; default 200.
#' @param d_signal Nonnegative signal-feature counts, one per ordering.
#' @param d_noise Number of noise-only features; default 20.
#' @param noise_sd Gaussian observation-noise standard deviation; default 0.3.
#' @param trajectory_family `"sinusoidal"`, `"linear"`, `"monotone"`, or
#'   `"quadratic"`, as one family or one per ordering.
#' @param seed Optional random seed; default 42.
#' @param control Named advanced settings list:
#'   `signal_range = c(0.5, 2)` (signal amplitudes),
#'   `linear_slope_range = c(0.8, 2)`, `monotone_power_range = c(0.7, 2.5)`,
#'   `quadratic_curvature_range = c(0.8, 2)`,
#'   `quadratic_center_range = c(0.25, 0.75)`, `intercept_sd = 0`,
#'   `sinusoid_freq = 1:4` (sampled frequencies), `latent_positions = NULL`.
#'   Non-NULL positions must be an `n`-by-ordering-count matrix. Otherwise
#'   positions are sampled uniformly on `[0,1]`. Parameter ranges affect their
#'   corresponding trajectory families. Unknown settings are rejected.
#' @return A list containing `X`, `true_assign`, `latent_positions`,
#'   `ordering_labels`, `original_order`, and `trajectory_family`.
#' @examples
#' simulate_intrinsic_trajectories(n = 30, d_signal = c(3, 3), d_noise = 2,
#'   trajectory_family = "linear", control = list(linear_slope_range = c(1, 2)))
#' @md
#' @export
simulate_intrinsic_trajectories <- function(n = 200, d_signal = c(30, 30),
    d_noise = 20, noise_sd = 0.3, trajectory_family = "sinusoidal",
    seed = 42L, control = NULL) {
  n <- .mpcurve_integer_setting(n, "n", 2L)
  d_noise <- .mpcurve_integer_setting(d_noise, "d_noise", 0L)
  if (!is.numeric(d_signal) || !length(d_signal) || any(!is.finite(d_signal)) ||
      any(d_signal < 0 | d_signal != floor(d_signal))) {
    stop("d_signal must contain nonnegative integers.", call. = FALSE)
  }
  .mpcurve_scalar_noise_sd(noise_sd)
  do.call(.simulate_intrinsic_trajectories, c(list(n = n, d_signal = d_signal,
    d_noise = d_noise, sigma = noise_sd, trajectory_family = trajectory_family,
    seed = seed), .mpcurve_trajectory_options(control)))
}

#' Simulate two latent trajectories
#'
#' Generate two signal-feature blocks, optionally with related sample orderings.
#' @inheritParams simulate_intrinsic_trajectories
#' @param d1,d2 Number of signal features in each of the two blocks; default 30.
#' @param crossing If `TRUE`, set the second position to a rescaled noisy
#'   reversal of the first, using noise SD 0.3. Otherwise positions are independent.
#' @param control The trajectory-shape settings listed in
#'   [simulate_intrinsic_trajectories()], except `latent_positions`, which this
#'   two-ordering generator constructs from `crossing`.
#' @return A trajectory simulation list with additional `t1` and `t2` vectors.
#' @examples
#' simulate_dual_trajectory(n = 30, d1 = 3, d2 = 3, d_noise = 2, seed = 1)
#' @md
#' @export
simulate_dual_trajectory <- function(n = 200, d1 = 30, d2 = 30, d_noise = 20,
    noise_sd = 0.3, crossing = FALSE,
    trajectory_family = c("sinusoidal", "sinusoidal"), seed = 42L, control = NULL) {
  n <- .mpcurve_integer_setting(n, "n", 2L)
  d1 <- .mpcurve_integer_setting(d1, "d1", 0L)
  d2 <- .mpcurve_integer_setting(d2, "d2", 0L)
  d_noise <- .mpcurve_integer_setting(d_noise, "d_noise", 0L)
  .mpcurve_scalar_noise_sd(noise_sd)
  crossing <- .mpcurve_logical_setting(crossing, "crossing")
  do.call(.simulate_dual_trajectory, c(list(n = n, d1 = d1, d2 = d2,
    d_noise = d_noise, sigma = noise_sd, crossing = crossing,
    trajectory_family = trajectory_family, seed = seed),
    .mpcurve_trajectory_options(control, allow_positions = FALSE)))
}

#' Simulate data from a single-ordering MPCurve model
#'
#' Draw trajectories from random-walk Gaussian priors on a discrete position
#' grid, then sample observations with feature-specific Gaussian noise.
#' @param n Number of samples; default 150.
#' @param d Number of features; default 20.
#' @param num_bins Number of position bins; default 8.
#' @param seed Optional random seed.
#' @param control Named advanced settings list: `rw_order = 2`,
#'   `lambda_range = c(0.5, 3)` (log-uniform precision draws),
#'   `noise_sd_range = c(0.08, 0.18)` (log-uniform noise-SD draws),
#'   `position_weights = NULL` (uniform position probabilities), `ridge = 0.001`.
#'   Supplied nonnegative position weights have one value per bin and positive
#'   total mass. The positive ridge makes the generating trajectory prior proper;
#'   it need not equal the intrinsic prior used by the fitting model.
#' @return A list with observations `X`, assignments, trajectories, true
#'   parameters, position probabilities, and random-walk precision matrices.
#' @examples
#' simulate_mpcurve(n = 30, d = 4, num_bins = 5, seed = 1,
#'   control = list(rw_order = 3))
#' @md
#' @export
simulate_mpcurve <- function(n = 150, d = 20, num_bins = 8, seed = NULL,
                             control = NULL) {
  n <- .mpcurve_integer_setting(n, "n", 2L)
  d <- .mpcurve_integer_setting(d, "d")
  num_bins <- .mpcurve_integer_setting(num_bins, "num_bins", 2L)
  opts <- .mpcurve_named_settings(control, list(rw_order = 2L,
    lambda_range = c(0.5, 3), noise_sd_range = c(0.08, 0.18),
    position_weights = NULL, ridge = 1e-3))
  do.call(.simulate_mpcurve, list(n = n, d = d, K = num_bins, seed = seed,
    rw_q = .mpcurve_integer_setting(opts$rw_order, "rw_order"),
    lambda_range = opts$lambda_range, sigma_range = opts$noise_sd_range,
    pi = opts$position_weights, ridge = opts$ridge))
}

#' Simulate two Gaussian-process feature blocks
#'
#' Draw one independent Matérn process per feature, with one latent sample
#' ordering per block, then optionally add Gaussian measurement noise.
#' @param n Number of samples; default 1000.
#' @param d Even number of features, split equally between the two blocks;
#'   default 16.
#' @param noise_sd Gaussian measurement-noise SD; default 0.05.
#' @param seed Optional random seed.
#' @param control Named advanced settings list: `t_range = c(0, 10)` (uniform
#'   latent-coordinate interval), `range = 5` (Matérn range), `smoothness = 2.5`,
#'   `variance = 3` (process variance), `shift_positive = TRUE`,
#'   `permute_cols = FALSE`, `permute_rows_block2 = TRUE`.
#'   Positive shifting adds a block-specific constant before measurement noise;
#'   permutations are recorded in the returned truth.
#' @return A list containing `X`, `t1`, `t2`, `true_group`, `permut_cols`,
#'   and `row_perm_block2`. When the second block's rows are permuted, its
#'   observation-aligned positions are `t2[row_perm_block2]`.
#' @examples
#' simulate_two_order_gp_dataset(n = 30, d = 4, seed = 1)
#' @md
#' @export
simulate_two_order_gp_dataset <- function(n = 1000, d = 16, noise_sd = 0.05,
                                          seed = NULL, control = NULL) {
  n <- .mpcurve_integer_setting(n, "n", 2L)
  d <- .mpcurve_integer_setting(d, "d", 2L)
  .mpcurve_scalar_noise_sd(noise_sd)
  opts <- .mpcurve_named_settings(control, list(t_range = c(0, 10),
    range = 5, smoothness = 2.5, variance = 3, shift_positive = TRUE,
    permute_cols = FALSE, permute_rows_block2 = TRUE))
  for (nm in c("shift_positive", "permute_cols", "permute_rows_block2")) {
    opts[[nm]] <- .mpcurve_logical_setting(opts[[nm]], nm)
  }
  do.call(.simulate_two_order_gp_dataset, c(list(N = n, D = d,
    noise_sd = noise_sd, seed = seed), opts))
}

.mpcurve_scalar_noise_sd <- function(noise_sd) {
  .mpcurve_positive_setting(noise_sd, "noise_sd", allow_zero = TRUE)
  if (length(noise_sd) != 1L) stop("noise_sd must be a scalar.", call. = FALSE)
  invisible(noise_sd)
}
