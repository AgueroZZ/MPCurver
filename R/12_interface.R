# Validated settings for the ordinary fitting interface.
.mpcurve_integer_setting <- function(value, name, minimum = 1L) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value != floor(value) || value < minimum || value > .Machine$integer.max) {
    stop(name, " must be a single integer >= ", minimum, ".", call. = FALSE)
  }
  as.integer(value)
}

.mpcurve_positive_setting <- function(value, name, allow_zero = FALSE) {
  if (!is.numeric(value) || !length(value) || any(!is.finite(value)) ||
      any(if (allow_zero) value < 0 else value <= 0)) {
    stop(name, " must contain finite ",
         if (allow_zero) "nonnegative" else "positive", " values.", call. = FALSE)
  }
  value
}

.mpcurve_bounds_setting <- function(value, name) {
  .mpcurve_positive_setting(value, name)
  if (length(value) != 2L || value[1] > value[2]) {
    stop(name, " must contain ordered lower and upper bounds.", call. = FALSE)
  }
  as.numeric(value)
}

.mpcurve_logical_setting <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  value
}

.mpcurve_interface_options <- function(value, constructor, name) {
  if (is.null(value)) return(constructor())
  if (!is.list(value) || (length(value) &&
      (is.null(names(value)) || any(!nzchar(names(value))) || anyDuplicated(names(value))))) {
    stop(name, " must be a list with unique named settings.", call. = FALSE)
  }
  unknown <- setdiff(names(value), names(formals(constructor)))
  if (length(unknown)) {
    stop("Unknown ", name, " option(s): ", paste(unknown, collapse = ", "),
         ". See ", deparse(substitute(constructor)), "().", call. = FALSE)
  }
  do.call(constructor, value)
}

#' Advanced model and fitting settings
#'
#' Configure trajectory priors, initial parameters, numerical bounds, and
#' feature-assignment annealing for [fit_mpcurve()]. Omitted settings retain
#' their defaults. The returned settings are checked again at fitting time.
#'
#' @param rw_order Positive integer random-walk order. Default: `2` (RW2).
#'   RW3 uses `3`. The realized number of position bins must exceed this order.
#' @param lambda_init Positive initial smoothness precision. Default: `1`.
#'   A larger value penalizes trajectory differences more strongly for the same
#'   random-walk order. Accepts a scalar or one value per feature. Partition
#'   fits also accept a feature-by-ordering matrix or its column-major vector;
#'   a feature vector is shared across orderings. Initial values are clipped
#'   to `lambda_bounds`, including when `fix_lambda = TRUE`.
#' @param fix_lambda Keep the initialized precision fixed throughout fitting?
#'   Default: `FALSE`, allowing empirical-Bayes updates. For precision fixed
#'   at 5, supply both `lambda_init = 5` and `fix_lambda = TRUE`.
#' @param sigma2_init Optional positive initial noise variances, as a scalar
#'   or one value per feature. Default: `NULL`, using initialization estimates.
#'   Noise variances are subsequently updated and shared across orderings for
#'   each feature. These are variances, not standard deviations. Cannot be
#'   supplied with known measurement standard deviations `S`.
#' @param lambda_bounds Positive `c(lower, upper)` precision limits.
#'   Default: `c(1e-10, 1e10)`. Constrain initialized and updated precisions;
#'   bounds do not fix a precision to its initial value.
#' @param sigma2_bounds Positive `c(lower, upper)` noise-variance limits.
#'   Default: `c(1e-10, 1e10)`. Apply to estimated variances, in squared data
#'   units. They do not constrain supplied measurement standard deviations `S`.
#' @param ridge Nonnegative diagonal addition to the base random-walk precision.
#'   Default: `0`, using an intrinsic prior. A positive value gives a proper
#'   prior and penalizes directions unpenalized by the intrinsic random walk;
#'   it changes the model, not just the numerical solver.
#' @param lambda_sd_prior_rate Optional positive exponential-prior rate on
#'   `1 / sqrt(lambda)`. Default: `NULL`; zero also disables the extra prior.
#'   Larger rates favor smaller smoothness scales and hence larger precision.
#'   This adds regularization to precision estimation. With `fix_lambda = TRUE`
#'   it does not update the fixed precision.
#' @param position_prior_weights Optional nonnegative vector with one weight
#'   per realized bin and positive total mass. Default: `NULL`. Weights are
#'   normalized to sum to one and shared across orderings. With
#'   `position_prior = "adaptive"` they set initial position probabilities;
#'   with `"fixed"` they set probabilities held fixed. Omitted fixed weights
#'   give a uniform distribution; omitted adaptive weights use initialization.
#' @param partition_prior_weights Optional nonnegative vector with one weight
#'   per ordering and positive total mass, normalized to sum to one.
#'   Default: `NULL`, giving uniform weights under `partition_prior = "fixed"`.
#'   Custom values require `partition_prior = "fixed"` and are held fixed.
#'   They specify ordering probabilities for features, not position-bin weights.
#' @param anneal_steps Nonnegative integer number of partition temperature steps.
#'   Default: `0`, skipping annealing and fitting directly at temperature 1.
#'   A positive value enables a geometric schedule from `anneal_start` to 1.
#'   With one step, that step uses temperature 1. Applies only to multi-ordering
#'   inference.
#' @param anneal_start Positive starting partition temperature. Default: `5`.
#'   Higher temperatures encourage more diffuse feature assignments early in
#'   fitting. Used only when `anneal_steps > 0`. The schedule endpoint and final
#'   convergence temperature are 1.
#' @param anneal_sweeps Positive integer joint coordinate sweeps per temperature.
#'   Default: `1`. Partition annealing performs `anneal_steps * anneal_sweeps`
#'   sweeps, in addition to initialization and the final `max_iter` budget.
#' @param convergence Stopping-rule scale. Default: `"normalized"`, checking
#'   `abs(delta_ELBO) / (nrow(X) * ncol(X))` against the main `tol` argument.
#'   `"relative"` retains historical scaling by the previous ELBO magnitude
#'   plus the backend's numerical offset. Stored objective values are unchanged.
#' @param effective_count_tol Numerical-zero tolerance for estimated ordering
#'   prior weights. Default: `1e-8`. With `intrinsic_dim = "auto"` and an
#'   adaptive partition prior, an ordering is removed when its estimated prior
#'   weight is at or below `effective_count_tol / P`, where `P = ncol(X)`.
#'   Empirical Bayes estimates this weight as `sum_j Pr(Z_j = m) / P`.
#'   The tolerance must be nonnegative and less than `P / intrinsic_dim` for
#'   the initial fit. Automatic adaptive pruning warns when the tolerance
#'   exceeds `1e-6`, since large values may remove supported orderings.
#'   This warning boundary applies before division by `P`. Pruning stops with
#'   an error if any feature has no positive remaining assignment probability.
#'   Remaining probabilities are normalized
#'   and the final objective is recomputed without further CAVI updates.
#'   An explicitly specified dimension or a fixed partition prior retains all
#'   fitted orderings.
#' @section Smoothness and noise:
#' `rw_order` specifies which trajectory differences are penalized;
#' `lambda_init` and `fix_lambda` determine whether precision is estimated or
#' held fixed. `sigma2_init` supplies a starting noise variance. The two bounds
#' constrain parameter values. `ridge` and `lambda_sd_prior_rate` change prior
#' regularization. Known observation errors are supplied through `S` in
#' [fit_mpcurve()].
#' @section Position and partition probabilities:
#' `position_prior_weights` describes sample positions along each ordering.
#' `partition_prior_weights` describes which ordering a feature follows.
#' Their update modes are selected by the main `position_prior` and
#' `partition_prior` arguments.
#' @section Optimization and reporting:
#' Fitting uses temperature 1 throughout by default. Optional annealing applies
#' only when the realized fit has multiple orderings and `anneal_steps > 0`.
#' It precedes the temperature-1 phase governed by the main `max_iter` and `tol`
#' arguments. `convergence` selects tolerance scaling;
#' `effective_count_tol` controls removal of empty orderings after automatic
#' dimension estimation.
#' @return A validated settings list for the `control` argument. Supply only
#'   settings to change; omitted settings retain their defaults. Inspect the
#'   complete default list with `mpcurve_control()`. Unknown options are rejected
#'   by [fit_mpcurve()]. Initialization settings belong in
#'   [mpcurve_init_control()].
#' @examples
#' # Change the smoothness prior to RW3; retain precision updates.
#' mpcurve_control(rw_order = 3)
#'
#' # Hold precision at 5 throughout fitting.
#' mpcurve_control(lambda_init = 5, fix_lambda = TRUE)
#'
#' # Constrain estimated precision and noise variance.
#' mpcurve_control(lambda_bounds = c(0.1, 100), sigma2_bounds = c(0.01, 10))
#'
#' # Budget 10 annealing sweeps before the main max_iter fitting phase.
#' mpcurve_control(anneal_steps = 5, anneal_start = 3, anneal_sweeps = 2)
#' @md
#' @export
mpcurve_control <- function(
    rw_order = 2L, lambda_init = 1, fix_lambda = FALSE, sigma2_init = NULL,
    lambda_bounds = c(1e-10, 1e10), sigma2_bounds = c(1e-10, 1e10),
    ridge = 0, lambda_sd_prior_rate = NULL,
    position_prior_weights = NULL, partition_prior_weights = NULL,
    anneal_steps = 0L, anneal_start = 5, anneal_sweeps = 1L,
    convergence = c("normalized", "relative"), effective_count_tol = 1e-8) {
  rw_order <- .mpcurve_integer_setting(rw_order, "rw_order")
  lambda_init <- .mpcurve_positive_setting(lambda_init, "lambda_init")
  fix_lambda <- .mpcurve_logical_setting(fix_lambda, "fix_lambda")
  if (!is.null(sigma2_init)) .mpcurve_positive_setting(sigma2_init, "sigma2_init")
  lambda_bounds <- .mpcurve_bounds_setting(lambda_bounds, "lambda_bounds")
  sigma2_bounds <- .mpcurve_bounds_setting(sigma2_bounds, "sigma2_bounds")
  .mpcurve_positive_setting(ridge, "ridge", allow_zero = TRUE)
  if (length(ridge) != 1L) stop("ridge must be a scalar.", call. = FALSE)
  lambda_sd_prior_rate <- .normalize_lambda_sd_prior_rate(lambda_sd_prior_rate)
  for (weights in list(position_prior_weights, partition_prior_weights)) {
    if (!is.null(weights)) {
      .mpcurve_positive_setting(weights, "prior weights", allow_zero = TRUE)
      if (sum(weights) <= 0) stop("Prior weights must have positive mass.", call. = FALSE)
    }
  }
  anneal_steps <- .mpcurve_integer_setting(anneal_steps, "anneal_steps", minimum = 0L)
  .mpcurve_positive_setting(anneal_start, "anneal_start")
  if (length(anneal_start) != 1L) stop("anneal_start must be a scalar.", call. = FALSE)
  anneal_sweeps <- .mpcurve_integer_setting(anneal_sweeps, "anneal_sweeps")
  convergence <- match.arg(convergence)
  .mpcurve_positive_setting(effective_count_tol, "effective_count_tol", allow_zero = TRUE)
  if (length(effective_count_tol) != 1L) stop("effective_count_tol must be a scalar.", call. = FALSE)
  list(rw_order = rw_order, lambda_init = lambda_init, fix_lambda = fix_lambda,
       sigma2_init = sigma2_init, lambda_bounds = lambda_bounds,
       sigma2_bounds = sigma2_bounds, ridge = ridge,
       lambda_sd_prior_rate = lambda_sd_prior_rate,
       position_prior_weights = position_prior_weights,
       partition_prior_weights = partition_prior_weights,
       anneal_steps = anneal_steps, anneal_start = anneal_start,
       anneal_sweeps = anneal_sweeps, convergence = convergence,
       effective_count_tol = effective_count_tol)
}

#' Advanced initialization settings
#'
#' Configure feature grouping, automatic ordering-count initialization, and
#' discretization for [fit_mpcurve()]. Each feature block uses the single
#' method selected by `initial_method`. Multiple candidate methods must be
#' fitted with separate calls.
#'
#' @param discretization Initial position discretization: `"quantile"`,
#'   `"equal"`, or `"kmeans"`.
#' @param similarity_metric Feature similarity: `"spearman"`, `"pearson"`,
#'   `"spline_r2"`, or `"smooth_fit"`. `NULL` retains Spearman for specified
#'   ordering counts and spline R-squared for automatic initialization.
#'   See the Feature similarity section for definitions.
#' @param max_intrinsic_dim Upper bound for automatic initialization.
#' @param min_cluster_size Minimum features per eligible automatic cluster.
#' @param spline_r2_df Natural-spline degrees of freedom for spline similarity.
#' @param cluster_linkage Hierarchical-clustering linkage on similarity distance;
#'   defaults to `"single"`.
#' @param smooth_fit_lambda_mode,smooth_fit_lambda_value Precision settings used
#'   only by smooth-fit similarity.
#' @param pca_components Optional PCA-component indices within feature groups,
#'   shared as a scalar or supplied per group. The default is each group's PC1.
#' @param responsibilities_init Optional single-ordering responsibility matrix.
#' @param fits_init Optional list of initial CAVI fits for partition inference.
#' @param method_args Named arguments forwarded to the selected ordering helper.
#'   For `initial_method = "fiedler"` or `"isomap"`, set the graph's neighbor
#'   count with `method_args = list(num_neighbors = 20)`. Available arguments
#'   are documented in [fiedler_ordering()], [isomap_ordering()], and the other
#'   ordering helpers.
#'   Isomap defaults to the smallest count connecting all samples, computed
#'   separately within each feature group. The realized count is recorded as
#'   `k_used` in the fit's initialization metadata.
#' @section Initializing multiple orderings:
#' Feature similarity scores are converted to distances as one minus similarity.
#' Hierarchical clustering uses
#' `cluster_linkage` (single linkage by default), and the tree is cut into the
#' specified number of feature groups. The selected `initial_method` is applied
#' independently within each group to initialize its sample ordering. With
#' `initial_method = "PCA"`, each group uses its own first principal component.
#' These groups initialize the model; feature assignments are updated during
#' fitting.
#'
#' Compare initialization methods through separate fits with common model and
#' fitting settings, using final temperature-one objectives and fitted results.
#' @section Feature similarity:
#' The available metrics are:
#' * `"spearman"`: absolute Spearman rank correlation. This is the default when
#'   `intrinsic_dim` is specified as an integer.
#' * `"pearson"`: absolute Pearson correlation. Taking the absolute value allows
#'   features with reversed profiles to share an ordering.
#' * `"spline_r2"`: order samples by one feature and regress the other on a
#'   natural cubic spline of the resulting sample ranks. Similarity is the
#'   larger R-squared from the two directions, clipped between 0 and 1. The spline
#'   degrees of freedom are set by `spline_r2_df` (default 5). Nonconstant
#'   features are centered and scaled to unit sample variance for this
#'   calculation only; fitting retains the original data and measurement
#'   errors. This is the default for `intrinsic_dim = "auto"`.
#' * `"smooth_fit"`: order samples by one feature and compare a random-walk
#'   smooth fit of the other with a constant-mean model. Take the larger
#'   positive log-evidence improvement from the two directions, then divide
#'   by the largest improvement across feature pairs to obtain scores between
#'   0 and 1. The smoothness precision is optimized or fixed according to
#'   `smooth_fit_lambda_mode` and `smooth_fit_lambda_value`. With the default
#'   intrinsic random-walk prior (`ridge = 0`), this uses pseudo-evidence;
#'   `ridge > 0` gives a proper Gaussian prior for marginal-evidence comparison.
#'
#' Constant features have zero similarity to other features. Self-similarity
#' is one. Metric calculations affect initialization, while the model is fitted
#' to the original data.
#' @param on_failure Action when the requested ordering initializer fails after
#'   argument validation. Default: `"pca"`, warning with the original error and
#'   retrying with PCA component 1. `"error"` stops immediately. Invalid settings
#'   always error. If PCA itself fails, fitting stops. The realized method and
#'   fallback reason are recorded in the fit's initialization metadata.
#' @return A validated settings list for the `init_control` argument.
#' @examples
#' mpcurve_init_control(max_intrinsic_dim = 6, min_cluster_size = 3)
#' mpcurve_init_control(method_args = list(num_neighbors = 10))
#' mpcurve_init_control(on_failure = "error")
#' @md
#' @export
mpcurve_init_control <- function(
    discretization = c("quantile", "equal", "kmeans"), similarity_metric = NULL,
    max_intrinsic_dim = 8L, min_cluster_size = 2L, spline_r2_df = 5L,
    cluster_linkage = "single", smooth_fit_lambda_mode = c("optimize", "fixed"),
    smooth_fit_lambda_value = 1, pca_components = NULL,
    responsibilities_init = NULL, fits_init = NULL, method_args = list(),
    on_failure = c("pca", "error")) {
  on_failure <- match.arg(on_failure)
  if (!is.null(pca_components)) {
    if (!is.numeric(pca_components) || !length(pca_components)) {
      stop("pca_components must contain positive integers.", call. = FALSE)
    }
    pca_components <- vapply(pca_components, .mpcurve_integer_setting,
      integer(1), name = "pca_components")
  }
  discretization <- match.arg(discretization)
  if (!is.null(similarity_metric)) {
    similarity_metric <- match.arg(similarity_metric,
      c("spearman", "pearson", "spline_r2", "smooth_fit"))
  }
  max_intrinsic_dim <- .mpcurve_integer_setting(max_intrinsic_dim, "max_intrinsic_dim")
  min_cluster_size <- .mpcurve_integer_setting(min_cluster_size, "min_cluster_size")
  spline_r2_df <- .mpcurve_integer_setting(spline_r2_df, "spline_r2_df")
  cluster_linkage <- .cavi_validate_cluster_linkage(cluster_linkage)
  smooth_fit_lambda_mode <- match.arg(smooth_fit_lambda_mode)
  .mpcurve_positive_setting(smooth_fit_lambda_value, "smooth_fit_lambda_value")
  if (length(smooth_fit_lambda_value) != 1L) stop("smooth_fit_lambda_value must be a scalar.", call. = FALSE)
  if (!is.list(method_args) || (length(method_args) &&
      (is.null(names(method_args)) || any(!nzchar(names(method_args))) || anyDuplicated(names(method_args))))) {
    stop("method_args must contain uniquely named ordering-helper arguments.", call. = FALSE)
  }
  if (any(names(method_args) %in% c("X", "method"))) {
    stop("method_args cannot override X or method.", call. = FALSE)
  }
  list(discretization = discretization, on_failure = on_failure,
       similarity_metric = similarity_metric, max_intrinsic_dim = max_intrinsic_dim,
       min_cluster_size = min_cluster_size, spline_r2_df = spline_r2_df,
       cluster_linkage = cluster_linkage, smooth_fit_lambda_mode = smooth_fit_lambda_mode,
       smooth_fit_lambda_value = smooth_fit_lambda_value, pca_components = pca_components,
       responsibilities_init = responsibilities_init, fits_init = fits_init,
       method_args = method_args)
}

#' Estimate sample orderings and smooth feature trajectories
#'
#' Fit the MPCurve model by coordinate-ascent variational inference. One
#' ordering describes all features when `intrinsic_dim = 1`. With several
#' orderings, structural inference jointly estimates sample positions,
#' conditional trajectories, and soft feature assignments.
#'
#' @param X Numeric sample-by-feature matrix or all-numeric data frame, with
#'   finite values only. Rows are samples and columns are features. Handle
#'   missing values before fitting; the function does not impute or standardize.
#'   Sample and feature names are retained in public results.
#' @param S Optional known measurement standard deviations, as a length-d feature
#'   vector or matrix matching `X`. Values must be finite and nonnegative, in
#'   the same units and positional order as `X`; names do not trigger reordering.
#'   Zero SDs are floored at `sqrt(.Machine$double.eps)` for computation. Otherwise
#'   feature noise variances are estimated.
#' @param num_bins Requested number of latent-position bins. `NULL` uses
#'   `max(control$rw_order + 1, min(50, floor(nrow(X) / 5)))`.
#'   The sample count must be at least `control$rw_order + 1`. Collapsed initialization cuts can
#'   reduce the realized count in a single-ordering fit.
#' @param intrinsic_dim Positive model intrinsic dimension (number of orderings),
#'   or `"auto"` to estimate it. Automatic fitting initializes the dimension
#'   from an eligible feature-similarity cut. With an adaptive partition prior,
#'   orderings with estimated prior weight at or below
#'   `control$effective_count_tol / ncol(X)` are removed without refitting.
#'   An integer retains the specified dimension.
#' @param initial_method One initialization method: `"PCA"`, `"fiedler"`,
#'   `"pcurve"`, `"tSNE"`, `"random"`, or `"isomap"`. Under similarity
#'   grouping this method is applied independently within each feature
#'   block. The default PCA uses each block's own PC1. Method vectors are rejected.
#' @param position_prior Estimate position probabilities (`"adaptive"`) or keep
#'   them fixed (`"fixed"`, uniform unless custom control weights are supplied).
#' @param partition_prior Estimate ordering probabilities by empirical Bayes
#'   (`"adaptive"`) or keep them fixed (`"fixed"`, uniform by default).
#' @param max_iter Maximum CAVI sweeps for one ordering, or maximum joint sweeps
#'   at temperature 1 for multiple orderings. Annealing is disabled by default;
#'   enabling it adds `anneal_steps * anneal_sweeps` sweeps before this phase,
#'   configured through [mpcurve_control()]. Zero skips this fitting phase.
#' @param tol Shared ELBO-change stopping tolerance. The default rule uses
#'   `abs(delta_ELBO) / (nrow(X) * ncol(X))`; partition stopping is checked
#'   at temperature 1, after any optional annealing. Existing numerical decrease
#'   checks apply.
#'   Zero disables early stopping. A small increment does not bound the
#'   remaining optimization gap.
#' @param verbose Print fitting progress?
#' @param init_control Settings from [mpcurve_init_control()] or a named list
#'   of its options. `NULL` uses defaults. Arguments for the selected ordering
#'   helper go in `method_args`; for example, `num_neighbors` controls the
#'   nearest-neighbor graphs in [fiedler_ordering()] and [isomap_ordering()].
#' @param control Settings from [mpcurve_control()] or a named list of its
#'   options. `NULL` uses defaults. Unknown settings are rejected.
#' @details
#' Smoothness uses RW2 by default, with feature-specific empirical-Bayes
#' precision updates. Use `control = mpcurve_control(rw_order = 3)` for RW3,
#' or `control = mpcurve_control(lambda_init = 5, fix_lambda = TRUE)` to fix
#' precision at 5. Setting a positive ridge changes the smoothness prior.
#'
#' The partition variational family is
#' \deqn{q(C)\prod_j q(Z_j)q(U_j\mid Z_j).}
#' Noise variance is shared across candidate orderings for each feature.
#' The result's `intrinsic_dim` is the actual number of orderings in the
#' returned model. The compatibility field `model_intrinsic_dim` has the same
#' value for current fits. The grid size is stored in `K`; see [mpcurve].
#' @return One `mpcurve` fit. Use [summary()], [plot()], [fitted_prior()], and
#'   [fitted_positions()], [fitted_trajectories()], [fitted_assignments()], and
#'   [do_mpcurve()] to inspect, extract results, or continue it. Automatic initialization also
#'   stores candidate silhouettes and selected clusters in
#'   `$dimension_initialization`. `$dimension_estimation` records the initial
#'   and final dimensions and any removed orderings. If removal occurs,
#'   preceding fitting traces are archived there, and the returned model's
#'   objective trace starts at its recomputed final objective.
#' @examples
#' sim <- simulate_mpcurve(n = 40, d = 4, num_bins = 5, seed = 1)
#' fit <- fit_mpcurve(sim$X, num_bins = 5, max_iter = 2)
#' fit_rw3 <- fit_mpcurve(sim$X, num_bins = 5, max_iter = 2,
#'                       control = mpcurve_control(rw_order = 3))
#' fit_fixed <- fit_mpcurve(sim$X, num_bins = 5, max_iter = 2,
#'                         control = mpcurve_control(lambda_init = 5,
#'                                                   fix_lambda = TRUE))
#' fit_fiedler <- fit_mpcurve(sim$X, initial_method = "fiedler",
#'   num_bins = 5, max_iter = 2,
#'   init_control = mpcurve_init_control(
#'     method_args = list(num_neighbors = 20)))
#' fit_isomap <- fit_mpcurve(sim$X, initial_method = "isomap",
#'   num_bins = 5, max_iter = 2,
#'   init_control = mpcurve_init_control(
#'     method_args = list(num_neighbors = 20)))
#' @md
#' @export
fit_mpcurve <- function(
    X, S = NULL, num_bins = NULL, intrinsic_dim = 1L, initial_method = "PCA",
    position_prior = c("adaptive", "fixed"),
    partition_prior = c("adaptive", "fixed"), max_iter = 100L, tol = 1e-6,
    verbose = FALSE, init_control = NULL, control = NULL) {
  X <- .mpcurve_validate_data(X, S)
  if (!is.character(initial_method) || length(initial_method) != 1L || is.na(initial_method)) {
    stop("initial_method must select exactly one initialization method.", call. = FALSE)
  }
  initial_method <- match.arg(initial_method,
    c("PCA", "fiedler", "pcurve", "tSNE", "random", "isomap"))
  position_prior <- match.arg(position_prior)
  partition_prior <- match.arg(partition_prior)
  max_iter <- .mpcurve_integer_setting(max_iter, "max_iter", minimum = 0L)
  .mpcurve_positive_setting(tol, "tol", allow_zero = TRUE)
  if (length(tol) != 1L) stop("tol must be a scalar.", call. = FALSE)
  verbose <- .mpcurve_logical_setting(verbose, "verbose")
  if (!is.null(num_bins)) num_bins <- .mpcurve_integer_setting(num_bins, "num_bins", 2L)
  automatic <- identical(intrinsic_dim, "auto")
  if (!automatic) intrinsic_dim <- .mpcurve_integer_setting(intrinsic_dim, "intrinsic_dim")
  init <- .mpcurve_interface_options(init_control, mpcurve_init_control, "init_control")
  ctrl <- .mpcurve_interface_options(control, mpcurve_control, "control")
  helper <- switch(initial_method, PCA = PCA_ordering, fiedler = fiedler_ordering,
                   pcurve = pcurve_ordering, tSNE = tSNE_ordering,
                   isomap = isomap_ordering, random = function() NULL)
  unknown <- setdiff(names(init$method_args), setdiff(names(formals(helper)), "X"))
  if (length(unknown)) {
    stop("Unknown initialization argument(s) for ", initial_method, ": ",
         paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  if (initial_method %in% c("fiedler", "tSNE", "pcurve", "isomap")) {
    .mpcurve_ordering_options(initial_method, init$method_args$control)
  }
  if (initial_method == "PCA" && "component" %in% names(init$method_args)) {
    stop("Use pca_components to choose PCA components.", call. = FALSE)
  }
  .mpcurve_validate_ordering_args(as.matrix(X), initial_method, init$method_args)
  if (!is.null(init$pca_components) && initial_method != "PCA") {
    stop("pca_components requires initial_method = 'PCA'.", call. = FALSE)
  }
  if (!is.null(ctrl$partition_prior_weights) && partition_prior != "fixed") {
    stop("partition_prior_weights requires partition_prior = 'fixed'.", call. = FALSE)
  }
  if (!is.null(S) && !is.null(ctrl$sigma2_init)) {
    stop("sigma2_init cannot be supplied when S is non-NULL.", call. = FALSE)
  }
  if (!is.null(init$responsibilities_init) && (automatic || intrinsic_dim != 1L)) {
    stop("responsibilities_init requires intrinsic_dim = 1.", call. = FALSE)
  }
  if (!is.null(init$fits_init) && (automatic || intrinsic_dim == 1L)) {
    stop("fits_init requires a specified intrinsic_dim >= 2.", call. = FALSE)
  }
  args <- list(
    X = X, S = S, K = num_bins, intrinsic_dim = intrinsic_dim, method = initial_method,
    init_on_failure = init$on_failure,
    iter = max_iter, max_converge_iter = max_iter, tol = tol, tol_outer = tol,
    position_prior = position_prior, partition_prior = partition_prior,
    rw_q = ctrl$rw_order, lambda = ctrl$lambda_init, fix_lambda = ctrl$fix_lambda,
    sigma2_init = ctrl$sigma2_init, ridge = ctrl$ridge,
    lambda_min = ctrl$lambda_bounds[1], lambda_max = ctrl$lambda_bounds[2],
    sigma_min = ctrl$sigma2_bounds[1], sigma_max = ctrl$sigma2_bounds[2],
    lambda_sd_prior_rate = ctrl$lambda_sd_prior_rate,
    position_prior_init = ctrl$position_prior_weights,
    partition_prior_init = ctrl$partition_prior_weights,
    n_outer = ctrl$anneal_steps, inner_iter = ctrl$anneal_sweeps,
    T_start = ctrl$anneal_start, T_end = 1, convergence = ctrl$convergence,
    effective_count_tol = ctrl$effective_count_tol, verbose = verbose,
    discretization = init$discretization,
    max_intrinsic_dim = init$max_intrinsic_dim,
    similarity_min_cluster_size = init$min_cluster_size, spline_r2_df = init$spline_r2_df,
    cluster_linkage = init$cluster_linkage,
    smooth_fit_lambda_mode = init$smooth_fit_lambda_mode,
    smooth_fit_lambda_value = init$smooth_fit_lambda_value
  )
  if (!is.null(init$similarity_metric)) args$similarity_metric <- init$similarity_metric
  if (!is.null(init$pca_components)) args$pca_components <- init$pca_components
  if (!is.null(init$responsibilities_init)) args$responsibilities_init <- init$responsibilities_init
  if (!is.null(init$fits_init)) args$fits_init <- init$fits_init
  if (automatic || intrinsic_dim >= 2L) {
    args$method_args <- init$method_args
  } else {
    if (!is.null(init$pca_components)) {
      args$pca_components <- NULL
      init$method_args$component <- .mpcurve_integer_setting(init$pca_components, "pca_components")
    }
    args$method_args <- init$method_args
  }
  do.call(.fit_mpcurve, args)
}

# Stable names also support selector tests that replace the fitting function.
.mpcurve_fit_argument_names <- names(formals(fit_mpcurve))

#' Advanced continuation settings
#'
#' Override selected settings when resuming with [do_mpcurve()]. Every default
#' is `NULL`, which inherits the value in the fitted model. This differs from
#' [mpcurve_control()], whose defaults configure a new fit.
#'
#' @param lambda_init Optional positive precision reset. Accepts a scalar or
#'   one value per feature; partition fits also accept a feature-by-ordering
#'   matrix or column-major vector. Values are clipped to `lambda_bounds`.
#'   Resetting precision does not change whether it is updated.
#' @param fix_lambda Optional logical update mode. `TRUE` fixes precision at
#'   its current or reset value; `FALSE` enables empirical-Bayes updates.
#' @param lambda_bounds Optional positive `c(lower, upper)` precision bounds.
#' @param sigma2_bounds Optional positive `c(lower, upper)` estimated noise
#'   variance bounds. Supplied measurement errors remain unchanged.
#' @param lambda_sd_prior_rate Optional positive exponential-prior rate on
#'   `1 / sqrt(lambda)`. Zero removes the penalty; `NULL` inherits it.
#' @param convergence Optional `"normalized"` or `"relative"` stopping scale.
#'   The tolerance is supplied through `tol` in [do_mpcurve()].
#' @return A validated list of continuation overrides.
#' @examples
#' mpcurve_continue_control(lambda_init = 5, fix_lambda = TRUE)
#' mpcurve_continue_control(fix_lambda = FALSE, lambda_bounds = c(0.1, 100))
#' @md
#' @export
mpcurve_continue_control <- function(lambda_init = NULL, fix_lambda = NULL,
    lambda_bounds = NULL, sigma2_bounds = NULL,
    lambda_sd_prior_rate = NULL, convergence = NULL) {
  if (!is.null(lambda_init)) .mpcurve_positive_setting(lambda_init, "lambda_init")
  if (!is.null(fix_lambda)) .mpcurve_logical_setting(fix_lambda, "fix_lambda")
  if (!is.null(lambda_bounds)) lambda_bounds <- .mpcurve_bounds_setting(lambda_bounds, "lambda_bounds")
  if (!is.null(sigma2_bounds)) sigma2_bounds <- .mpcurve_bounds_setting(sigma2_bounds, "sigma2_bounds")
  if (!is.null(lambda_sd_prior_rate)) {
    .mpcurve_positive_setting(lambda_sd_prior_rate, "lambda_sd_prior_rate", allow_zero = TRUE)
    if (length(lambda_sd_prior_rate) != 1L) stop("lambda_sd_prior_rate must be a scalar.", call. = FALSE)
  }
  if (!is.null(convergence)) convergence <- match.arg(convergence, c("normalized", "relative"))
  list(lambda_init = lambda_init, fix_lambda = fix_lambda,
       lambda_bounds = lambda_bounds, sigma2_bounds = sigma2_bounds,
       lambda_sd_prior_rate = lambda_sd_prior_rate, convergence = convergence)
}
