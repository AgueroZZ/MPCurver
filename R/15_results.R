# Public data and result contracts, independent of the number of orderings.
.mpcurve_validate_data <- function(X, S = NULL) {
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a numeric sample-by-feature matrix or data frame.", call. = FALSE)
  }
  if (is.data.frame(X)) {
    numeric_columns <- vapply(X, function(x) is.numeric(x) && !is.complex(x), logical(1))
    if (any(!numeric_columns)) {
      j <- which(!numeric_columns)[1L]
      stop(sprintf("X must be numeric; feature %d ('%s') is not numeric.",
                   j, names(X)[j]), call. = FALSE)
    }
    X <- as.matrix(X)
  }
  if (!is.numeric(X) || is.complex(X)) {
    stop("X must be numeric and real-valued; character, factor, logical, and complex data are not supported.",
         call. = FALSE)
  }
  if (nrow(X) < 2L || ncol(X) < 1L) {
    stop("X must have at least two samples and one feature.", call. = FALSE)
  }
  bad <- which(!is.finite(X), arr.ind = TRUE)
  if (nrow(bad)) {
    i <- bad[1L, 1L]
    j <- bad[1L, 2L]
    label <- function(index, labels) {
      if (is.null(labels)) as.character(index) else sprintf("%d ('%s')", index, labels[index])
    }
    stop(sprintf(paste0("X must contain only finite values; sample %s, feature %s ",
                        "contains NA, NaN, or Inf. Handle missing values before fitting."),
                 label(i, rownames(X)), label(j, colnames(X))), call. = FALSE)
  }
  if (!is.null(S)) {
    if (!is.numeric(S) || is.complex(S) || any(!is.finite(S)) || any(S < 0)) {
      stop("S must contain only finite nonnegative numeric standard deviations.", call. = FALSE)
    }
    valid_shape <- if (is.matrix(S)) identical(dim(S), dim(X)) else
      is.null(dim(S)) && length(S) == ncol(X)
    if (!valid_shape) {
      stop("S must be NULL, a length-d feature vector, or an n x d matrix matching X.",
           call. = FALSE)
    }
  }
  X
}

.mpcurve_feature_label <- function(x, j) {
  labels <- colnames(x$data)
  if (is.null(labels)) sprintf("dim %d", j) else labels[j]
}

.mpcurve_plot_dimensions <- function(dims, d, feature_names = NULL) {
  if (is.null(dims)) return(seq_len(min(2L, d)))
  if (length(dims) < 1L || length(dims) > 2L || anyNA(dims)) {
    stop("dims must contain one or two feature indices or names, without missing values.",
         call. = FALSE)
  }
  if (is.character(dims)) {
    if (is.null(feature_names)) {
      stop("Named dims require feature names in the fitted data.", call. = FALSE)
    }
    ambiguous <- dims[vapply(dims, function(nm) sum(feature_names == nm, na.rm = TRUE) > 1L,
                            logical(1))]
    if (length(ambiguous)) {
      stop("Ambiguous feature name in dims: ", paste(unique(ambiguous), collapse = ", "),
           ". Use column indices.", call. = FALSE)
    }
    indices <- match(dims, feature_names)
    if (anyNA(indices)) {
      stop("Unknown feature name in dims: ", paste(dims[is.na(indices)], collapse = ", "),
           ".", call. = FALSE)
    }
    return(indices)
  }
  if (!is.numeric(dims) || any(!is.finite(dims)) ||
      any(dims != floor(dims)) || any(dims < 1L) || any(dims > d)) {
    stop("dims out of range or invalid: use integer feature indices between 1 and ",
         d, ".", call. = FALSE)
  }
  as.integer(dims)
}

.mpcurve_label_results <- function(x) {
  samples <- rownames(x$data)
  features <- colnames(x$data)
  label_gamma <- function(gamma) {
    if (!is.null(gamma)) rownames(gamma) <- samples
    gamma
  }
  label_mu <- function(mu) {
    if (!is.null(mu)) rownames(mu) <- features
    mu
  }
  label_locations <- function(locations) {
    if (is.null(locations)) return(NULL)
    for (type in c("mean", "map")) {
      for (quantity in c("index", "pseudotime")) {
        names(locations[[type]][[quantity]]) <- samples
      }
    }
    locations
  }
  if (.mpcurve_is_partition(x)) {
    x$gamma <- lapply(x$gamma, label_gamma)
    x$params$mu <- lapply(x$params$mu, label_mu)
    x$locations <- lapply(x$locations, label_locations)
    rownames(x$partition$pi_weights) <- features
    names(x$partition$assign) <- features
    if (!is.null(x$lambda_mat)) rownames(x$lambda_mat) <- features
    if (!is.null(x$conditional_posterior)) {
      for (type in c("mean", "var", "diag")) {
        x$conditional_posterior[[type]] <- lapply(x$conditional_posterior[[type]], label_mu)
      }
      x$conditional_posterior$cov <- lapply(x$conditional_posterior$cov,
        function(cov) stats::setNames(cov, features))
      x$posterior <- x$conditional_posterior
    }
  } else {
    x$gamma <- label_gamma(x$gamma)
    x$params$mu <- label_mu(x$params$mu)
    x$locations <- label_locations(x$locations)
    if (!is.null(x$priors$position)) {
      label <- .mpcurve_result_orderings(x)
      x$priors$position$ordering_labels <- label
      names(x$priors$position$values) <- label
      names(x$priors$position$init) <- label
    }
  }
  if (!is.null(x$params$sigma2)) names(x$params$sigma2) <- features
  x
}

.mpcurve_result_orderings <- function(x) {
  if (!inherits(x, "mpcurve")) stop("x must be an 'mpcurve' object.", call. = FALSE)
  if (.mpcurve_is_partition(x)) return(names(x$gamma))
  names(x$dimension_estimation$expected_feature_counts) %||%
    .mpcurve_single_ordering_label()
}

#' Extract fitted sample positions
#'
#' Summarize each sample's variational position distribution on the normalized
#' grid from zero to one. The ordering dimension is retained for every fit,
#' including single-ordering and automatic fits.
#'
#' @param x An `mpcurve` fit.
#' @param type Posterior mean (`"mean"`), most probable grid position (`"map"`),
#'   posterior standard deviation (`"sd"`), or full grid probabilities
#'   (`"probability"`). MAP ties use the first grid position.
#' @return For mean, MAP, and SD, an `n x M` sample-by-ordering matrix. For
#'   probabilities, an `n x K x M` sample-by-bin-by-ordering array. Row names
#'   preserve input sample names; ordering labels are always supplied. Bin
#'   labels give the normalized grid positions. Dimensions are never dropped.
#' @details Pseudotime is a relative position, not elapsed time. The direction
#'   and labels of inferred orderings are arbitrary. SD describes variational
#'   position uncertainty in normalized pseudotime units.
#' @seealso [fitted_trajectories()], [fitted_assignments()], [fitted_prior()]
#' @examples
#' sim <- simulate_mpcurve(n = 30, d = 4, num_bins = 5, seed = 1)
#' fit <- fit_mpcurve(sim$X, num_bins = 5, max_iter = 2)
#' head(fitted_positions(fit))
#' dim(fitted_positions(fit, type = "probability"))
#' @md
#' @export
fitted_positions <- function(x, type = c("mean", "map", "sd", "probability")) {
  labels <- .mpcurve_result_orderings(x)
  type <- match.arg(type)
  gamma <- if (.mpcurve_is_partition(x)) x$gamma else list(x$gamma)
  grid <- seq(0, 1, length.out = x$K)
  if (type == "probability") {
    return(array(unlist(gamma, use.names = FALSE), c(x$n, x$K, length(labels)),
                 dimnames = list(rownames(x$data), as.character(grid), labels)))
  }
  values <- lapply(gamma, function(g) {
    if (type == "map") return(grid[max.col(g, ties.method = "first")])
    mean <- as.numeric(g %*% grid)
    if (type == "mean") return(mean)
    sqrt(pmax(0, as.numeric(g %*% grid^2) - mean^2))
  })
  matrix(unlist(values, use.names = FALSE), nrow = x$n,
         dimnames = list(rownames(x$data), labels))
}

#' Extract fitted feature trajectories
#'
#' Extract trajectory means or uncertainty on the normalized position grid.
#' For multiple orderings these are conditional on the feature following the
#' specified ordering, so interpret them alongside [fitted_assignments()].
#'
#' @param x An `mpcurve` fit.
#' @param type Posterior mean (`"mean"`), pointwise posterior standard deviation
#'   (`"sd"`), or full within-trajectory posterior covariance (`"covariance"`).
#' @return For mean and SD, a `d x K x M` feature-by-bin-by-ordering array.
#'   For covariance, a `K x K x d x M` array. Feature names preserve input
#'   column names, bin labels give normalized positions, and ordering labels
#'   match [fitted_positions()]. Dimensions are never dropped, including when
#'   `M = 1`. Means and SDs use the input feature units; covariances use squared
#'   units. For multiple orderings the distributions are `q(U_j | Z_j = m)`.
#' @seealso [fitted_positions()], [fitted_assignments()]
#' @md
#' @export
fitted_trajectories <- function(x, type = c("mean", "sd", "covariance")) {
  labels <- .mpcurve_result_orderings(x)
  type <- match.arg(type)
  partition <- .mpcurve_is_partition(x)
  posterior <- if (partition) x$conditional_posterior else x$fit$posterior
  field <- switch(type, mean = "mean", sd = "var", covariance = "cov")
  values <- if (partition) posterior[[field]] else list(posterior[[field]])
  if (is.null(posterior[[field]])) {
    stop("Trajectory posterior unavailable; refit with fit_mpcurve().", call. = FALSE)
  }
  grid <- as.character(seq(0, 1, length.out = x$K))
  if (type == "covariance") {
    return(array(unlist(values, use.names = FALSE), c(x$K, x$K, x$d, length(labels)),
                 dimnames = list(grid, grid, colnames(x$data), labels)))
  }
  out <- array(unlist(values, use.names = FALSE), c(x$d, x$K, length(labels)),
               dimnames = list(colnames(x$data), grid, labels))
  if (type == "sd") sqrt(pmax(out, 0)) else out
}

#' Extract feature-assignment probabilities
#'
#' Return the soft posterior probability that each feature follows each
#' ordering. These probabilities differ from the estimated partition prior
#' weights returned by [fitted_prior()].
#'
#' @param x An `mpcurve` fit.
#' @return A `d x M` feature-by-ordering matrix of probabilities, whose rows
#'   sum to one. With one ordering, returns one column of ones. Input feature
#'   names are preserved and ordering labels match [fitted_positions()].
#' @seealso [fitted_positions()], [fitted_trajectories()], [fitted_prior()]
#' @md
#' @export
fitted_assignments <- function(x) {
  labels <- .mpcurve_result_orderings(x)
  values <- if (.mpcurve_is_partition(x)) x$partition$pi_weights else rep(1, x$d)
  matrix(values, nrow = x$d, dimnames = list(colnames(x$data), labels))
}
