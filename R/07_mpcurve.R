# ============================================================
# mpcurve S3 class - unified model object for MPCurve fits
# ============================================================
# MPCurve is the statistical model (GMM + GMRF prior on trajectories), fitted
# by coordinate-ascent variational inference (CAVI).
#
# as_mpcurve()   - convert a raw algorithm fit to an mpcurve object
# fit_mpcurve()  - run the full pipeline and return an mpcurve
#
# Single-ordering objects expose matrix-valued params and gamma. Fixed-M
# structural partition objects instead expose named ordering lists for pi, mu,
# gamma, and locations, one shared sigma2 vector, conditional_posterior,
# lambda_mat, and partition$pi_weights. Their $fits list is derived for
# compatibility and is never continuation state. See ?mpcurve for the complete
# public schema.
# ============================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

.mpcurve_single_ordering_label <- function() "ordering1"

.mpcurve_ordering_labels <- function(M) {
  if (M <= 26L) LETTERS[seq_len(M)] else paste0("ord", seq_len(M))
}

.mpcurve_is_partition <- function(x) {
  if (inherits(x, "mpcurve")) {
    return(!is.null(x$partition) || inherits(x$fit, "soft_partition_cavi"))
  }
  inherits(x, "soft_partition_cavi")
}

.mpcurve_position_prior_mode <- function(control = NULL, default = "estimated") {
  ctl <- control %||% list()
  as.character(ctl$position_prior %||% default)[1]
}

.mpcurve_position_prior_bundle <- function(mode,
                                           values,
                                           init = NULL,
                                           ordering_labels = NULL) {
  if (is.null(ordering_labels)) {
    if (is.list(values) && !is.null(names(values))) {
      ordering_labels <- names(values)
    } else {
      ordering_labels <- .mpcurve_single_ordering_label()
    }
  }
  ordering_labels <- as.character(ordering_labels)

  if (!is.list(values)) {
    values <- stats::setNames(list(as.numeric(values)), ordering_labels[1L])
  } else {
    values <- stats::setNames(lapply(values, as.numeric), ordering_labels)
  }

  if (is.null(init)) {
    init <- values
  } else if (!is.list(init)) {
    init <- stats::setNames(list(as.numeric(init)), ordering_labels[1L])
  } else {
    init <- stats::setNames(lapply(init, as.numeric), ordering_labels)
  }

  list(
    mode = as.character(mode)[1],
    values = values,
    init = init,
    ordering_labels = ordering_labels
  )
}

.mpcurve_single_fit_position_priors <- function(fit,
                                                ordering_label = .mpcurve_single_ordering_label(),
                                                default_mode = "estimated") {
  current <- as.numeric((fit$params %||% list())$pi %||% numeric(0))
  init <- if (length(fit$pi_trace %||% list()) > 0L) {
    as.numeric(fit$pi_trace[[1L]])
  } else {
    current
  }
  .mpcurve_position_prior_bundle(
    mode = .mpcurve_position_prior_mode(fit$control %||% list(), default = default_mode),
    values = current,
    init = init,
    ordering_labels = ordering_label
  )
}

.mpcurve_position_priors_from_fits <- function(fits,
                                               ordering_labels,
                                               default_mode = "adaptive") {
  modes <- vapply(
    fits,
    function(fit) .mpcurve_position_prior_mode(fit$control %||% list(), default = default_mode),
    character(1)
  )
  mode <- unique(modes)
  mode <- if (length(mode) == 1L) mode else "mixed"
  values <- stats::setNames(
    lapply(fits, function(fit) as.numeric((fit$params %||% list())$pi %||% numeric(0))),
    ordering_labels
  )
  init <- stats::setNames(
    lapply(fits, function(fit) {
      if (length(fit$pi_trace %||% list()) > 0L) {
        as.numeric(fit$pi_trace[[1L]])
      } else {
        as.numeric((fit$params %||% list())$pi %||% numeric(0))
      }
    }),
    ordering_labels
  )
  .mpcurve_position_prior_bundle(
    mode = mode,
    values = values,
    init = init,
    ordering_labels = ordering_labels
  )
}

.mpcurve_partition_prior_mode <- function(control = NULL) {
  as.character((control %||% list())$partition_prior %||% "adaptive")[1]
}

.mpcurve_partition_prior_init <- function(control = NULL, M = NULL) {
  ctl <- control %||% list()
  init <- ctl$partition_prior_init %||% NULL
  mode <- .mpcurve_partition_prior_mode(ctl)
  if (is.null(init) && identical(mode, "fixed") && !is.null(M)) {
    init <- rep(1 / M, M)
  }
  if (is.null(init)) NULL else as.numeric(init)
}

.mpcurve_partition_prior_bundle <- function(mode,
                                            omega,
                                            init = NULL,
                                            ordering_labels = NULL) {
  ordering_labels <- as.character(ordering_labels %||% names(omega) %||% seq_along(omega))
  omega <- as.numeric(omega)
  if (length(omega) == length(ordering_labels)) {
    names(omega) <- ordering_labels
  }
  if (!is.null(init)) {
    init <- as.numeric(init)
    if (length(init) == length(ordering_labels)) {
      names(init) <- ordering_labels
    }
  }
  list(
    mode = as.character(mode)[1],
    omega = omega,
    init = init,
    ordering_labels = ordering_labels
  )
}

.mpcurve_validate_effective_count_tol <- function(tol, M, num_features) {
  if (!is.numeric(num_features) || length(num_features) != 1L ||
      !is.finite(num_features) || num_features < 1) {
    stop("X must contain at least one feature column.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) ||
      tol < 0 || tol >= num_features / M) {
    stop(
      "effective_count_tol must be a single finite number in [0, ncol(X) / intrinsic_dim).",
      call. = FALSE
    )
  }
  as.numeric(tol)
}

.mpcurve_effective_count_tol <- function(fit, M, num_features) {
  tol <- (fit$control %||% list())$effective_count_tol %||% 1e-8
  .mpcurve_validate_effective_count_tol(tol, M, num_features)
}

.mpcurve_partition_assignment_info_from_fit <- function(fit) {
  ctl <- fit$control %||% list()
  .cavi_partition_assignment_info(
    weights = fit$pi_weights,
    T_now = 1,
    partition_prior = ctl$partition_prior %||% "adaptive",
    partition_prior_init = ctl$partition_prior_init %||% NULL
  )
}

.mpcurve_priors_from_cavi_fit <- function(fit) {
  priors <- fit$priors %||% NULL
  if (!is.null(priors$position)) {
    return(priors)
  }
  list(
    position = .mpcurve_single_fit_position_priors(fit, default_mode = "adaptive"),
    partition = NULL
  )
}


.mpcurve_priors_from_partition_fit <- function(fit) {
  priors <- fit$priors %||% list()
  if (!is.null(priors$position) && !is.null(priors$partition)) {
    return(priors)
  }

  M <- fit$M %||% length(fit$fits)
  ordering_labels <- colnames(fit$pi_weights) %||% .mpcurve_ordering_labels(M)
  assignment_info <- .mpcurve_partition_assignment_info_from_fit(fit)
  partition_mode <- .mpcurve_partition_prior_mode(fit$control %||% list())

  list(
    position = .mpcurve_position_priors_from_fits(
      fits = fit$fits,
      ordering_labels = ordering_labels,
      default_mode = "adaptive"
    ),
    partition = .mpcurve_partition_prior_bundle(
      mode = partition_mode,
      omega = assignment_info$omega,
      init = .mpcurve_partition_prior_init(fit$control %||% list(), M = M),
      ordering_labels = ordering_labels
    )
  )
}

.mpcurve_get_priors <- function(x) {
  if (inherits(x, "mpcurve")) return(x$priors)
  if (inherits(x, "soft_partition_cavi")) return(.mpcurve_priors_from_partition_fit(x))
  if (inherits(x, "cavi")) return(.mpcurve_priors_from_cavi_fit(x))
  NULL
}


.mpcurve_resolve_ordering <- function(ordering, ordering_labels, caller = "fitted_prior()") {
  if (is.null(ordering)) return(NULL)
  if (length(ordering) != 1L || is.na(ordering)) {
    stop(caller, ": `ordering` must be a single numeric index or ordering label.", call. = FALSE)
  }
  if (is.numeric(ordering)) {
    idx <- as.integer(ordering)[1]
    if (!is.finite(idx) || idx < 1L || idx > length(ordering_labels)) {
      stop(caller, ": `ordering` index is out of range.", call. = FALSE)
    }
    return(ordering_labels[idx])
  }
  ordering <- as.character(ordering)[1]
  if (!ordering %in% ordering_labels) {
    stop(
      caller, ": unknown ordering `", ordering, "`. Available orderings are: ",
      paste(ordering_labels, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  ordering
}

.mpcurve_extract_fitted_prior <- function(priors,
                                          type = c("position", "partition"),
                                          ordering = NULL) {
  type <- match.arg(type)
  if (is.null(priors)) return(NULL)

  if (identical(type, "position")) {
    position <- priors$position %||% NULL
    if (is.null(position)) return(NULL)
    ordering_labels <- position$ordering_labels %||% names(position$values %||% list())
    selected <- .mpcurve_resolve_ordering(ordering, ordering_labels, caller = "fitted_prior()")
    if (is.null(selected)) {
      if (length(position$values) == 1L) {
        return(as.numeric(position$values[[1L]]))
      }
      return(position$values)
    }
    return(as.numeric(position$values[[selected]]))
  }

  partition <- priors$partition %||% NULL
  if (is.null(partition)) return(NULL)
  ordering_labels <- partition$ordering_labels %||% names(partition$omega %||% numeric(0))
  selected <- .mpcurve_resolve_ordering(ordering, ordering_labels, caller = "fitted_prior()")
  if (is.null(selected)) {
    return(partition)
  }
  omega <- partition$omega %||% numeric(0)
  if (!length(omega)) return(NULL)
  if (!is.null(names(omega))) {
    return(as.numeric(omega[[selected]]))
  }
  as.numeric(omega[match(selected, ordering_labels)])
}


.mpcurve_locations_from_gamma <- function(gamma) {
  if (is.null(gamma)) return(NULL)

  gamma <- as.matrix(gamma)
  if (!nrow(gamma) || !ncol(gamma)) return(NULL)

  K <- ncol(gamma)
  comp_grid <- seq_len(K)
  pseudo_grid <- if (K <= 1L) rep(0, K) else (comp_grid - 1L) / (K - 1L)

  map_index <- max.col(gamma, ties.method = "first")
  mean_index <- as.numeric(gamma %*% comp_grid)

  list(
    mean = list(
      index = mean_index,
      pseudotime = as.numeric(gamma %*% pseudo_grid)
    ),
    map = list(
      index = map_index,
      pseudotime = pseudo_grid[map_index]
    ),
    component_grid = comp_grid,
    pseudotime_grid = pseudo_grid
  )
}


# ---- as_mpcurve generic ----------------------------------------

#' Convert an algorithm fit object to an \code{mpcurve} model object
#'
#' @param x A \code{cavi} or \code{soft_partition_cavi} fit object.
#' @param ... Ignored.
#' @return An object of class \code{"mpcurve"}.
#'
#' For single-ordering fits, the returned object stores a unified \code{$params}
#' block, the responsibility matrix \code{$gamma}, and a \code{$locations} field
#' containing inferred cell locations derived from \code{$gamma}:
#' \itemize{
#'   \item \code{$locations$mean$index}: posterior-mean component index in \code{[1, K]}
#'   \item \code{$locations$mean$pseudotime}: responsibility-weighted pseudotime in \code{[0, 1]}
#'   \item \code{$locations$map$index}: MAP component index from \code{max.col(gamma)}
#'   \item \code{$locations$map$pseudotime}: pseudotime corresponding to the MAP component
#' }
#'
#' For partition fits (\code{intrinsic_dim >= 2}), \code{$locations} is a named
#' list with one such location object per ordering. Fits store their model
#' dimension in \code{$intrinsic_dim}; automatic adaptive fits remove empty
#' orderings before conversion.
#' @noRd
as_mpcurve <- function(x, ...) {
  ns <- asNamespace("MPCurver")
  classes <- class(x)

  for (cls in classes) {
    method <- get0(
      paste0("as_mpcurve.", cls),
      envir = ns,
      mode = "function",
      inherits = FALSE
    )
    if (is.function(method)) {
      return(.mpcurve_label_results(method(x, ...)))
    }
  }

  stop(
    "No as_mpcurve() method for object of class ",
    paste(sprintf("'%s'", classes), collapse = ", "),
    ".",
    call. = FALSE
  )
}


#' Extract fitted prior values
#'
#' @description
#' Extracts the fitted probability distribution over sample positions or
#' feature orderings. Position probabilities \eqn{\pi_k} describe the
#' population distribution across the latent grid. Partition probabilities
#' \eqn{\omega_m} describe the relative prevalence of the orderings among
#' features. These are the priors used to calculate sample-specific position
#' responsibilities and feature-specific assignment probabilities.
#'
#' Single-ordering position priors are length-\code{K} vectors. Multi-ordering
#' fits return a named list of \code{M} such vectors unless a specific ordering
#' is requested. For partition priors,
#' \code{ordering = NULL} returns the full prior record, including the prior
#' mode and effective \eqn{\omega} vector; a specific \code{ordering} returns
#' the corresponding scalar \eqn{\omega_m} mass.
#'
#' @param x A fitted \code{mpcurve}, raw \code{cavi}, or raw structural
#'   \code{soft_partition_cavi} object.
#' @param type One of \code{"position"} or \code{"partition"}.
#' @param ordering Optional ordering label or 1-based ordering index. For a
#'   single-ordering position prior it may identify the sole ordering; for a
#'   partition fit use the labels in \code{names(x$locations)}.
#' @param ... No additional arguments are accepted.
#'
#' @return If \code{type = "position"}, a numeric vector for a
#'   single-ordering fit, a named list of vectors for an unfiltered partition
#'   fit, or one numeric vector when \code{ordering} is supplied. If
#'   \code{type = "partition"}, \code{NULL} for a single-ordering fit, the
#'   full partition-prior record when \code{ordering = NULL}, or one unnamed
#'   numeric mass when an ordering is supplied.
#' @export
fitted_prior <- function(x,
                         type = c("position", "partition"),
                         ordering = NULL,
                         ...) {
  if (length(list(...))) stop("fitted_prior() does not accept additional arguments.", call. = FALSE)
  UseMethod("fitted_prior")
}


#' @rdname fitted_prior
#' @export
fitted_prior.mpcurve <- function(x,
                                 type = c("position", "partition"),
                                 ordering = NULL,
                                 ...) {
  .mpcurve_extract_fitted_prior(
    priors = .mpcurve_get_priors(x),
    type = type,
    ordering = ordering
  )
}


#' @rdname fitted_prior
#' @export
fitted_prior.cavi <- function(x,
                              type = c("position", "partition"),
                              ordering = NULL,
                              ...) {
  .mpcurve_extract_fitted_prior(
    priors = .mpcurve_get_priors(x),
    type = type,
    ordering = ordering
  )
}


#' @rdname fitted_prior
#' @export
fitted_prior.soft_partition_cavi <- function(x,
                                             type = c("position", "partition"),
                                             ordering = NULL,
                                             ...) {
  .mpcurve_extract_fitted_prior(
    priors = .mpcurve_get_priors(x),
    type = type,
    ordering = ordering
  )
}


#' @noRd
as_mpcurve.cavi <- function(x, ...) {
  params <- x$params
  mu_mat <- x$posterior$mean

  K <- length(params$pi)
  n <- if (!is.null(x$data)) nrow(x$data) else
    if (!is.null(x$gamma)) nrow(x$gamma) else NA_integer_
  d <- nrow(mu_mat)

  structure(
    list(
      data = x$data %||% NULL,
      params = list(
        pi     = as.numeric(params$pi),
        mu     = mu_mat,
        sigma2 = params$sigma2
      ),
      gamma        = x$gamma,
      measurement_sd = x$measurement_sd %||% NULL,
      locations    = .mpcurve_locations_from_gamma(x$gamma),
      elbo_trace   = x$elbo_trace   %||% numeric(0),
      loglik_trace = x$loglik_trace %||% numeric(0),
      lambda_trace = x$lambda_trace %||% list(),
      iter         = as.integer(x$iter %||% length(x$elbo_trace %||% numeric(0))),
      K            = as.integer(K),
      n            = as.integer(n),
      d            = as.integer(d),
      algorithm      = "cavi",
      modelName      = x$control$modelName %||% "homoskedastic",
      intrinsic_dim  = 1L,
      priors        = .mpcurve_priors_from_cavi_fit(x),
      init_info     = x$init_info %||% NULL,
      ordering_similarity = x$ordering_similarity %||% NULL,
      similarity_init = x$similarity_init %||% NULL,
      dimension_initialization = x$dimension_initialization %||% NULL,
      dimension_estimation = x$dimension_estimation %||% NULL,
      converged     = x$converged %||% NULL,
      convergence_info = x$convergence_info %||% NULL,
      control       = x$control %||% list(),
      fit            = x
    ),
    class = "mpcurve"
  )
}


#' @noRd
as_mpcurve.soft_partition_cavi <- function(x, ...) {
  M <- as.integer(x$M)
  ord_labels <- x$ordering_labels %||% colnames(x$pi_weights) %||%
    .mpcurve_ordering_labels(M)
  fits_mp <- stats::setNames(lapply(x$fits, as_mpcurve), ord_labels)
  locations <- stats::setNames(lapply(fits_mp, `[[`, "locations"), ord_labels)
  priors <- .mpcurve_priors_from_partition_fit(x)
  effective_count_tol <- .mpcurve_effective_count_tol(x, M, x$d)
  structure(
    list(
      data = x$data,
      params = list(
        pi = stats::setNames(x$position_pi, ord_labels),
        mu = stats::setNames(x$params$mu, ord_labels),
        sigma2 = x$params$sigma2
      ),
      gamma = stats::setNames(x$gamma, ord_labels),
      conditional_posterior = x$conditional_posterior,
      posterior = x$conditional_posterior,
      lambda_mat = x$lambda_mat,
      fits = fits_mp,
      measurement_sd = x$measurement_sd,
      locations = locations,
      elbo_trace = x$objective_history,
      loglik_trace = numeric(0),
      lambda_trace = x$lambda_trace,
      sigma2_trace = x$sigma2_trace,
      partition = list(
        pi_weights = x$pi_weights,
        assign = x$assign,
        variational_family = "structured"
      ),
      objective_history = x$objective_history,
      temperature_history = x$temperature_history,
      iter = x$iter,
      K = x$K,
      n = x$n,
      d = x$d,
      algorithm = "cavi",
      modelName = "partition_cavi_structured",
      intrinsic_dim = M,
      effective_count_tol = effective_count_tol,
      variational_family = "structured",
      priors = priors,
      init_info = x$init_info,
      ordering_similarity = x$ordering_similarity,
      similarity_init = x$similarity_init,
      dimension_initialization = x$dimension_initialization %||% NULL,
      dimension_estimation = x$dimension_estimation %||% NULL,
      converged = x$converged,
      convergence_info = x$convergence_info,
      control = x$control,
      fit = x
    ),
    class = "mpcurve"
  )
}


# ---- S3 methods ------------------------------------------------

#' Print an \code{mpcurve} fit
#'
#' Displays the model dimensions, iteration count, convergence status, and
#' final variational objective. Multi-ordering fits also show the fitted and
#' effective numbers of orderings and the feature-assignment counts. These
#' counts use the most probable assignment for each feature; inspect
#' \code{x$partition$pi_weights} to assess assignment uncertainty.
#'
#' @param x An \code{mpcurve} object.
#' @param ... Ignored.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.mpcurve <- function(x, ...) {
  idim <- as.integer(x$intrinsic_dim)
  is_partition <- .mpcurve_is_partition(x)
  selection_info <- x$dimension_selection %||% NULL
  initialization_info <- x$dimension_initialization %||% NULL

  if (is_partition) {
    cat("MPCurve partition fit\n")
    cat(sprintf("  Backend        : %s\n", x$algorithm))
    cat("  Variational    : structural q(C) q(Z) q(U | Z)\n")
  } else {
    cat("MPCurve fit\n")
    cat(sprintf("  Backend        : %s\n",  x$algorithm))
  }
  if (!is.null(selection_info)) {
    cat(sprintf(
      "  M selection   : %s (upper bound %d -> selected %d)\n",
      selection_info$direction,
      selection_info$max_intrinsic_dim,
      selection_info$selected_M
    ))
  }
  if (!is.null(initialization_info) && is.null(x$dimension_estimation)) {
    cat(sprintf(
      "  Auto M init    : mean silhouette (upper bound %d -> selected %d; min size %d)\n",
      initialization_info$requested_max_intrinsic_dim,
      initialization_info$selected_M,
      initialization_info$min_cluster_size
    ))
  }

  if (is_partition) {
    cat(sprintf("  Intrinsic dim  : %d\n", idim))
    cat(sprintf("  n / d / K      : %d / %d / %d\n", x$n, x$d, x$K))
    cat(sprintf("  Iterations     : %d\n", x$iter))
    part <- x$partition
    if (!is.null(part)) {
      tbl <- table(part$assign)
      part_str <- paste(sprintf("%s=%d", names(tbl), as.integer(tbl)), collapse = ", ")
      cat(sprintf("  Partition      : %s\n", part_str))
    }
    if (length(x$objective_history) > 0L)
      cat(sprintf("  Objective (last): %.6f\n",
                  tail(x$objective_history, 1L)))
    if (!is.null(x$converged)) {
      cat(sprintf("  Converged      : %s\n", x$converged))
    }
  } else {
    cat(sprintf("  Model          : %s\n",  x$modelName))
    cat(sprintf("  Intrinsic dim  : %d\n", idim))
    cat(sprintf("  n / d / K      : %d / %d / %d\n", x$n, x$d, x$K))
    cat(sprintf("  Iterations     : %d\n",  x$iter))
    if (!is.null(x$converged)) {
      cat(sprintf("  Converged      : %s\n", if (isTRUE(x$converged)) "yes" else "no"))
    }
    if (length(x$elbo_trace) > 0L)
      cat(sprintf("  ELBO (last)    : %.6f\n", tail(x$elbo_trace, 1L)))
  }
  invisible(x)
}


#' Summarise an \code{mpcurve} model fit
#'
#' Summarises model dimensions, fitted priors, and convergence diagnostics.
#' Multi-ordering summaries also contain feature-assignment probabilities,
#' the most probable assignment for each feature, and the noise specification.
#' The returned list can be used to construct comparison tables across fits.
#'
#' @param object An \code{mpcurve} object.
#' @param ... Passed to the underlying summary method for single-ordering fits;
#'   ignored for structural partition summaries.
#' @return An object of class \code{summary.mpcurve}. Every result contains
#'   \code{$algorithm}, \code{$modelName}, \code{$intrinsic_dim},
#'   \code{$dimension_selection} when returned by
#'   \code{select_mpcurve_dimension()},
#'   \code{$dimension_initialization} when \code{intrinsic_dim = "auto"} was
#'   used,
#'   \code{$dimension_estimation} for automatic dimension finalization,
#'   \code{$K}, \code{$n}, \code{$d}, \code{$priors}, and
#'   \code{$converged}. A single-ordering result also contains
#'   \code{$underlying}. A structural partition result instead contains
#'   \code{$partition}, \code{$objective_history},
#'   \code{$convergence_info}, \code{$control}, the shared
#'   \code{$sigma2} vector (or \code{NULL}), \code{$measurement_sd}, and
#'   \code{$variational_family}.
#' @export
summary.mpcurve <- function(object, ...) {
  idim <- as.integer(object$intrinsic_dim)
  priors <- .mpcurve_get_priors(object)

  if (.mpcurve_is_partition(object)) {
    result <- list(
      algorithm     = object$algorithm,
      intrinsic_dim = idim,
      effective_count_tol = object$effective_count_tol,
      modelName     = object$modelName,
      K             = object$K,
      n             = object$n,
      d             = object$d,
      priors        = priors,
      partition     = object$partition,
      objective_history = object$objective_history,
      converged     = object$converged,
      convergence_info = object$convergence_info,
      control       = object$control,
      sigma2        = object$params$sigma2,
      measurement_sd = object$measurement_sd,
      variational_family = object$variational_family,
      dimension_selection = object$dimension_selection %||% NULL,
      dimension_initialization = object$dimension_initialization %||% NULL,
      dimension_estimation = object$dimension_estimation %||% NULL
    )
  } else {
    underlying <- summary(object$fit, ...)
    result <- list(
      algorithm     = object$algorithm,
      modelName     = object$modelName,
      intrinsic_dim = idim,
      K             = object$K,
      n             = object$n,
      d             = object$d,
      priors        = priors,
      converged     = object$converged %||% underlying$converged %||% NULL,
      underlying    = underlying,
      dimension_selection = object$dimension_selection %||% NULL,
      dimension_initialization = object$dimension_initialization %||% NULL,
      dimension_estimation = object$dimension_estimation %||% NULL
    )
  }
  class(result) <- "summary.mpcurve"
  result
}


#' @rdname summary.mpcurve
#' @param x A \code{summary.mpcurve} object.
#' @export
print.summary.mpcurve <- function(x, ...) {
  idim <- as.integer(x$intrinsic_dim)
  selection_info <- x$dimension_selection %||% NULL
  initialization_info <- x$dimension_initialization %||% NULL

  if (!is.null(x$partition)) {
    cat("MPCurve Partition Summary\n")
    if (!is.null(selection_info)) {
      cat(sprintf("M selection : %s  |  upper bound = %d  |  selected = %d\n",
                  selection_info$direction,
                  selection_info$max_intrinsic_dim,
                  selection_info$selected_M))
    }
    if (!is.null(initialization_info) && is.null(x$dimension_estimation)) {
      cat(sprintf(
        "Auto M init : mean silhouette  |  upper bound = %d  |  selected = %d  |  min size = %d\n",
        initialization_info$requested_max_intrinsic_dim,
        initialization_info$selected_M,
        initialization_info$min_cluster_size
      ))
    }
    cat(sprintf("Algorithm : %s  |  structural VI  |  n=%d  d=%d  K=%d\n",
                x$algorithm, x$n, x$d, x$K))
    cat(sprintf("Intrinsic dimension : %d\n", idim))
    if (is.null(x$measurement_sd)) {
      cat(sprintf("Shared sigma2 range : [%.4g, %.4g]\n", min(x$sigma2), max(x$sigma2)))
    } else {
      cat("Noise model : known measurement SD\n")
    }
    tbl <- table(x$partition$assign)
    part_str <- paste(sprintf("%s=%d", names(tbl), as.integer(tbl)), collapse = "  ")
    cat(sprintf("Partition : %s\n", part_str))
    priors <- x$priors %||% NULL
    if (!is.null(priors$position)) {
      cat(sprintf("Position prior mode : %s\n", priors$position$mode))
    }
    if (!is.null(priors$partition)) {
      cat(sprintf("Partition prior mode : %s\n", priors$partition$mode))
      omega <- priors$partition$omega %||% numeric(0)
      if (length(omega) > 0L) {
        omega_str <- paste(sprintf("%s=%.3f", names(omega), omega), collapse = "  ")
        cat(sprintf("Effective omega : %s\n", omega_str))
      }
    }
    if (length(x$objective_history) > 0L)
      cat(sprintf("Objective (last) : %.6f\n", tail(x$objective_history, 1L)))
    cat("Objective type : fixed-M structural ELBO\n")
    if (!is.null(x$converged))
      cat(sprintf("Converged : %s\n", x$converged))
    conv_info <- x$convergence_info
    if (!is.null(conv_info)) {
      if (!is.null(conv_info$last_rel_delta) && is.finite(conv_info$last_rel_delta)) {
        cat(sprintf("Phase-2 rel_delta : %.3e\n", conv_info$last_rel_delta))
      }
      if (!is.null(conv_info$consecutive_small_steps)) {
        cat(sprintf("Small-step streak : %d\n", conv_info$consecutive_small_steps))
      }
      if (!is.null(conv_info$reason) && nzchar(conv_info$reason)) {
        cat(sprintf("Reason : %s\n", conv_info$reason))
      }
    }
  } else {
    cat("MPCurve Model Summary\n")
    if (!is.null(selection_info)) {
      cat(sprintf("M selection : %s  |  upper bound = %d  |  selected = %d\n",
                  selection_info$direction,
                  selection_info$max_intrinsic_dim,
                  selection_info$selected_M))
    }
    if (!is.null(initialization_info) && is.null(x$dimension_estimation)) {
      cat(sprintf(
        "Auto M init : mean silhouette  |  upper bound = %d  |  selected = %d  |  min size = %d\n",
        initialization_info$requested_max_intrinsic_dim,
        initialization_info$selected_M,
        initialization_info$min_cluster_size
      ))
    }
    cat(sprintf("Algorithm : %s  |  Model : %s  |  n=%d  d=%d  K=%d\n",
                x$algorithm, x$modelName, x$n, x$d, x$K))
    priors <- x$priors %||% NULL
    if (!is.null(priors$position)) {
      cat(sprintf("Position prior : %s\n", priors$position$mode))
      pi_now <- fitted_prior(structure(list(priors = priors), class = "mpcurve"), type = "position")
      if (length(pi_now) > 0L) {
        cat(sprintf("Current pi range : [%.4g, %.4g]\n", min(pi_now), max(pi_now)))
      }
    }
    if (!is.null(x$converged)) {
      cat(sprintf("Converged : %s\n", if (isTRUE(x$converged)) "yes" else "no"))
    }
    cat("Underlying fit summary\n")
    print(x$underlying, ...)
  }
  invisible(x)
}


# ---- internal gradient colour-bar legend ------------------------------------

# Draws a vertical gradient bar in the top-right corner of the current plot.
# Must be called after the main plot so par("usr") reflects real axis limits.
# R's arrow renderer cannot determine a direction below 0.001 inches.
.mpcurve_visible_segment <- function(from, to) {
  dx <- diff(graphics::grconvertX(c(from[1], to[1]), from = "user", to = "inches"))
  dy <- diff(graphics::grconvertY(c(from[2], to[2]), from = "user", to = "inches"))
  is.finite(dx) && is.finite(dy) && sqrt(dx^2 + dy^2) >= 1e-3
}

.mpcurve_pseudotime_colors <- function(pseudotime, pal) {
  index <- 1L + floor(pmax(0, pmin(1, pseudotime)) * (length(pal) - 1L))
  pal[index]
}

.mpcurve_plot_call <- function(fun, defaults, ...) {
  overrides <- list(...)
  if (length(overrides)) {
    if (is.null(names(overrides)) || any(!nzchar(names(overrides)))) {
      stop("Additional plot arguments must be named.", call. = FALSE)
    }
    defaults[names(overrides)] <- overrides
  }
  do.call(fun, defaults)
}

.draw_gradient_legend <- function(
    pal,
    title    = "pseudotime",
    lo_label = "0",
    hi_label = "1",
    n_rect   = 200L
) {
  usr <- par("usr")
  pw  <- usr[2] - usr[1]
  ph  <- usr[4] - usr[3]

  # bar: 3 % wide, 35 % tall, top-right with a 2 % margin
  bar_w  <- 0.03 * pw
  bar_h  <- 0.35 * ph
  margin <- 0.02

  xr <- usr[2] - margin * pw
  xl <- xr - bar_w
  yt <- usr[4] - margin * ph
  yb <- yt - bar_h

  # draw gradient (bottom = pseudotime 0, top = pseudotime 1)
  ys   <- seq(yb, yt, length.out = n_rect + 1L)
  cols <- pal[pmax(1L, round(seq(1L, length(pal), length.out = n_rect)))]

  old_xpd <- par(xpd = FALSE)
  on.exit(par(old_xpd), add = TRUE)

  for (i in seq_len(n_rect))
    rect(xl, ys[i], xr, ys[i + 1L], col = cols[i], border = NA)
  rect(xl, yb, xr, yt, border = "black", lwd = 0.5)

  # endpoint labels and title
  gap <- 0.008 * pw
  text(xr + gap, yb, lo_label, adj = c(0, 0.5), cex = 0.65)
  text(xr + gap, yt, hi_label, adj = c(0, 0.5), cex = 0.65)
  text((xl + xr) / 2, yt + 0.018 * ph, title, adj = c(0.5, 0), cex = 0.7)
}


#' Plot a 2D embedding with trajectory component means overlaid
#'
#' @description
#' Plot an \eqn{n \times 2} embedding (e.g., two selected dimensions, PCA/t-SNE/UMAP)
#' and overlay the trajectory component means (\code{mu_list}) as points connected by arrows.
#' Optionally color observations by a continuous vector (e.g., a Fiedler vector).
#'
#' @param mu_list A length-\eqn{K} list of mean vectors. Each element should be a numeric vector
#'   of length at least 2, where the first two entries correspond to the coordinates in \code{X2}.
#' @param X2 Numeric matrix of dimension \eqn{n \times 2} giving the embedding coordinates for observations.
#' @param t_vec Optional numeric vector of length \eqn{n}. If provided, points are colored by \code{t_vec}
#'   using \code{pal}.
#' @param mu_order Optional integer permutation of \code{1:K} specifying the plotting order of components.
#'   Defaults to \code{1:K}.
#' @param pch Point character for observations (passed to \code{\link[graphics]{points}} / \code{\link[graphics]{plot}}).
#' @param col Observation color used when \code{t_vec} is \code{NULL}.
#' @param cex Point size for observations.
#' @param mu_pch Point character for component means.
#' @param mu_col Color for component means.
#' @param mu_cex Point size for component means.
#' @param arrow_col Color for arrows connecting consecutive component means.
#' @param arrow_lwd Line width for arrows.
#' @param arrow_len Arrow head length (see \code{\link[graphics]{arrows}}).
#' @param add Logical; if \code{TRUE}, add to an existing plot. If \code{FALSE}, create a new plot.
#' @param xlab,ylab Axis labels. If \code{NULL}, defaults to \code{"dim 1"} and \code{"dim 2"}.
#' @param main Plot title. If \code{NULL}, uses a default title.
#' @param add_legend Logical; if \code{TRUE} and \code{t_vec} is provided, add a small legend (quantile ticks).
#' @param legend_title Title for the legend when \code{t_vec} is provided.
#' @param pal A vector of colors used to map \code{t_vec} to point colors. Defaults to a blue-white-red palette.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}} when \code{add = FALSE}.
#'
#' @return An invisible list with elements:
#' \itemize{
#'   \item \code{X2}: the embedding coordinates used for plotting.
#'   \item \code{mu}: the \eqn{K \times 2} matrix of plotted component means (after \code{mu_order}).
#'   \item \code{mu_order}: the component order used.
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 50
#' X2 <- cbind(rnorm(n), rnorm(n))
#' mu_list <- list(c(-1, -1), c(0, 0), c(1, 1))
#' t_vec <- rnorm(n)
#'
#' .mpcurve_plot_embedding_2d(mu_list, X2)
#' .mpcurve_plot_embedding_2d(mu_list, X2, t_vec = t_vec, legend_title = "Fiedler")
#'
#' @noRd
.mpcurve_plot_embedding_2d <- function(
    mu_list,
    X2,
    t_vec = NULL,
    mu_order = NULL,
    pch = 19,
    col = "grey60",
    cex = 0.7,
    mu_pch = 8,
    mu_col = "orange",
    mu_cex = 1.0,
    arrow_col = "orange",
    arrow_lwd = 2.5,
    arrow_len = 0.08,
    add = FALSE,
    xlab = NULL,
    ylab = NULL,
    main = NULL,
    add_legend = !is.null(t_vec),
    legend_title = "t",
    pal = grDevices::colorRampPalette(c("#2b8cbe", "white", "#de2d26"))(256),
    ...
) {
  X2 <- as.matrix(X2)
  if (!is.matrix(X2) || ncol(X2) != 2) stop("X2 must be an n x 2 matrix.")
  n <- nrow(X2)

  K <- length(mu_list)
  if (K < 1) stop("mu_list must be a non-empty list.")
  mu_mat <- do.call(rbind, mu_list)
  if (!is.matrix(mu_mat) || nrow(mu_mat) != K) stop("Cannot rbind mu_list into a K x d matrix.")
  if (ncol(mu_mat) < 2) stop("Each mu must have at least 2 entries.")

  mup <- mu_mat[, 1:2, drop = FALSE]

  if (is.null(mu_order)) mu_order <- seq_len(K)
  mu_order <- as.integer(mu_order)
  if (length(mu_order) != K || any(sort(mu_order) != seq_len(K))) {
    stop("mu_order must be a permutation of 1:K.")
  }
  mup <- mup[mu_order, , drop = FALSE]

  if (!is.null(t_vec)) {
    if (length(t_vec) != n) stop("t_vec must have length nrow(X2).")
    rng <- range(t_vec, na.rm = TRUE)
    z <- (t_vec - rng[1]) / (rng[2] - rng[1] + 1e-12)
    idx <- pmax(1L, pmin(256L, 1L + floor(z * 255)))
    pt_col <- pal[idx]
  } else {
    pt_col <- col
  }

  if (is.null(xlab)) xlab <- "dim 1"
  if (is.null(ylab)) ylab <- "dim 2"
  if (is.null(main)) main <- "2D embedding"

  if (!add) {
    graphics::plot(X2[, 1], X2[, 2],
                   pch = pch, col = pt_col, cex = cex,
                   xlab = xlab, ylab = ylab, main = main, ...
    )
  } else {
    graphics::points(X2[, 1], X2[, 2], pch = pch, col = pt_col, cex = cex, ...)
  }

  if (K >= 2) {
    for (k in 1:(K - 1)) {
      if (!.mpcurve_visible_segment(mup[k, ], mup[k + 1, ])) next
      graphics::arrows(
        mup[k, 1], mup[k, 2],
        mup[k + 1, 1], mup[k + 1, 2],
        col = arrow_col, lwd = arrow_lwd, length = arrow_len
      )
    }
  }
  graphics::points(mup[, 1], mup[, 2], pch = mu_pch, col = mu_col, cex = mu_cex)

  if (!is.null(t_vec) && isTRUE(add_legend)) {
    qs <- stats::quantile(t_vec, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)
    rng <- range(t_vec, na.rm = TRUE)
    zq <- (qs - rng[1]) / (rng[2] - rng[1] + 1e-12)
    idxq <- pmax(1L, pmin(256L, 1L + floor(zq * 255)))
    graphics::legend(
      "topright",
      legend = sprintf("%.3g", qs),
      col = pal[idxq], pch = 16, pt.cex = 1,
      bty = "n", title = legend_title
    )
  }

  invisible(list(X2 = X2, mu = mup, mu_order = mu_order))
}

#' @noRd
.plot_structural_mu_mpcurve <- function(x, dims, ...) {
  means <- x$conditional_posterior$mean
  if (!is.list(means) || length(means) < 1L) {
    stop(
      "x$conditional_posterior$mean must be a non-empty list of trajectory mean matrices.",
      call. = FALSE
    )
  }

  M <- as.integer(x$intrinsic_dim)
  if (length(M) != 1L || is.na(M) || M < 2L || length(means) != M) {
    stop(
      "x$conditional_posterior$mean must contain one matrix per fixed ordering.",
      call. = FALSE
    )
  }

  valid_matrix <- vapply(
    means,
    function(mu) is.matrix(mu) && is.numeric(mu) && all(is.finite(mu)),
    logical(1)
  )
  if (!all(valid_matrix)) {
    stop(
      "Every entry of x$conditional_posterior$mean must be a finite numeric matrix.",
      call. = FALSE
    )
  }

  d <- as.integer(x$d %||% nrow(means[[1L]]))
  K <- as.integer(x$K %||% ncol(means[[1L]]))
  if (length(d) != 1L || is.na(d) || d < 1L ||
      length(K) != 1L || is.na(K) || K < 2L) {
    stop("The structural state must contain positive d and K >= 2.", call. = FALSE)
  }
  expected_dims <- vapply(
    means,
    function(mu) identical(dim(mu), c(d, K)),
    logical(1)
  )
  if (!all(expected_dims) || K < 2L) {
    stop(
      sprintf(
        "Each conditional posterior mean must be a d x K matrix (%d x %d) with K >= 2.",
        d,
        K
      ),
      call. = FALSE
    )
  }

  dims <- as.integer(dims)
  if (length(dims) < 1L || length(dims) > 2L) {
    stop("dims must have length 1 or 2.", call. = FALSE)
  }
  if (anyNA(dims) || any(dims < 1L) || any(dims > d)) {
    stop("dims out of range for the structural trajectory means.", call. = FALSE)
  }

  pi_weights <- x$partition$pi_weights
  if (!is.matrix(pi_weights) || !is.numeric(pi_weights) ||
      !identical(dim(pi_weights), c(d, M)) ||
      any(!is.finite(pi_weights)) ||
      any(pi_weights < -1e-10) || any(pi_weights > 1 + 1e-10) ||
      any(abs(rowSums(pi_weights) - 1) > 1e-8)) {
    stop(
      sprintf(
        "x$partition$pi_weights must be a finite row-stochastic %d x %d matrix.",
        d,
        M
      ),
      call. = FALSE
    )
  }

  ordering_labels <- names(means)
  if (is.null(ordering_labels) || length(ordering_labels) != M ||
      anyNA(ordering_labels) || any(!nzchar(ordering_labels))) {
    ordering_labels <- colnames(pi_weights)
  }
  if (is.null(ordering_labels) || length(ordering_labels) != M ||
      anyNA(ordering_labels) || any(!nzchar(ordering_labels))) {
    ordering_labels <- .mpcurve_ordering_labels(M)
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  panel_cols <- if (M <= 3L) M else ceiling(sqrt(M))
  panel_rows <- ceiling(M / panel_cols)
  graphics::par(
    mfrow = c(panel_rows, panel_cols),
    mar = c(4, 4, 6, 1),
    cex.main = 0.85
  )

  positions <- (seq_len(K) - 1L) / (K - 1L)
  for (m in seq_len(M)) {
    mu_m <- means[[m]]
    weight_text <- paste(
      sprintf(
        "q(Z%d=%s)=%.2f",
        dims,
        ordering_labels[m],
        pi_weights[dims, m]
      ),
      collapse = "\n"
    )
    panel_title <- sprintf("Ordering %s\n%s", ordering_labels[m], weight_text)

    if (length(dims) == 2L) {
      .mpcurve_plot_call(graphics::plot, list(
        mu_m[dims[1L], ], mu_m[dims[2L], ],
        type = "o", pch = 16, col = "orange", lwd = 2,
        xlab = .mpcurve_feature_label(x, dims[1L]),
        ylab = .mpcurve_feature_label(x, dims[2L]),
        main = panel_title), ...)
    } else {
      .mpcurve_plot_call(graphics::plot, list(
        positions, mu_m[dims[1L], ],
        type = "o", pch = 16, col = "orange", lwd = 2,
        xlab = "pseudotime",
        ylab = .mpcurve_feature_label(x, dims[1L]),
        main = panel_title), ...)
    }
  }

  invisible(x)
}


#' Plot an \code{mpcurve} model fit
#'
#' @description
#' Visualise an \code{mpcurve} fit. Pseudotime is mapped to \eqn{[0,1]}, with
#' grid component \eqn{k} at \eqn{(k-1)/(K-1)}. For a responsibility row
#' \eqn{\gamma_i}, the posterior-mean pseudotime is
#' \deqn{t_i = \sum_{k=1}^{K} \gamma_{ik} \cdot \frac{k-1}{K-1}.}
#'
#' For a multi-ordering fit, scatterplots contain one
#' panel per ordering. Panel \eqn{m} uses that ordering's
#' \eqn{q(C^{(m)})} responsibilities for pseudotime and overlays the
#' conditional trajectory mean from \eqn{q(U_j\mid Z_j=m)}. The fitted
#' feature-assignment probabilities for the selected dimensions are shown in
#' each panel title. Point colors use the same absolute `[0,1]` pseudotime
#' scale in every panel, without stretching each panel's observed range.
#'
#' \describe{
#'   \item{\code{"scatterplot"}}{Scatter of two chosen dimensions, colored by
#'     pseudotime, with posterior means overlaid as an orange path. For a
#'     multi-ordering fit, draws one panel per ordering. With one dimension,
#'     displays the observed feature values against inferred pseudotime.}
#'   \item{\code{"elbo"}}{Variational objective over fitting iterations.
#'     Fitting uses \eqn{T=1} throughout by default. If optional annealing is
#'     enabled, assess convergence within the final \eqn{T=1} segment, since
#'     changing the temperature changes the objective.}
#'   \item{\code{"mu"}}{Posterior-mean trajectory only. With one selected
#'     dimension, plots the feature mean against normalized latent position;
#'     with two dimensions, plots the posterior-mean path in that feature
#'     plane. Multi-ordering fits draw one panel per ordering, showing
#'     \eqn{E[U_j\mid Z_j=m]}. Interpret each conditional trajectory alongside
#'     its feature-assignment probability, shown in the panel title.}
#' }
#'
#' @param x An \code{mpcurve} object.
#' @param plot_type One of \code{"scatterplot"} (default), \code{"elbo"},
#'   \code{"mu"}.
#' @param dims One or two feature column indices or names. Determines the
#'   scatterplot axes and trajectory coordinates for \code{plot_type = "mu"}.
#'   \code{NULL} selects the first two features, or the only feature in a
#'   one-feature fit. Available feature names are used as axis labels.
#' @param data Optional numeric \code{n x d} matrix used by scatterplots. If
#'   \code{NULL}, data stored on \code{x} are used. It is not required for
#'   \code{plot_type = "elbo"} or \code{plot_type = "mu"}.
#' @param pal Colour palette (length-256 character vector) used to map
#'   pseudotime to point colours. Defaults to a rainbow palette
#'   (blue \eqn{\to} cyan \eqn{\to} green \eqn{\to} yellow \eqn{\to} red).
#' @param add_legend Logical; draw a gradient color-bar legend?  Default
#'   \code{TRUE}.
#' @param ... Additional base-graphics arguments. For single-ordering
#'   non-scatter plots these are forwarded to the underlying CAVI plot method.
#'
#' @return Invisibly returns \code{x}.
#' @examples
#' \dontrun{
#' sim <- simulate_dual_trajectory(n = 100, d1 = 3, d2 = 3, seed = 1)
#' fit <- fit_mpcurve(
#'   sim$X,
#'   intrinsic_dim = 2,
#'   num_bins = 8,
#'   max_iter = 3
#' )
#' plot(fit, plot_type = "scatterplot", dims = c(1, 4))
#' plot(fit, plot_type = "mu", dims = 1)
#' plot(fit, plot_type = "mu", dims = c(1, 4))
#' plot(fit, plot_type = "elbo")
#' }
#' @export
plot.mpcurve <- function(
    x,
    plot_type  = c("scatterplot", "elbo", "mu"),
    dims       = NULL,
    data       = NULL,
    pal        = grDevices::colorRampPalette(
                   c("#0000FF", "#00FFFF", "#00FF00",
                     "#FFFF00", "#FF0000"))(256L),
    add_legend = TRUE,
    ...
) {
  if (!inherits(x, "mpcurve")) stop("x must be an 'mpcurve' object.")
  plot_type <- match.arg(plot_type)
  if (plot_type != "elbo") {
    dims <- .mpcurve_plot_dimensions(dims, x$d, colnames(x$data))
  }
  is_partition <- .mpcurve_is_partition(x)

  # ---- Multi-ordering partition plot ----
  if (is_partition && plot_type == "scatterplot") {
    M <- x$intrinsic_dim
    fits_list <- x$fits
    pi_w <- x$partition$pi_weights   # d x M

    # Resolve data
    if (is.null(data)) {
      data <- x$data %||% fits_list[[1]]$data %||% fits_list[[1]]$fit$data
      if (is.null(data))
        stop("No data found. Please supply `data` explicitly.")
    }
    data <- as.matrix(data)

    dims <- as.integer(dims)
    if (length(dims) < 1L || length(dims) > 2L) stop("dims must have length 1 or 2.")
    if (any(dims < 1L) || any(dims > ncol(data)))
      stop("dims out of range for the number of columns in data.")

    ord_labels <- colnames(pi_w)
    if (is.null(ord_labels)) {
      ord_labels <- if (M <= 26L) LETTERS[seq_len(M)] else paste0("ord", seq_len(M))
    }

    old_par <- par(mfrow = c(1, M), mar = c(4, 4, 3.5, 1))
    on.exit(par(old_par), add = TRUE)

    for (m in seq_len(M)) {
      sub_fit <- fits_list[[m]]
      K_m <- sub_fit$K
      gamma_m <- sub_fit$gamma
      positions_m <- (seq_len(K_m) - 1L) / (K_m - 1L)
      t_pseudo_m <- as.numeric(gamma_m %*% positions_m)

      # Colour by pseudotime
      pt_col <- .mpcurve_pseudotime_colors(t_pseudo_m, pal)

      # Weight annotation for plotted dims
      w_dims <- pi_w[dims, m]
      w_str <- paste(sprintf("d%d:w=%.2f", dims, w_dims), collapse = ", ")
      main_m <- sprintf("Ordering %s  (%s)", ord_labels[m], w_str)

      mu_m <- sub_fit$params$mu   # d x K

      if (length(dims) == 2L) {
        .mpcurve_plot_call(plot, list(data[, dims[1]], data[, dims[2]],
             pch = 19, col = pt_col, cex = 0.6,
             main = main_m,
             xlab = .mpcurve_feature_label(x, dims[1]),
             ylab = .mpcurve_feature_label(x, dims[2])), ...)
        # Overlay component means with arrows
        for (k in 2:K_m) {
          if (!.mpcurve_visible_segment(mu_m[dims, k - 1], mu_m[dims, k])) next
          arrows(mu_m[dims[1], k - 1], mu_m[dims[2], k - 1],
                 mu_m[dims[1], k], mu_m[dims[2], k],
                 col = "orange", lwd = 1.5, length = 0.06)
        }
        points(mu_m[dims[1], ], mu_m[dims[2], ],
               pch = 8, col = "orange", cex = 0.9)
      } else {
        j <- dims[1]
        .mpcurve_plot_call(plot, list(t_pseudo_m, data[, j],
             pch = 19, col = pt_col, cex = 0.6,
             main = main_m,
             xlab = "pseudotime",
             ylab = .mpcurve_feature_label(x, j)), ...)
        mu_j <- mu_m[j, ]
        lines(positions_m, mu_j, col = "orange", lwd = 2)
        points(positions_m, mu_j, pch = 8, col = "orange", cex = 0.9)
      }

      if (add_legend)
        .draw_gradient_legend(pal, title = "pseudotime",
                              lo_label = "0", hi_label = "1")
    }

    return(invisible(x))
  }

  if (is_partition && plot_type == "mu") {
    return(.plot_structural_mu_mpcurve(x, dims = dims, ...))
  }

  # ---- plot objective traces or delegate single-ordering diagnostics ----
  if (plot_type != "scatterplot") {
    if (is_partition) {
      # The partition ELBO is stored only in the canonical structural state.
      obj <- x$objective_history
      if (length(obj) > 0L) {
        .mpcurve_plot_call(plot, list(obj, type = "b", pch = 19, cex = 0.7,
             xlab = "Iteration", ylab = "Fixed-M structural objective",
             main = "Fixed-M structural objective trace"), ...)
      }
    } else {
      .mpcurve_plot_call(plot, list(
        x$fit,
        plot_type = plot_type,
        dims = dims,
        data = data,
        pal = pal,
        add_legend = add_legend), ...)
    }
    return(invisible(x))
  }

  # ---- resolve data ----
  if (is.null(data)) {
    data <- x$data %||% x$fit$data
    if (is.null(data))
      stop("No data found in x$data. Please supply `data` explicitly.")
  }
  data <- as.matrix(data)

  dims <- as.integer(dims)
  if (length(dims) < 1L || length(dims) > 2L) stop("dims must have length 1 or 2.")
  if (any(dims < 1L) || any(dims > ncol(data)))
    stop("dims out of range for the number of columns in data.")

  # ---- pseudotime: responsibility-weighted component position in [0, 1] ----
  K     <- x$K
  gamma <- x$gamma   # n x K
  if (is.null(gamma)) stop("x$gamma is NULL; cannot compute pseudotime.")

  positions <- (seq_len(K) - 1L) / (K - 1L)   # 0, 1/(K-1), ..., 1
  t_pseudo  <- as.numeric(gamma %*% positions)  # length n

  # ---- component means projected onto dims ----
  mu_mat <- x$params$mu   # d x K

  # ---- title ----
  main <- sprintf("MPCurve pseudotime  (algorithm: %s, K=%d)",
                  x$algorithm, K)

  # ---- map pseudotime -> color ----
  pt_col <- .mpcurve_pseudotime_colors(t_pseudo, pal)

  if (length(dims) == 2L) {
    # ===== 2D scatter =====
    mu_list_dims <- lapply(seq_len(K), function(k) mu_mat[dims, k])

    .mpcurve_plot_call(.mpcurve_plot_embedding_2d, list(
      mu_list    = mu_list_dims,
      X2         = data[, dims, drop = FALSE],
      t_vec      = NULL,
      col        = pt_col,
      pal        = pal,
      add_legend = FALSE,
      main       = main,
      xlab       = .mpcurve_feature_label(x, dims[1]),
      ylab       = .mpcurve_feature_label(x, dims[2])), ...)

  } else {
    # ===== 1D: pseudotime (x) vs selected dimension (y) =====
    j    <- dims[1]
    xlab <- "pseudotime"
    ylab <- .mpcurve_feature_label(x, j)

    mu_j <- mu_mat[j, ]               # length K
    mu_t <- positions                  # pseudotime of each component (0..1)

    .mpcurve_plot_call(graphics::plot, list(t_pseudo, data[, j],
                   pch = 19, col = pt_col, cex = 0.7,
                   xlab = xlab, ylab = ylab, main = main), ...)
    graphics::lines(mu_t, mu_j, col = "orange", lwd = 2)
    graphics::points(mu_t, mu_j, pch = 8, col = "orange", cex = 1)
  }

  if (add_legend)
    .draw_gradient_legend(pal, title = "pseudotime", lo_label = "0", hi_label = "1")

  invisible(x)
}


# ---- do_mpcurve ------------------------------------------------

#' Continue a fitted MPCurve model
#'
#' Resume coordinate-ascent variational inference from the stored data,
#' posterior, and parameters. Measurement errors, trajectory priors, and
#' initialization are inherited. Partition inference resumes at temperature 1.
#'
#' @param object An `mpcurve` object returned by [fit_mpcurve()].
#' @param max_iter Positive integer maximum number of additional CAVI sweeps.
#' @param tol Optional nonnegative stopping tolerance. `NULL` inherits the
#'   stored tolerance; zero disables early stopping. See [fit_mpcurve()].
#' @param verbose Print fitting progress?
#' @param control Settings from [mpcurve_continue_control()] or a named list
#'   of its options. Omitted or `NULL` settings inherit the fitted values.
#' @details
#' Objective and parameter traces are extended. Changing precision, its prior,
#' or parameter bounds can change the objective or feasible parameter range;
#' assess convergence within the new run. Use [fit_mpcurve()] to change the
#' random-walk order, ridge, data, known measurement errors, or prior modes.
#' Model selection records are retained and marked as continued. Automatic
#' dimension estimates retain their estimation records and remove any new
#' empty orderings after continuation without additional fitting sweeps.
#' @return An updated `mpcurve` object.
#' @examples
#' sim <- simulate_spiral2d(n = 40, seed = 1)
#' fit <- fit_mpcurve(sim$obs, num_bins = 5, max_iter = 2)
#' fit <- do_mpcurve(fit, max_iter = 3)
#' fixed <- do_mpcurve(fit, max_iter = 3,
#'   control = mpcurve_continue_control(lambda_init = 5, fix_lambda = TRUE))
#' @md
#' @export
do_mpcurve <- function(object, max_iter = 1L, tol = NULL,
                       verbose = FALSE, control = NULL) {
  if (!inherits(object, "mpcurve")) stop("object must be an 'mpcurve' object.")
  max_iter <- .mpcurve_integer_setting(max_iter, "max_iter")
  if (!is.null(tol)) {
    .mpcurve_positive_setting(tol, "tol", allow_zero = TRUE)
    if (length(tol) != 1L) stop("tol must be a scalar.", call. = FALSE)
  }
  verbose <- .mpcurve_logical_setting(verbose, "verbose")
  ctrl <- .mpcurve_interface_options(control, mpcurve_continue_control, "control")
  raw <- object$fit
  if (!inherits(raw, c("cavi", "soft_partition_cavi"))) {
    stop("object$fit must be a cavi or soft_partition_cavi object.")
  }
  partition <- inherits(raw, "soft_partition_cavi")
  if (!is.null(ctrl$fix_lambda)) raw$control$fix_lambda <- ctrl$fix_lambda
  bounds <- ctrl$lambda_bounds %||%
    c(raw$control$lambda_min %||% 1e-10, raw$control$lambda_max %||% 1e10)
  if (partition) {
    updated <- .continue_structural_partition_cavi(
      fit = raw, iter = max_iter, lambda = ctrl$lambda_init, tol_outer = tol,
      convergence = ctrl$convergence,
      lambda_sd_prior_rate = ctrl$lambda_sd_prior_rate,
      lambda_min = bounds[1], lambda_max = bounds[2],
      sigma_min = if (!is.null(ctrl$sigma2_bounds)) ctrl$sigma2_bounds[1],
      sigma_max = if (!is.null(ctrl$sigma2_bounds)) ctrl$sigma2_bounds[2],
      verbose = verbose)
  } else {
    # Precision overrides reset the value; only fix_lambda changes its update mode.
    if (!is.null(ctrl$lambda_init)) {
      d <- ncol(raw$data)
      if (!(length(ctrl$lambda_init) %in% c(1L, d)) ||
          !is.null(dim(ctrl$lambda_init))) {
        stop("lambda_init must be a scalar or one value per feature.", call. = FALSE)
      }
      raw$lambda_vec <- rep(as.numeric(ctrl$lambda_init), length.out = d)
    }
    raw$lambda_vec <- pmax(bounds[1], pmin(bounds[2], raw$lambda_vec))
    updated <- do_cavi(
      object = raw, iter = max_iter, tol = tol,
      convergence = ctrl$convergence,
      lambda_sd_prior_rate = ctrl$lambda_sd_prior_rate,
      lambda_min = bounds[1], lambda_max = bounds[2],
      sigma_min = if (!is.null(ctrl$sigma2_bounds)) ctrl$sigma2_bounds[1],
      sigma_max = if (!is.null(ctrl$sigma2_bounds)) ctrl$sigma2_bounds[2],
      verbose = verbose)
  }
  updated$dimension_initialization <- object$dimension_initialization %||%
    object$fit$dimension_initialization %||% NULL
  updated$dimension_estimation <- object$dimension_estimation %||%
    object$fit$dimension_estimation %||% NULL
  updated$control$effective_count_tol <- raw$control$effective_count_tol
  updated <- .mpcurve_finalize_automatic_dimension(updated)
  out <- as_mpcurve(updated)
  out$dimension_selection <- object$dimension_selection %||% NULL
  if (!is.null(out$dimension_selection)) {
    out$dimension_selection$continued_after_selection <- TRUE
  }
  out$dimension_initialization <- updated$dimension_initialization
  out$continuation_history <- c(object$continuation_history %||% list(), list(list(
    start_trace_length = length(object$elbo_trace), max_iter = max_iter,
    tol = tol, control = ctrl)))
  out
}


# ---- fit_mpcurve -----------------------------------------------

.fit_mpcurve <- function(
    X,
    method = "PCA",
    K = NULL,
    rw_q = 2,
    lambda = 1,
    S = NULL,
    sigma2_init = NULL,
    fix_lambda = FALSE,
    iter = 100,
    tol = 1e-6,
    intrinsic_dim = 1L,
    discretization = NULL,
    ridge = 0,
    lambda_sd_prior_rate = NULL,
    lambda_min = 1e-10,
    lambda_max = 1e10,
    sigma_min = 1e-10,
    sigma_max = 1e10,
    position_prior = c("adaptive", "fixed"),
    position_prior_init = NULL,
    partition_prior = c("adaptive", "fixed"),
    partition_prior_init = NULL,
    effective_count_tol = 1e-8,
    similarity_metric = c("spearman", "pearson", "smooth_fit", "spline_r2"),
    smooth_fit_lambda_mode = c("optimize", "fixed"),
    smooth_fit_lambda_value = 1,
    cluster_linkage = "single",
    similarity_min_feature_sd = 1e-8,
    T_start = 5,
    T_end = 1,
    n_outer = 0L,
    inner_iter = 1L,
    max_converge_iter = NULL,
    tol_outer = 1e-6,
    verbose = FALSE,
    convergence = c("normalized", "relative"),
    max_intrinsic_dim = 8L,
    spline_r2_df = 5L,
    similarity_min_cluster_size = 2L,
    ...
) {
  convergence <- match.arg(convergence)
  similarity_metric_missing <- missing(similarity_metric)
  automatic_dimension <- is.character(intrinsic_dim) &&
    length(intrinsic_dim) == 1L && !is.na(intrinsic_dim) &&
    identical(intrinsic_dim, "auto")
  if (is.character(intrinsic_dim) && !automatic_dimension) {
    stop("intrinsic_dim must be a positive integer or \"auto\".", call. = FALSE)
  }
  if (!automatic_dimension) {
    intrinsic_dim <- .cavi_validate_positive_integer(
      intrinsic_dim, "intrinsic_dim"
    )
  }
  position_prior <- match.arg(position_prior)
  partition_prior <- match.arg(partition_prior)
  similarity_metric <- if (automatic_dimension && similarity_metric_missing) {
    "spline_r2"
  } else {
    match.arg(similarity_metric)
  }
  smooth_fit_lambda_mode <- match.arg(smooth_fit_lambda_mode)
  lambda_sd_prior_rate <- .normalize_lambda_sd_prior_rate(lambda_sd_prior_rate)
  dots <- list(...)
  reserved_initialization_args <- intersect(
    names(dots),
    c("similarity_precomputed", "initial_partition_probabilities",
      "dimension_initialization")
  )
  if (length(reserved_initialization_args)) {
    stop(
      paste(reserved_initialization_args, collapse = ", "),
      " are internal initialization arguments and cannot be supplied through the ellipsis.",
      call. = FALSE
    )
  }
  max_converge_iter <- max_converge_iter %||% as.integer(iter)

  similarity_precomputed <- NULL
  initial_partition_probabilities <- NULL
  dimension_initialization <- NULL
  if (automatic_dimension) {
    if ("fits_init" %in% names(dots)) {
      stop(
        "fits_init cannot be supplied when intrinsic_dim = \"auto\" because ",
        "the selected feature cut defines the ordering initializations.",
        call. = FALSE
      )
    }
    max_intrinsic_dim <- .cavi_validate_positive_integer(
      max_intrinsic_dim, "max_intrinsic_dim"
    )
    similarity_min_cluster_size <- .cavi_validate_positive_integer(
      similarity_min_cluster_size, "similarity_min_cluster_size"
    )
    cluster_linkage <- .cavi_validate_cluster_linkage(cluster_linkage)

    X_similarity <- as.matrix(X)
    discretization_similarity <- match.arg(
      discretization %||% "quantile",
      c("quantile", "equal", "kmeans")
    )
    K_similarity <- .cavi_resolve_K(X_similarity, K, rw_q)
    similarity_precomputed <- .compute_same_ordering_similarity(
      X = X_similarity,
      S = S,
      metric = similarity_metric,
      min_feature_sd = similarity_min_feature_sd,
      spline_r2_df = spline_r2_df,
      K = K_similarity,
      rw_q = rw_q,
      ridge = ridge,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      smooth_fit_lambda_mode = smooth_fit_lambda_mode,
      smooth_fit_lambda_value = smooth_fit_lambda_value,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      discretization = discretization_similarity
    )
    dimension_initialization <- .cavi_select_similarity_dimension(
      distance = similarity_precomputed$distance,
      cluster_linkage = cluster_linkage,
      max_intrinsic_dim = max_intrinsic_dim,
      min_cluster_size = similarity_min_cluster_size
    )
    intrinsic_dim <- dimension_initialization$selected_M
    cluster_linkage <- dimension_initialization$cluster_linkage
    if (intrinsic_dim >= 2L && identical(partition_prior, "adaptive")) {
      initial_partition_probabilities <-
        dimension_initialization$cluster_proportions
    }
    actual_initial_probabilities <- if (intrinsic_dim == 1L) {
      stats::setNames(1, .cavi_partition_order_labels(1L))
    } else if (is.null(initial_partition_probabilities)) {
      stats::setNames(
        rep(1 / intrinsic_dim, intrinsic_dim),
        .cavi_partition_order_labels(intrinsic_dim)
      )
    } else {
      initial_partition_probabilities
    }
    dimension_initialization$automatic <- TRUE
    dimension_initialization$similarity_metric <- similarity_metric
    dimension_initialization$spline_r2_df <- if (
      identical(similarity_metric, "spline_r2")
    ) {
      similarity_precomputed$spline_r2_df %||% as.integer(spline_r2_df)[1]
    } else NULL
    dimension_initialization$initial_partition_probabilities <-
      actual_initial_probabilities
    dimension_initialization$similarity <- similarity_precomputed$S
    dimension_initialization$distance <- similarity_precomputed$distance

  }
  effective_count_tol <- .mpcurve_validate_effective_count_tol(
    effective_count_tol, intrinsic_dim, ncol(X)
  )

  # ---- Partition model (intrinsic_dim >= 2) ----
  if (intrinsic_dim >= 2L) {
    M <- intrinsic_dim

    init_methods <- rep(method, M)
    pca_components <- NULL

    if ("pca_components" %in% names(dots)) {
      pca_components <- dots$pca_components
      dots$pca_components <- NULL
    }

    sp_args <- c(
      list(
        X = X,
        S = S,
        M = M,
        init_methods = init_methods,
        pca_components = pca_components,
        similarity_metric = similarity_metric,
        spline_r2_df = spline_r2_df,
        smooth_fit_lambda_mode = smooth_fit_lambda_mode,
        smooth_fit_lambda_value = smooth_fit_lambda_value,
        cluster_linkage = cluster_linkage,
        similarity_min_feature_sd = similarity_min_feature_sd,
        similarity_min_cluster_size = similarity_min_cluster_size,
        similarity_precomputed = similarity_precomputed,
        initial_partition_probabilities = initial_partition_probabilities,
        K = K,
        rw_q = rw_q,
        lambda_init = lambda,
        fix_lambda = fix_lambda,
        discretization = discretization %||% c("quantile", "equal", "kmeans")[1L],
        T_start = T_start,
        T_end = T_end,
        n_outer = n_outer,
        inner_iter = inner_iter,
        max_converge_iter = max_converge_iter,
        tol_outer = tol_outer,
        convergence = convergence,
        ridge = ridge,
        lambda_sd_prior_rate = lambda_sd_prior_rate,
        lambda_min = lambda_min,
        lambda_max = lambda_max,
        sigma_min = sigma_min,
        sigma_max = sigma_max,
        sigma2_init = sigma2_init,
        position_prior = position_prior,
        position_prior_init = position_prior_init,
        partition_prior = partition_prior,
        partition_prior_init = partition_prior_init,
        verbose = verbose
      ),
      dots
    )
    raw <- do.call(soft_partition_cavi, sp_args)
    raw$control$effective_count_tol <- effective_count_tol
    if (automatic_dimension) {
      raw$control$max_intrinsic_dim <- max_intrinsic_dim
      raw$dimension_initialization <- dimension_initialization
      raw <- .mpcurve_finalize_automatic_dimension(raw)
    }
    out <- as_mpcurve(raw)
    if (automatic_dimension) {
      out$dimension_initialization <- dimension_initialization
      out$fit$dimension_initialization <- dimension_initialization
    }
    return(out)
  }

  # Automatic initialization can fall back to one ordering.
  if ("pca_components" %in% names(dots)) {
    dots$method_args$component <- .mpcurve_integer_setting(dots$pca_components, "pca_components")
    dots$pca_components <- NULL
  }
  run_one_cavi <- function(method_i) {
    cavi_args <- c(
      list(
        X = X,
        K = K,
        method = method_i,
        rw_q = rw_q,
        ridge = ridge,
        discretization = discretization,
        S = S,
        sigma2_init = sigma2_init,
        fix_lambda = fix_lambda,
        lambda_sd_prior_rate = lambda_sd_prior_rate,
        lambda_min = lambda_min,
        lambda_max = lambda_max,
        sigma_min = sigma_min,
        sigma_max = sigma_max,
        position_prior = position_prior,
        position_prior_init = position_prior_init,
        max_iter = iter,
        tol = tol,
        convergence = convergence,
        verbose = verbose
      ),
      dots
    )
    if (!("lambda_init" %in% names(cavi_args))) {
      cavi_args$lambda_init <- lambda
    }
    fit_error <- NULL
    fit_result <- tryCatch(
      do.call(cavi, cavi_args),
      error = function(e) {
        fit_error <<- e
        NULL
      }
    )
    list(
      result = fit_result,
      error = fit_error,
      method = method_i
    )
  }

  fit_info <- run_one_cavi(method[[1L]])
  if (is.null(fit_info$result)) {
    if (!is.null(fit_info$error)) {
      stop(
        sprintf(
          "fit_mpcurve() failed while fitting cavi with method='%s': %s",
          fit_info$method,
          conditionMessage(fit_info$error)
        ),
        call. = FALSE
      )
    }
    stop("fit_mpcurve() failed to produce a valid cavi fit for unknown reasons.", call. = FALSE)
  }

  if (automatic_dimension) {
    fit_info$result$dimension_initialization <- dimension_initialization
    fit_info$result$control$effective_count_tol <- effective_count_tol
    fit_info$result$control$max_intrinsic_dim <- max_intrinsic_dim
    fit_info$result$control$similarity_metric <- similarity_metric
    fit_info$result$control$spline_r2_df <- as.integer(spline_r2_df)[1]
    fit_info$result$control$similarity_min_cluster_size <-
      as.integer(similarity_min_cluster_size)[1]
    fit_info$result <- .mpcurve_finalize_automatic_dimension(fit_info$result)
  }
  out <- as_mpcurve(fit_info$result)
  if (automatic_dimension) out$dimension_initialization <- dimension_initialization
  out
}
