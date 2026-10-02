# ============================================================
# mpcurve S3 class - unified model object for MPCurve fits
# ============================================================
# MPCurve is the statistical model (GMM + GMRF prior on trajectories).
# The public fitting path is CAVI. csmooth_em and smooth_em remain available
# only for internal normalization, compatibility, and benchmarking.
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
  ctl <- control %||% list()
  if (identical(ctl$assignment_mode %||% NULL, "legacy_dirichlet") ||
      identical(ctl$assignment_prior %||% NULL, "dirichlet")) {
    return("legacy_dirichlet")
  }
  if (!is.null(ctl$partition_prior)) {
    return(as.character(ctl$partition_prior)[1])
  }
  if (identical(ctl$assignment_prior %||% NULL, "uniform")) {
    return("fixed")
  }
  "adaptive"
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
                                            active_idx = NULL,
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
    active_idx = as.integer(active_idx %||% integer(0)),
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
  ctl <- fit$control %||% list()
  # Preserve the meaning of explicit fraction thresholds in older saved fits.
  tol <- ctl$effective_count_tol %||%
    if (!is.null(ctl$effective_weight_tol)) {
      num_features * ctl$effective_weight_tol
    } else 1e-8
  .mpcurve_validate_effective_count_tol(tol, M, num_features)
}

.mpcurve_effective_intrinsic_dim <- function(priors, M, tol, feature_counts) {
  if (is.null(priors$partition) ||
      !identical(priors$partition$mode, "adaptive")) {
    return(as.integer(M))
  }
  if (!is.numeric(feature_counts) || length(feature_counts) != M ||
      any(!is.finite(feature_counts)) || any(feature_counts < 0)) {
    stop("The fitted partition has invalid expected feature counts.", call. = FALSE)
  }
  as.integer(sum(feature_counts > tol))
}

.mpcurve_partition_assignment_info_from_fit <- function(fit) {
  ctl <- fit$control %||% list()
  partition_prior <- ctl$partition_prior %||%
    if (identical(ctl$assignment_prior %||% NULL, "uniform")) "fixed" else "adaptive"
  assignment_prior <- ctl$assignment_prior %||% NULL
  if (!is.null(assignment_prior) && !assignment_prior %in% c("uniform", "dirichlet")) {
    assignment_prior <- NULL
  }
  .cavi_partition_assignment_info(
    weights = fit$pi_weights,
    T_now = 1,
    partition_prior = partition_prior,
    partition_prior_init = ctl$partition_prior_init %||% NULL,
    assignment_prior = assignment_prior,

    assignment_M = fit$M %||% length(fit$fits),
    active_orderings = fit$active_orderings,
    active_feature_pairs = fit$active_feature_pairs,
    drop_unused_ordering = isTRUE(ctl$drop_unused_ordering)
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

.mpcurve_priors_from_legacy_fit <- function(fit) {
  list(
    position = .mpcurve_single_fit_position_priors(fit, default_mode = "estimated"),
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
      mode = if (identical(partition_mode, "legacy_dirichlet")) {
        "legacy_dirichlet"
      } else {
        partition_mode
      },
      omega = assignment_info$omega,
      init = .mpcurve_partition_prior_init(fit$control %||% list(), M = M),
      active_idx = assignment_info$active_idx,
      ordering_labels = ordering_labels
    )
  )
}

.mpcurve_get_priors <- function(x) {
  if (inherits(x, "mpcurve")) {
    if (!is.null(x$priors)) return(x$priors)
    if (.mpcurve_is_partition(x)) {
      return(.mpcurve_priors_from_partition_fit(x$fit))
    }
    if (inherits(x$fit, "cavi")) {
      return(.mpcurve_priors_from_cavi_fit(x$fit))
    }
    return(.mpcurve_priors_from_legacy_fit(x$fit))
  }
  if (inherits(x, "soft_partition_cavi")) return(.mpcurve_priors_from_partition_fit(x))
  if (inherits(x, "cavi")) return(.mpcurve_priors_from_cavi_fit(x))
  if (inherits(x, c("csmooth_em", "smooth_em"))) return(.mpcurve_priors_from_legacy_fit(x))
  NULL
}

.mpcurve_partition_param_bundle <- function(fits, ordering_labels) {
  list(
    pi = stats::setNames(lapply(fits, function(fit) fit$params$pi), ordering_labels),
    mu = stats::setNames(lapply(fits, function(fit) fit$params$mu), ordering_labels),
    sigma2 = stats::setNames(lapply(fits, function(fit) fit$params$sigma2), ordering_labels)
  )
}

.mpcurve_partition_trace_bundle <- function(fits, field, ordering_labels) {
  stats::setNames(lapply(fits, `[[`, field), ordering_labels)
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

.mpcurve_greedy_provenance_fields <- c(
  "similarity_init",
  "init_info",
  "ordering_similarity"
)

.mpcurve_check_legacy_public_args <- function(dots,
                                              caller = "fit_mpcurve()") {
  # Reject removed arguments before the ellipsis can forward them to a backend.
  removed_args <- c(
    "algorithm",
    "freeze_unused_ordering",
    "freeze_unused_ordering_threshold",
    "freeze_feature",
    "freeze_feature_weight_threshold",
    "drop_unused_ordering"
  )
  removed <- intersect(names(dots), removed_args)
  if (length(removed)) {
    stop(
      caller, " no longer accepts: ",
      paste(sprintf("`%s`", removed), collapse = ", "),
      ". Omit these arguments; fitting uses CAVI and structural VI has no ",
      "freeze/drop controls.",
      call. = FALSE
    )
  }

  legacy_args <- c(
    "relative_lambda",
    "adaptive",
    "sigma_update",
    "check_decrease",
    "tol_decrease"
  )
  bad <- intersect(names(dots), legacy_args)
  if (!length(bad)) return(invisible(NULL))

  stop(
    caller,
    " is now a CAVI-only public wrapper. Legacy smoothEM/csmoothEM controls ",
    "are no longer supported here: ",
    paste(sprintf("`%s`", bad), collapse = ", "),
    ". Use the lower-level legacy functions directly if you need those paths.",
    call. = FALSE
  )
}

.mpcurve_update_named_list <- function(x, updates) {
  if (!length(updates)) return(x)
  for (nm in names(updates)) x[[nm]] <- updates[[nm]]
  x
}

.mpcurve_last_objective <- function(x) {
  obj <- if (inherits(x, "mpcurve")) {
    if (!is.null(x$partition) || inherits(x$fit, "soft_partition_cavi")) {
      x$objective_history %||% numeric(0)
    } else {
      x$elbo_trace %||% numeric(0)
    }
  } else if (inherits(x, "soft_partition_cavi")) {
    x$objective_history %||% numeric(0)
  } else if (inherits(x, "cavi")) {
    x$elbo_trace %||% numeric(0)
  } else {
    numeric(0)
  }
  if (!length(obj)) -Inf else as.numeric(tail(obj, 1L))
}

.mpcurve_greedy_method_for_dim <- function(method,
                                           target_dim,
                                           method_missing = FALSE) {
  target_dim <- as.integer(target_dim)[1]
  if (target_dim <= 1L) {
    if (isTRUE(method_missing)) return("PCA")
    return(method[[1L]])
  }
  if (isTRUE(method_missing)) return("PCA")
  if (length(method) == 1L) return(method)
  if (length(method) >= target_dim) return(method[seq_len(target_dim)])
  method
}

.mpcurve_fit_once <- function(fit_args,
                              dots = list(),
                              fit_overrides = list(),
                              dot_overrides = list()) {
  fit_args_use <- .mpcurve_update_named_list(fit_args, fit_overrides)
  dots_use <- .mpcurve_update_named_list(dots, dot_overrides)
  do.call(
    fit_mpcurve,
    c(fit_args_use, dots_use, list(greedy = "none"))
  )
}

# Stored ordering slots determine array sizes, independently of the reported
# intrinsic dimension. Fall back to the old schema for saved pre-0.4 objects.
.mpcurve_model_intrinsic_dim <- function(x) {
  if (inherits(x, "soft_partition_cavi")) {
    return(as.integer(x$M %||% length(x$fits)))
  }
  as.integer(x$model_intrinsic_dim %||% x$num_orderings %||% x$fit$M %||%
               x$requested_intrinsic_dim %||% x$intrinsic_dim %||% 1L)
}

.mpcurve_reported_intrinsic_dim <- function(x) {
  as.integer(x$effective_intrinsic_dim %||% x$intrinsic_dim %||% 1L)
}

.mpcurve_requested_intrinsic_dim <- function(x) {
  .mpcurve_model_intrinsic_dim(x)
}

.mpcurve_active_intrinsic_dim <- function(x) {
  if (inherits(x, "mpcurve") && !is.null(x$active_intrinsic_dim)) {
    return(as.integer(x$active_intrinsic_dim))
  }
  if (inherits(x, "mpcurve") && !is.null(x$fit)) {
    x <- x$fit
  }
  if (inherits(x, "soft_partition_cavi")) {
    M <- x$M %||% length(x$fits)
    return(as.integer(sum(x$active_orderings %||% rep(TRUE, M))))
  }
  1L
}

.mpcurve_extract_greedy_provenance <- function(x) {
  out <- setNames(vector("list", length(.mpcurve_greedy_provenance_fields)),
                  .mpcurve_greedy_provenance_fields)
  sources <- list(x)
  if (inherits(x, "mpcurve") && !is.null(x$fit)) {
    sources[[length(sources) + 1L]] <- x$fit
  }

  for (source in sources) {
    source_items <- list(source)
    source_provenance <- source$greedy_provenance %||% NULL
    if (!is.null(source_provenance)) {
      source_items[[length(source_items) + 1L]] <- source_provenance
    }
    for (item in source_items) {
      for (nm in .mpcurve_greedy_provenance_fields) {
        if (is.null(out[[nm]]) && !is.null(item[[nm]])) {
          out[[nm]] <- item[[nm]]
        }
      }
    }
  }

  out[!vapply(out, is.null, logical(1))]
}

.mpcurve_merge_greedy_provenance <- function(existing = NULL, x = NULL) {
  out <- existing %||% list()
  if (is.null(x)) return(out)

  incoming <- .mpcurve_extract_greedy_provenance(x)
  for (nm in names(incoming)) {
    if (is.null(out[[nm]]) && !is.null(incoming[[nm]])) {
      out[[nm]] <- incoming[[nm]]
    }
  }
  out
}

.mpcurve_apply_greedy_provenance_raw <- function(raw_fit, provenance = NULL) {
  if (is.null(raw_fit) || !length(provenance %||% list())) {
    return(raw_fit)
  }
  for (nm in names(provenance)) {
    raw_fit[[nm]] <- provenance[[nm]]
  }
  raw_fit$greedy_provenance <- provenance
  raw_fit
}

.mpcurve_apply_greedy_provenance <- function(x, provenance = NULL) {
  if (!length(provenance %||% list())) {
    return(x)
  }
  for (nm in names(provenance)) {
    x[[nm]] <- provenance[[nm]]
  }
  x$greedy_provenance <- provenance
  if (inherits(x, "mpcurve") && !is.null(x$fit)) {
    x$fit <- .mpcurve_apply_greedy_provenance_raw(x$fit, provenance)
  }
  x
}

.mpcurve_partition_active_idx <- function(raw_fit) {
  M <- raw_fit$M %||% length(raw_fit$fits)
  which(raw_fit$active_orderings %||% rep(TRUE, M))
}

.mpcurve_partition_ordering_contributions <- function(raw_fit) {
  if (!inherits(raw_fit, "soft_partition_cavi")) {
    stop("raw_fit must inherit from 'soft_partition_cavi'.")
  }
  X <- raw_fit$fits[[1]]$data
  ctl <- raw_fit$control %||% list()
  partition_prior <- ctl$partition_prior %||% if (identical(ctl$assignment_prior %||% NULL, "uniform")) "fixed" else "adaptive"
  partition_prior_init <- ctl$partition_prior_init %||% NULL
  assignment_prior <- ctl$assignment_prior %||% NULL
  if (!is.null(assignment_prior) && !assignment_prior %in% c("uniform", "dirichlet")) {
    assignment_prior <- NULL
  }
  ordering_alpha <- ctl$ordering_alpha %||% NULL
  obj_terms <- .cavi_partition_objective_from_fits(
    fits = raw_fit$fits,
    X = X,
    weights = raw_fit$pi_weights,
    T_now = 1,
    active_orderings = raw_fit$active_orderings,
    active_feature_pairs = raw_fit$active_feature_pairs,
    partition_prior = partition_prior,
    partition_prior_init = partition_prior_init,
    assignment_prior = assignment_prior,
    ordering_alpha = ordering_alpha,
    drop_unused_ordering = FALSE
  )
  contrib <- colSums(raw_fit$pi_weights * obj_terms$like_mat) +
    colSums(obj_terms$prior_entropy_mat) +
    obj_terms$cell_terms_by_fit
  ord_labels <- colnames(raw_fit$pi_weights) %||%
    .mpcurve_ordering_labels(length(contrib))
  stats::setNames(as.numeric(contrib), ord_labels)
}

.mpcurve_partition_gamma_correlation <- function(raw_fit) {
  active_idx <- .mpcurve_partition_active_idx(raw_fit)
  if (length(active_idx) < 2L) {
    return(list(
      active_idx = active_idx,
      corr_mat = matrix(NA_real_, nrow = length(active_idx), ncol = length(active_idx))
    ))
  }
  gamma_vecs <- lapply(raw_fit$fits[active_idx], function(fit) as.numeric(fit$gamma))
  corr_mat <- matrix(NA_real_, nrow = length(active_idx), ncol = length(active_idx))
  diag(corr_mat) <- 1
  for (i in seq_along(active_idx)) {
    for (j in seq_len(i - 1L)) {
      cij <- suppressWarnings(stats::cor(gamma_vecs[[i]], gamma_vecs[[j]]))
      if (!is.finite(cij)) cij <- 0
      corr_mat[i, j] <- abs(cij)
      corr_mat[j, i] <- abs(cij)
    }
  }
  list(active_idx = active_idx, corr_mat = corr_mat)
}

.mpcurve_partition_backward_drop <- function(raw_fit) {
  corr_info <- .mpcurve_partition_gamma_correlation(raw_fit)
  active_idx <- corr_info$active_idx
  if (length(active_idx) < 2L) return(NULL)

  corr_mat <- corr_info$corr_mat
  diag(corr_mat) <- -Inf
  max_corr <- max(corr_mat, na.rm = TRUE)
  if (!is.finite(max_corr)) return(NULL)

  pair_rows <- which(corr_mat == max_corr, arr.ind = TRUE)
  pair_rows <- pair_rows[pair_rows[, 1] < pair_rows[, 2], , drop = FALSE]
  if (!nrow(pair_rows)) return(NULL)

  contrib <- .mpcurve_partition_ordering_contributions(raw_fit)
  ord_labels <- names(contrib)
  pair_tbl <- data.frame(
    i = integer(0),
    j = integer(0),
    drop_idx = integer(0),
    drop_contrib = numeric(0),
    stringsAsFactors = FALSE
  )
  for (r in seq_len(nrow(pair_rows))) {
    i_loc <- pair_rows[r, 1]
    j_loc <- pair_rows[r, 2]
    i_idx <- active_idx[i_loc]
    j_idx <- active_idx[j_loc]
    ci <- contrib[[ord_labels[i_idx]]]
    cj <- contrib[[ord_labels[j_idx]]]
    drop_idx <- if (ci <= cj) i_idx else j_idx
    drop_contrib <- min(ci, cj)
    pair_tbl <- rbind(
      pair_tbl,
      data.frame(
        i = i_idx,
        j = j_idx,
        drop_idx = drop_idx,
        drop_contrib = drop_contrib,
        stringsAsFactors = FALSE
      )
    )
  }
  best <- pair_tbl[order(pair_tbl$drop_contrib, pair_tbl$drop_idx), , drop = FALSE][1L, , drop = FALSE]
  keep_idx <- setdiff(active_idx, best$drop_idx)
  list(
    drop_idx = as.integer(best$drop_idx),
    drop_label = ord_labels[best$drop_idx],
    pair_idx = c(as.integer(best$i), as.integer(best$j)),
    pair_label = paste(ord_labels[c(best$i, best$j)], collapse = " vs "),
    corr = as.numeric(max_corr),
    keep_idx = keep_idx,
    contributions = contrib
  )
}

.mpcurve_compact_partition_fit <- function(current_fit,
                                           fit_args,
                                           dots,
                                           method_missing = FALSE,
                                           drop_unused_ordering = FALSE) {
  raw_fit <- current_fit$fit
  active_idx <- .mpcurve_partition_active_idx(raw_fit)
  active_dim <- length(active_idx)
  current_requested <- current_fit$requested_intrinsic_dim %||% current_fit$intrinsic_dim

  if (!inherits(raw_fit, "soft_partition_cavi") ||
      active_dim == current_requested ||
      active_dim < 1L) {
    return(list(fit = current_fit, compacted = FALSE, active_dim = active_dim))
  }

  if (active_dim == 1L) {
    compact_fit <- as_mpcurve(raw_fit$fits[[active_idx]])
    return(list(fit = compact_fit, compacted = TRUE, active_dim = 1L))
  }

  compact_fit <- .mpcurve_fit_once(
    fit_args = fit_args,
    dots = dots,
    fit_overrides = list(
      intrinsic_dim = active_dim,
      method = .mpcurve_greedy_method_for_dim(
        fit_args$method,
        active_dim,
        method_missing = method_missing
      ),
      drop_unused_ordering = isTRUE(drop_unused_ordering)
    ),
    dot_overrides = list(
      fits_init = raw_fit$fits[active_idx]
    )
  )
  list(fit = compact_fit, compacted = TRUE, active_dim = active_dim)
}

.mpcurve_prepare_greedy_fit <- function(fit,
                                        fit_args,
                                        dots,
                                        method_missing = FALSE,
                                        greedy_provenance = NULL,
                                        drop_unused_ordering = FALSE) {
  requested_dim <- .mpcurve_requested_intrinsic_dim(fit)
  greedy_provenance <- .mpcurve_merge_greedy_provenance(greedy_provenance, fit)
  compact_info <- .mpcurve_compact_partition_fit(
    current_fit = fit,
    fit_args = fit_args,
    dots = dots,
    method_missing = method_missing,
    drop_unused_ordering = drop_unused_ordering
  )
  normalized_fit <- .mpcurve_apply_greedy_provenance(
    compact_info$fit,
    greedy_provenance
  )
  list(
    fit = normalized_fit,
    requested_dim = as.integer(requested_dim),
    active_dim = .mpcurve_active_intrinsic_dim(normalized_fit),
    compacted = isTRUE(compact_info$compacted),
    greedy_provenance = greedy_provenance
  )
}

.mpcurve_empty_greedy_history <- function() {
  data.frame(
    step = integer(0),
    direction = character(0),
    current_M = integer(0),
    candidate_M = integer(0),
    current_requested_M = integer(0),
    current_active_M = integer(0),
    candidate_requested_M = integer(0),
    candidate_active_M = integer(0),
    current_objective = numeric(0),
    candidate_objective = numeric(0),
    accepted = logical(0),
    current_compacted = logical(0),
    candidate_compacted = logical(0),
    candidate_collapsed = logical(0),
    dropped_label = character(0),
    pair_label = character(0),
    comparison_basis = character(0),
    note = character(0),
    stringsAsFactors = FALSE
  )
}

.mpcurve_greedy_compare_candidate <- function(direction,
                                              step,
                                              current_info,
                                              candidate_info,
                                              current_objective,
                                              candidate_objective,
                                              dropped_label = NA_character_,
                                              pair_label = NA_character_) {
  direction <- match.arg(direction, c("forward", "backward"))
  candidate_collapsed <- candidate_info$active_dim < candidate_info$requested_dim

  improves_dimension <- switch(
    direction,
    forward = candidate_info$active_dim > current_info$active_dim,
    backward = candidate_info$active_dim < current_info$active_dim
  )
  objective_improves <- is.finite(candidate_objective) &&
    candidate_objective > current_objective
  accepted <- improves_dimension && objective_improves

  note <- switch(
    direction,
    forward = if (accepted) {
      if (candidate_collapsed) "larger_model_preferred_after_compaction" else "larger_model_preferred"
    } else if (!improves_dimension) {
      "candidate_not_larger_after_compaction"
    } else {
      "smaller_model_preferred"
    },
    backward = if (accepted) {
      if (candidate_collapsed) "smaller_model_preferred_after_compaction" else "smaller_model_preferred"
    } else if (!improves_dimension) {
      "candidate_not_smaller_after_compaction"
    } else {
      "larger_model_preferred"
    }
  )

  stop_reason <- switch(
    direction,
    forward = if (accepted) "" else if (!improves_dimension) {
      "candidate_not_larger_after_compaction"
    } else {
      "smaller_model_preferred"
    },
    backward = if (accepted) "" else if (!improves_dimension) {
      "candidate_not_smaller_after_compaction"
    } else {
      "larger_model_preferred"
    }
  )

  history_row <- data.frame(
    step = as.integer(step),
    direction = as.character(direction),
    current_M = as.integer(current_info$active_dim),
    candidate_M = as.integer(candidate_info$active_dim),
    current_requested_M = as.integer(current_info$requested_dim),
    current_active_M = as.integer(current_info$active_dim),
    candidate_requested_M = as.integer(candidate_info$requested_dim),
    candidate_active_M = as.integer(candidate_info$active_dim),
    current_objective = as.numeric(current_objective),
    candidate_objective = as.numeric(candidate_objective),
    accepted = isTRUE(accepted),
    current_compacted = isTRUE(current_info$compacted),
    candidate_compacted = isTRUE(candidate_info$compacted),
    candidate_collapsed = isTRUE(candidate_collapsed),
    dropped_label = as.character(dropped_label),
    pair_label = as.character(pair_label),
    comparison_basis = "active_intrinsic_dim",
    note = note,
    stringsAsFactors = FALSE
  )

  list(
    accepted = accepted,
    note = note,
    stop_reason = stop_reason,
    candidate_collapsed = candidate_collapsed,
    history_row = history_row
  )
}

.mpcurve_finalize_greedy_return <- function(selected_fit,
                                            fit_args,
                                            dots,
                                            method_missing = FALSE,
                                            greedy_info = NULL,
                                            greedy_provenance = NULL) {
  want_drop <- isTRUE(fit_args$drop_unused_ordering)
  final_fit <- selected_fit

  if (want_drop && inherits(selected_fit$fit, "soft_partition_cavi")) {
    selected_requested <- selected_fit$requested_intrinsic_dim %||% selected_fit$intrinsic_dim
    final_fit <- .mpcurve_fit_once(
      fit_args = fit_args,
      dots = dots,
      fit_overrides = list(
        intrinsic_dim = selected_requested,
        method = .mpcurve_greedy_method_for_dim(
          fit_args$method,
          selected_requested,
          method_missing = method_missing
        ),
        drop_unused_ordering = TRUE
      ),
      dot_overrides = list(
        fits_init = selected_fit$fit$fits
      )
    )
  }

  final_compact <- .mpcurve_prepare_greedy_fit(
    fit = final_fit,
    fit_args = fit_args,
    dots = dots,
    method_missing = method_missing,
    greedy_provenance = greedy_provenance,
    drop_unused_ordering = want_drop
  )
  final_fit <- final_compact$fit
  if (!is.null(greedy_info)) {
    greedy_info$selected_intrinsic_dim <- as.integer(final_compact$active_dim)
  }
  final_fit$requested_intrinsic_dim <- as.integer(final_fit$intrinsic_dim %||% final_compact$active_dim)
  final_fit$greedy_selection <- greedy_info
  final_fit
}

.mpcurve_greedy_forward <- function(fit_args, dots, method_missing = FALSE) {
  upper_bound <- as.integer(fit_args$intrinsic_dim)[1]
  history <- .mpcurve_empty_greedy_history()

  current_fit <- .mpcurve_fit_once(
    fit_args = fit_args,
    dots = dots,
    fit_overrides = list(
      intrinsic_dim = 1L,
      method = .mpcurve_greedy_method_for_dim(
        fit_args$method,
        1L,
        method_missing = method_missing
      ),
      num_cores = 1L,
      drop_unused_ordering = FALSE
    )
  )
  greedy_provenance <- .mpcurve_extract_greedy_provenance(current_fit)
  current_info <- .mpcurve_prepare_greedy_fit(
    fit = current_fit,
    fit_args = fit_args,
    dots = dots,
    method_missing = method_missing,
    greedy_provenance = greedy_provenance,
    drop_unused_ordering = FALSE
  )
  greedy_provenance <- current_info$greedy_provenance
  current_fit <- current_info$fit
  current_obj <- .mpcurve_last_objective(current_fit)
  current_M <- current_info$active_dim
  stop_reason <- if (upper_bound <= 1L) "reached_upper_bound" else "smaller_model_preferred"

  if (upper_bound >= 2L) {
    step_id <- 0L
    while (current_M < upper_bound) {
      candidate_requested_M <- current_M + 1L
      candidate_fit <- .mpcurve_fit_once(
        fit_args = fit_args,
        dots = dots,
        fit_overrides = list(
          intrinsic_dim = candidate_requested_M,
          method = .mpcurve_greedy_method_for_dim(
            fit_args$method,
            candidate_requested_M,
            method_missing = method_missing
          ),
          drop_unused_ordering = FALSE
        )
      )
      candidate_info <- .mpcurve_prepare_greedy_fit(
        fit = candidate_fit,
        fit_args = fit_args,
        dots = dots,
        method_missing = method_missing,
        greedy_provenance = greedy_provenance,
        drop_unused_ordering = FALSE
      )
      greedy_provenance <- candidate_info$greedy_provenance
      candidate_fit <- candidate_info$fit
      candidate_obj <- .mpcurve_last_objective(candidate_fit)
      step_id <- step_id + 1L
      decision <- .mpcurve_greedy_compare_candidate(
        direction = "forward",
        step = step_id,
        current_info = current_info,
        candidate_info = candidate_info,
        current_objective = current_obj,
        candidate_objective = candidate_obj
      )
      history <- rbind(history, decision$history_row)
      if (!decision$accepted) {
        stop_reason <- decision$stop_reason
        break
      }
      current_fit <- candidate_fit
      current_info <- candidate_info
      current_obj <- candidate_obj
      current_M <- current_info$active_dim
      stop_reason <- if (current_M >= upper_bound) "reached_upper_bound" else "smaller_model_preferred"
    }
  }

  greedy_info <- list(
    mode = "forward",
    requested_upper_bound = upper_bound,
    selected_intrinsic_dim = current_M,
    stop_reason = stop_reason,
    history = history
  )
  .mpcurve_finalize_greedy_return(
    selected_fit = current_fit,
    fit_args = fit_args,
    dots = dots,
    method_missing = method_missing,
    greedy_info = greedy_info,
    greedy_provenance = greedy_provenance
  )
}

.mpcurve_greedy_backward <- function(fit_args, dots, method_missing = FALSE) {
  upper_bound <- as.integer(fit_args$intrinsic_dim)[1]
  if (upper_bound <= 1L) {
    fit1 <- .mpcurve_fit_once(
      fit_args = fit_args,
      dots = dots,
      fit_overrides = list(
        intrinsic_dim = 1L,
        method = .mpcurve_greedy_method_for_dim(
          fit_args$method,
          1L,
          method_missing = method_missing
        ),
        num_cores = 1L,
        drop_unused_ordering = FALSE
      )
    )
    greedy_info <- list(
      mode = "backward",
      requested_upper_bound = upper_bound,
      selected_intrinsic_dim = 1L,
      starting_active_intrinsic_dim = 1L,
      stop_reason = "reached_dimension_1",
      history = .mpcurve_empty_greedy_history()
    )
    return(.mpcurve_finalize_greedy_return(
      selected_fit = fit1,
      fit_args = fit_args,
      dots = dots,
      method_missing = method_missing,
      greedy_info = greedy_info,
      greedy_provenance = .mpcurve_extract_greedy_provenance(fit1)
    ))
  }

  history <- .mpcurve_empty_greedy_history()

  current_fit <- .mpcurve_fit_once(
    fit_args = fit_args,
    dots = dots,
    fit_overrides = list(
      intrinsic_dim = upper_bound,
      method = .mpcurve_greedy_method_for_dim(
        fit_args$method,
        upper_bound,
        method_missing = method_missing
      ),
      drop_unused_ordering = FALSE
    )
  )
  greedy_provenance <- .mpcurve_extract_greedy_provenance(current_fit)
  start_active_dim <- if (inherits(current_fit$fit, "soft_partition_cavi")) {
    length(.mpcurve_partition_active_idx(current_fit$fit))
  } else {
    1L
  }
  current_info <- .mpcurve_prepare_greedy_fit(
    fit = current_fit,
    fit_args = fit_args,
    dots = dots,
    method_missing = method_missing,
    greedy_provenance = greedy_provenance,
    drop_unused_ordering = FALSE
  )
  greedy_provenance <- current_info$greedy_provenance
  current_fit <- current_info$fit
  current_obj <- .mpcurve_last_objective(current_fit)
  current_M <- current_info$active_dim
  stop_reason <- if (current_M <= 1L) "reached_dimension_1" else "larger_model_preferred"

  step_id <- 0L
  while (current_M > 1L && inherits(current_fit$fit, "soft_partition_cavi")) {
    current_info <- .mpcurve_prepare_greedy_fit(
      fit = current_fit,
      fit_args = fit_args,
      dots = dots,
      method_missing = method_missing,
      greedy_provenance = greedy_provenance,
      drop_unused_ordering = FALSE
    )
    greedy_provenance <- current_info$greedy_provenance
    current_fit <- current_info$fit
    current_obj <- .mpcurve_last_objective(current_fit)
    current_M <- current_info$active_dim
    if (!inherits(current_fit$fit, "soft_partition_cavi") || current_M <= 1L) {
      stop_reason <- "reached_dimension_1"
      break
    }

    drop_info <- .mpcurve_partition_backward_drop(current_fit$fit)
    if (is.null(drop_info)) {
      stop_reason <- "reached_dimension_1"
      break
    }

    candidate_requested_M <- current_M - 1L
    if (candidate_requested_M <= 1L) {
      candidate_fit <- .mpcurve_fit_once(
        fit_args = fit_args,
        dots = dots,
        fit_overrides = list(
          intrinsic_dim = 1L,
          method = .mpcurve_greedy_method_for_dim(
            fit_args$method,
            1L,
            method_missing = method_missing
          ),
          num_cores = 1L,
          drop_unused_ordering = FALSE
        ),
        dot_overrides = list(
          responsibilities_init = current_fit$fit$fits[[drop_info$keep_idx]]$gamma
        )
      )
    } else {
      candidate_fit <- .mpcurve_fit_once(
        fit_args = fit_args,
        dots = dots,
        fit_overrides = list(
          intrinsic_dim = candidate_requested_M,
          method = .mpcurve_greedy_method_for_dim(
            fit_args$method,
            candidate_requested_M,
            method_missing = method_missing
          ),
          drop_unused_ordering = FALSE
        ),
        dot_overrides = list(
          fits_init = current_fit$fit$fits[drop_info$keep_idx]
        )
      )
    }
    candidate_info <- .mpcurve_prepare_greedy_fit(
      fit = candidate_fit,
      fit_args = fit_args,
      dots = dots,
      method_missing = method_missing,
      greedy_provenance = greedy_provenance,
      drop_unused_ordering = FALSE
    )
    greedy_provenance <- candidate_info$greedy_provenance
    candidate_fit <- candidate_info$fit
    candidate_obj <- .mpcurve_last_objective(candidate_fit)
    step_id <- step_id + 1L
    decision <- .mpcurve_greedy_compare_candidate(
      direction = "backward",
      step = step_id,
      current_info = current_info,
      candidate_info = candidate_info,
      current_objective = current_obj,
      candidate_objective = candidate_obj,
      dropped_label = drop_info$drop_label,
      pair_label = drop_info$pair_label
    )
    history <- rbind(history, decision$history_row)
    if (!decision$accepted) {
      stop_reason <- decision$stop_reason
      break
    }
    current_fit <- candidate_fit
    current_info <- candidate_info
    current_obj <- candidate_obj
    current_M <- current_info$active_dim
    stop_reason <- if (current_M <= 1L) "reached_dimension_1" else "larger_model_preferred"
  }

  greedy_info <- list(
    mode = "backward",
    requested_upper_bound = upper_bound,
    selected_intrinsic_dim = current_M,
    starting_active_intrinsic_dim = start_active_dim,
    stop_reason = stop_reason,
    history = history
  )
  .mpcurve_finalize_greedy_return(
    selected_fit = current_fit,
    fit_args = fit_args,
    dots = dots,
    method_missing = method_missing,
    greedy_info = greedy_info,
    greedy_provenance = greedy_provenance
  )
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
#' @param x A \code{cavi}, \code{smooth_em}, or \code{csmooth_em} object.
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
#' For partition fits (\code{model_intrinsic_dim >= 2}), \code{$locations} is a named
#' list with one such location object per ordering. All current fits store
#' their actual model dimension in \code{$intrinsic_dim}. Automatic adaptive
#' fits remove empty orderings before conversion. Legacy objects can retain
#' an effective count distinct from their stored model slots.
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
as_mpcurve.csmooth_em <- function(x, ...) {
  params <- x$params

  # Unified mu: d x K matrix
  mu_mat <- do.call(cbind, params$mu)

  K <- length(params$pi)
  n <- if (!is.null(x$data)) nrow(x$data) else
    if (!is.null(x$gamma)) nrow(x$gamma) else NA_integer_
  d <- nrow(mu_mat)

  modelName <- x$control$modelName %||% "homoskedastic"

  structure(
    list(
      data = x$data %||% NULL,
      params = list(
        pi     = as.numeric(params$pi),
        mu     = mu_mat,
        sigma2 = params$sigma2   # d-vec (homo) or d x K matrix (hetero)
      ),
      gamma        = x$gamma,
      measurement_sd = NULL,
      locations    = .mpcurve_locations_from_gamma(x$gamma),
      elbo_trace   = x$elbo_trace   %||% numeric(0),
      loglik_trace = x$loglik_trace %||% numeric(0),
      lambda_trace = x$lambda_trace %||% list(),
      iter         = as.integer(x$iter %||% length(x$elbo_trace %||% numeric(0))),
      K            = as.integer(K),
      n            = as.integer(n),
      d            = as.integer(d),
      algorithm      = "csmooth_em",
      modelName      = modelName,
      intrinsic_dim  = 1L,
      model_intrinsic_dim = 1L,
      requested_intrinsic_dim = 1L,
      active_intrinsic_dim = 1L,
      displayed_intrinsic_dim = 1L,
      priors         = .mpcurve_priors_from_legacy_fit(x),
      converged      = x$converged %||% NULL,
      convergence_info = x$convergence_info %||% NULL,
      control        = x$control %||% list(),
      fit            = x
    ),
    class = "mpcurve"
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
      model_intrinsic_dim = 1L,
      requested_intrinsic_dim = 1L,
      active_intrinsic_dim = 1L,
      displayed_intrinsic_dim = 1L,
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
as_mpcurve.smooth_em <- function(x, ...) {
  params <- x$params

  # Unified mu: d x K matrix
  mu_mat <- do.call(cbind, params$mu)

  K <- length(params$pi)
  n <- if (!is.null(x$data)) nrow(x$data) else
    if (!is.null(x$gamma)) nrow(x$gamma) else NA_integer_
  d <- nrow(mu_mat)

  # Extract sigma2: diagonal of each sigma[[k]] -> d x K matrix
  sigma2_mat <- if (!is.null(params$sigma) && length(params$sigma) == K) {
    vapply(params$sigma, function(S) diag(as.matrix(S)), numeric(d))
  } else {
    matrix(NA_real_, d, K)
  }

  modelName <- x$control$modelName %||% (x$modelName %||% "unknown")

  structure(
    list(
      data = x$data %||% NULL,
      params = list(
        pi     = as.numeric(params$pi),
        mu     = mu_mat,
        sigma2 = sigma2_mat   # d x K  (one column per cluster)
      ),
      gamma        = x$gamma,
      measurement_sd = NULL,
      locations    = .mpcurve_locations_from_gamma(x$gamma),
      elbo_trace   = x$elbo_trace   %||% numeric(0),
      loglik_trace = x$loglik_trace %||% numeric(0),
      lambda_trace = x$lambda_trace %||% list(),
      iter         = as.integer(x$iter %||% length(x$elbo_trace %||% numeric(0))),
      K            = as.integer(K),
      n            = as.integer(n),
      d            = as.integer(d),
      algorithm      = "smooth_em",
      modelName      = modelName,
      intrinsic_dim  = 1L,
      model_intrinsic_dim = 1L,
      requested_intrinsic_dim = 1L,
      active_intrinsic_dim = 1L,
      displayed_intrinsic_dim = 1L,
      priors         = .mpcurve_priors_from_legacy_fit(x),
      converged      = x$converged %||% NULL,
      convergence_info = x$convergence_info %||% NULL,
      control        = x$control %||% list(),
      fit            = x
    ),
    class = "mpcurve"
  )
}


#' @noRd
as_mpcurve.soft_partition_cavi <- function(x, ...) {
  if (identical(x$variational_family %||% (x$control %||% list())$variational_family,
                "structured")) {
    M <- as.integer(x$M)
    ord_labels <- x$ordering_labels %||% colnames(x$pi_weights) %||%
      .mpcurve_ordering_labels(M)
    fits_mp <- stats::setNames(lapply(x$fits, as_mpcurve), ord_labels)
    locations <- stats::setNames(lapply(fits_mp, `[[`, "locations"), ord_labels)
    priors <- .mpcurve_priors_from_partition_fit(x)
    effective_count_tol <- .mpcurve_effective_count_tol(x, M, x$d)
    return(structure(
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
        model_intrinsic_dim = M,
        requested_intrinsic_dim = M,
        active_intrinsic_dim = M,
        displayed_intrinsic_dim = M,
        intrinsic_dim = if (identical(x$control$intrinsic_dim_semantics, "model")) M else
          .mpcurve_effective_intrinsic_dim(
            priors, M, effective_count_tol, colSums(x$pi_weights)),
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
    ))
  }
  M <- x$M %||% 2L
  ord_labels_full <- colnames(x$pi_weights) %||% .mpcurve_ordering_labels(M)
  active_full <- x$active_orderings %||% rep(TRUE, M)
  frozen_full <- x$frozen_orderings %||% !active_full
  drop_view <- isTRUE((x$control %||% list())$drop_unused_ordering)
  keep_idx <- if (drop_view) which(active_full) else seq_len(M)
  if (!length(keep_idx)) {
    keep_idx <- which.max(ifelse(active_full, 1, 0))
  }
  ord_labels <- ord_labels_full[keep_idx]

  # Support unified $fits list
  raw_fits <- x$fits
  if (is.null(raw_fits)) {
    # Legacy fallback (should not happen with new code)
    raw_fits <- list(x$fit1, x$fit2)
  }

  all_fits_mp <- lapply(raw_fits, as_mpcurve)
  fits_mp <- all_fits_mp[keep_idx]
  ref <- all_fits_mp[[keep_idx[1L]]]
  visible_weights <- x$pi_weights[, keep_idx, drop = FALSE]
  colnames(visible_weights) <- ord_labels
  visible_active <- if (drop_view) rep(TRUE, length(keep_idx)) else active_full[keep_idx]
  visible_frozen <- if (drop_view) rep(FALSE, length(keep_idx)) else frozen_full[keep_idx]
  priors <- .mpcurve_priors_from_partition_fit(x)
  effective_count_tol <- .mpcurve_effective_count_tol(x, M, ref$d)

  structure(
    list(
      data = ref$data %||% x$data %||% raw_fits[[1]]$data %||% NULL,
      params = .mpcurve_partition_param_bundle(all_fits_mp, ord_labels_full),
      gamma = stats::setNames(lapply(all_fits_mp, `[[`, "gamma"), ord_labels_full),
      fits      = fits_mp,
      measurement_sd = ref$measurement_sd %||% NULL,
      locations = stats::setNames(lapply(all_fits_mp, `[[`, "locations"), ord_labels_full),
      elbo_trace = .mpcurve_partition_trace_bundle(all_fits_mp, "elbo_trace", ord_labels_full),
      loglik_trace = .mpcurve_partition_trace_bundle(all_fits_mp, "loglik_trace", ord_labels_full),
      lambda_trace = .mpcurve_partition_trace_bundle(all_fits_mp, "lambda_trace", ord_labels_full),
      partition = list(
        pi_weights = visible_weights,
        assign     = x$assign,
        active_orderings = visible_active,
        frozen_orderings = visible_frozen,
        frozen_labels = ord_labels_full[frozen_full],
        ordering_events = x$ordering_events %||% list(),
        assignment_posterior = x$assignment_posterior,
        requested_intrinsic_dim = as.integer(M),
        active_intrinsic_dim = as.integer(sum(active_full)),
        displayed_intrinsic_dim = as.integer(length(keep_idx)),
        dropped_labels = if (drop_view) ord_labels_full[frozen_full] else character(0),
        kept_labels = ord_labels,
        compacted = drop_view
      ),
      objective_history = x$objective_history %||% numeric(0),
      iter          = length(x$objective_history %||% numeric(0)),
      K             = ref$K,
      n             = ref$n,
      d             = ref$d,
      algorithm     = "cavi",
      modelName     = "partition_cavi",
      model_intrinsic_dim = as.integer(M),
      requested_intrinsic_dim = as.integer(M),
      active_intrinsic_dim = as.integer(sum(active_full)),
      displayed_intrinsic_dim = as.integer(length(keep_idx)),
      intrinsic_dim = .mpcurve_effective_intrinsic_dim(
        priors, M, effective_count_tol, colSums(x$pi_weights)
      ),
      effective_count_tol = effective_count_tol,
      priors        = priors,
      init_info     = x$init_info,
      ordering_similarity = x$ordering_similarity,
      similarity_init = x$similarity_init %||% NULL,
      converged     = x$converged,
      convergence_info = x$convergence_info,
      ordering_events = x$ordering_events %||% list(),
      control       = x$control %||% list(),
      fit           = x
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
  idim <- .mpcurve_reported_intrinsic_dim(x)
  M <- .mpcurve_model_intrinsic_dim(x)
  active_dim <- x$active_intrinsic_dim %||% .mpcurve_active_intrinsic_dim(x)
  displayed_dim <- x$displayed_intrinsic_dim %||% if (.mpcurve_is_partition(x)) length(x$fits %||% list()) else 1L
  is_partition <- .mpcurve_is_partition(x)

  if (is_partition) {
    structured <- identical(
      x$variational_family %||% (x$control %||% list())$variational_family,
      "structured"
    )
    cat("MPCurve partition fit\n")
    cat(sprintf("  Backend        : %s\n", x$algorithm))
    if (structured) cat("  Variational    : structural q(C) q(Z) q(U | Z)\n")
    selection_info <- x$dimension_selection %||% NULL
    if (!is.null(selection_info)) {
      cat(sprintf(
        "  M selection   : %s (upper bound %d -> selected %d)\n",
        selection_info$direction,
        (selection_info$max_intrinsic_dim %||% selection_info$max_num_orderings),
        selection_info$selected_M
      ))
    }
    initialization_info <- x$dimension_initialization %||% NULL
    if (!is.null(initialization_info) && is.null(x$dimension_estimation)) {
      cat(sprintf(
        "  Auto M init    : mean silhouette (upper bound %d -> selected %d; min size %d)\n",
        initialization_info$requested_max_intrinsic_dim,
        initialization_info$selected_M,
        initialization_info$min_cluster_size
      ))
    }
    greedy_info <- x$greedy_selection %||% NULL
    if (!is.null(greedy_info)) {
      cat(sprintf(
        "  Greedy search  : %s (upper bound %d -> selected %d)\n",
        greedy_info$mode,
        greedy_info$requested_upper_bound,
        greedy_info$selected_intrinsic_dim
      ))
    }
    if (structured) {
      if (M != idim) cat(sprintf("  Model dim      : %d\n", M))
    } else {
      cat(sprintf("  Dim (req/act/view): %d / %d / %d\n", M, active_dim, displayed_dim))
    }
    cat(sprintf("  Intrinsic dim  : %d\n", idim))
    cat(sprintf("  n / d / K      : %d / %d / %d\n", x$n, x$d, x$K))
    cat(sprintf("  Iterations     : %d\n", x$iter))
    part <- x$partition
    if (!is.null(part)) {
      tbl <- table(part$assign)
      part_str <- paste(sprintf("%s=%d", names(tbl), as.integer(tbl)), collapse = ", ")
      cat(sprintf("  Partition      : %s\n", part_str))
      if (!structured) {
        active <- part$active_orderings
        ord_labels <- colnames(part$pi_weights) %||% .mpcurve_ordering_labels(M)
        if (!is.null(active)) {
          cat(sprintf("  Active         : %s\n",
                      paste(ord_labels[active], collapse = ", ")))
        }
        frozen <- part$frozen_labels
        if (!is.null(frozen) && length(frozen) > 0L) {
          if (isTRUE(part$compacted)) {
            cat(sprintf("  Dropped from view: %s\n",
                        paste(frozen, collapse = ", ")))
          } else {
            cat(sprintf("  Frozen         : %s\n",
                        paste(frozen, collapse = ", ")))
          }
        }
      }
    }
    if (length(x$objective_history) > 0L)
      cat(sprintf("  Objective (last): %.6f\n",
                  tail(x$objective_history, 1L)))
    if (!is.null(x$converged)) {
      cat(sprintf("  Converged      : %s\n", x$converged))
    }
  } else {
    cat("MPCurve fit\n")
    cat(sprintf("  Backend        : %s\n",  x$algorithm))
    selection_info <- x$dimension_selection %||% NULL
    if (!is.null(selection_info)) {
      cat(sprintf(
        "  M selection   : %s (upper bound %d -> selected %d)\n",
        selection_info$direction,
        (selection_info$max_intrinsic_dim %||% selection_info$max_num_orderings),
        selection_info$selected_M
      ))
    }
    initialization_info <- x$dimension_initialization %||% NULL
    if (!is.null(initialization_info) && is.null(x$dimension_estimation)) {
      cat(sprintf(
        "  Auto M init    : mean silhouette (upper bound %d -> selected %d; min size %d)\n",
        initialization_info$requested_max_intrinsic_dim,
        initialization_info$selected_M,
        initialization_info$min_cluster_size
      ))
    }
    greedy_info <- x$greedy_selection %||% NULL
    if (!is.null(greedy_info)) {
      cat(sprintf(
        "  Greedy search  : %s (upper bound %d -> selected %d)\n",
        greedy_info$mode,
        greedy_info$requested_upper_bound,
        greedy_info$selected_intrinsic_dim
      ))
    }
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
#'   \code{$model_intrinsic_dim},
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
  idim <- .mpcurve_reported_intrinsic_dim(object)
  M <- .mpcurve_model_intrinsic_dim(object)
  active_dim <- object$active_intrinsic_dim %||% .mpcurve_active_intrinsic_dim(object)
  displayed_dim <- object$displayed_intrinsic_dim %||% if (.mpcurve_is_partition(object)) length(object$fits %||% list()) else 1L
  priors <- .mpcurve_get_priors(object)
  is_partition <- .mpcurve_is_partition(object)

  if (is_partition) {
    result <- list(
      algorithm     = object$algorithm,
      intrinsic_dim = idim,
      model_intrinsic_dim = M,
      requested_intrinsic_dim = object$requested_intrinsic_dim %||% M,
      active_intrinsic_dim = active_dim,
      displayed_intrinsic_dim = displayed_dim,
      effective_count_tol = object$effective_count_tol %||%
        .mpcurve_effective_count_tol(object$fit, M, object$d),
      modelName     = object$modelName %||% "partition_cavi",
      K             = object$K,
      n             = object$n,
      d             = object$d,
      priors        = priors,
      partition     = object$partition,
      objective_history = object$objective_history,
      converged     = object$converged,
      convergence_info = object$convergence_info %||% object$fit$convergence_info,
      ordering_events = object$ordering_events %||% list(),
      control       = object$control %||% list(),
      sigma2        = object$params$sigma2,
      measurement_sd = object$measurement_sd,
      variational_family = object$variational_family %||%
        (object$control %||% list())$variational_family,
      greedy_selection = object$greedy_selection %||% NULL,
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
      model_intrinsic_dim = M,
      active_intrinsic_dim = active_dim,
      displayed_intrinsic_dim = displayed_dim,
      K             = object$K,
      n             = object$n,
      d             = object$d,
      priors        = priors,
      converged     = object$converged %||% underlying$converged %||% NULL,
      underlying    = underlying,
      greedy_selection = object$greedy_selection %||% NULL,
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
  idim <- .mpcurve_reported_intrinsic_dim(x)
  M <- .mpcurve_model_intrinsic_dim(x)
  active_dim <- x$active_intrinsic_dim %||% M
  displayed_dim <- x$displayed_intrinsic_dim %||% if (!is.null(x$partition)) length(x$partition$kept_labels %||% character(0)) else 1L
  is_partition <- !is.null(x$partition)

  if (is_partition) {
    structured <- identical(
      x$variational_family %||% (x$control %||% list())$variational_family,
      "structured"
    )
    cat("MPCurve Partition Summary\n")
    selection_info <- x$dimension_selection %||% NULL
    if (!is.null(selection_info)) {
      cat(sprintf("M selection : %s  |  upper bound = %d  |  selected = %d\n",
                  selection_info$direction,
                  (selection_info$max_intrinsic_dim %||% selection_info$max_num_orderings),
                  selection_info$selected_M))
    }
    initialization_info <- x$dimension_initialization %||% NULL
    if (!is.null(initialization_info) && is.null(x$dimension_estimation)) {
      cat(sprintf(
        "Auto M init : mean silhouette  |  upper bound = %d  |  selected = %d  |  min size = %d\n",
        initialization_info$requested_max_intrinsic_dim,
        initialization_info$selected_M,
        initialization_info$min_cluster_size
      ))
    }
    greedy_info <- x$greedy_selection %||% NULL
    if (!is.null(greedy_info)) {
      cat(sprintf("Greedy search : %s  |  upper bound = %d  |  selected = %d\n",
                  greedy_info$mode,
                  greedy_info$requested_upper_bound,
                  greedy_info$selected_intrinsic_dim))
    }
    if (structured) {
      cat(sprintf("Algorithm : %s  |  structural VI  |  n=%d  d=%d  K=%d\n",
                  x$algorithm, x$n, x$d, x$K))
      if (M != idim) cat(sprintf("Retained model dimension : %d\n", M))
      cat(sprintf("Intrinsic dimension : %d\n", idim))
      if (is.null(x$measurement_sd)) {
        cat(sprintf("Shared sigma2 range : [%.4g, %.4g]\n", min(x$sigma2), max(x$sigma2)))
      } else {
        cat("Noise model : known measurement SD\n")
      }
    } else {
      cat(sprintf("Algorithm : %s  |  requested=%d  active=%d  displayed=%d  |  n=%d  d=%d  K=%d\n",
                  x$algorithm, M, active_dim, displayed_dim, x$n, x$d, x$K))
      cat(sprintf("Intrinsic dimension : %d\n", idim))
    }
    if (!is.null(x$partition)) {
      tbl <- table(x$partition$assign)
      part_str <- paste(sprintf("%s=%d", names(tbl), as.integer(tbl)), collapse = "  ")
      cat(sprintf("Partition : %s\n", part_str))
      if (!structured) {
        active <- x$partition$active_orderings
        ord_labels <- colnames(x$partition$pi_weights) %||% .mpcurve_ordering_labels(M)
        if (!is.null(active)) {
          cat(sprintf("Active    : %s\n", paste(ord_labels[active], collapse = ", ")))
        }
        frozen <- x$partition$frozen_labels
        if (!is.null(frozen) && length(frozen) > 0L) {
          if (isTRUE(x$partition$compacted)) {
            cat(sprintf("Dropped from view : %s\n", paste(frozen, collapse = ", ")))
          } else {
            cat(sprintf("Frozen    : %s\n", paste(frozen, collapse = ", ")))
          }
        }
      }
    }
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
    events <- x$ordering_events %||% list()
    if (!structured && length(events) > 0L) {
      event_types <- vapply(events, `[[`, character(1), "event")
      cat(sprintf("Events    : %d freeze\n",
                  sum(event_types == "freeze")))
    }
    if (length(x$objective_history) > 0L)
      cat(sprintf("Objective (last) : %.6f\n", tail(x$objective_history, 1L)))
    if (structured) {
      cat("Objective type : fixed-M structural ELBO\n")
    } else {
      drop_view <- isTRUE((x$control %||% list())$drop_unused_ordering)
      cat(sprintf("Objective type : %s\n",
                  if (drop_view) {
                    "post-drop fit (not comparable across requested dimensions)"
                  } else {
                    "fixed-requested-M comparison"
                  }))
    }
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
    selection_info <- x$dimension_selection %||% NULL
    if (!is.null(selection_info)) {
      cat(sprintf("M selection : %s  |  upper bound = %d  |  selected = %d\n",
                  selection_info$direction,
                  (selection_info$max_intrinsic_dim %||% selection_info$max_num_orderings),
                  selection_info$selected_M))
    }
    initialization_info <- x$dimension_initialization %||% NULL
    if (!is.null(initialization_info) && is.null(x$dimension_estimation)) {
      cat(sprintf(
        "Auto M init : mean silhouette  |  upper bound = %d  |  selected = %d  |  min size = %d\n",
        initialization_info$requested_max_intrinsic_dim,
        initialization_info$selected_M,
        initialization_info$min_cluster_size
      ))
    }
    greedy_info <- x$greedy_selection %||% NULL
    if (!is.null(greedy_info)) {
      cat(sprintf("Greedy search : %s  |  upper bound = %d  |  selected = %d\n",
                  greedy_info$mode,
                  greedy_info$requested_upper_bound,
                  greedy_info$selected_intrinsic_dim))
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


#' @noRd
.plot_structural_mu_mpcurve <- function(x, dims, ...) {
  structured <- identical(
    x$variational_family %||% (x$control %||% list())$variational_family,
    "structured"
  )
  if (!structured) {
    stop(
      "plot_type = 'mu' for a partition fit requires a canonical structural state.",
      call. = FALSE
    )
  }

  means <- x$conditional_posterior$mean
  if (!is.list(means) || length(means) < 1L) {
    stop(
      "x$conditional_posterior$mean must be a non-empty list of trajectory mean matrices.",
      call. = FALSE
    )
  }

  M <- .mpcurve_model_intrinsic_dim(x)
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
  displayed_dim <- x$displayed_intrinsic_dim %||% if (is_partition) length(x$fits %||% list()) else 1L

  # ---- Multi-ordering partition plot ----
  if (is_partition && plot_type == "scatterplot") {
    M <- displayed_dim
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
    structured <- identical(
      x$variational_family %||% (x$control %||% list())$variational_family,
      "structured"
    )
    if (structured) {
      return(.plot_structural_mu_mpcurve(x, dims = dims, ...))
    }

    # Legacy partition objects remain readable through their derived fit views.
    legacy_fits <- x$fits
    if (!is.list(legacy_fits) || length(legacy_fits) < 1L) {
      stop(
        "Legacy partition trajectory plotting requires a non-empty x$fits list.",
        call. = FALSE
      )
    }
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = c(1L, length(legacy_fits)), mar = c(4, 4, 4, 1))
    for (m in seq_along(legacy_fits)) {
      .mpcurve_plot_call(plot, list(legacy_fits[[m]], plot_type = "mu", dims = dims), ...)
      graphics::mtext(
        sprintf("Ordering %s", names(legacy_fits)[m] %||% m),
        side = 3,
        line = 0.25,
        cex = 0.8
      )
    }
    return(invisible(x))
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

    .mpcurve_plot_call(plot_EM_embedding2D, list(
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
#' @param object An `mpcurve` object returned by [fit_mpcurve()]. Legacy
#'   augmented partition fits are read-only and require refitting.
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
  if (!identical(object$algorithm, "cavi")) {
    stop("Only CAVI fits can be continued; refit with fit_mpcurve().", call. = FALSE)
  }
  raw <- object$fit
  if (!inherits(raw, c("cavi", "soft_partition_cavi"))) {
    stop("object$fit must be a cavi or soft_partition_cavi object.")
  }
  partition <- inherits(raw, "soft_partition_cavi")
  if (partition && !identical(raw$variational_family %||%
                              raw$control$variational_family, "structured")) {
    stop("Legacy augmented partition fits are read-only. Refit with fit_mpcurve().",
         call. = FALSE)
  }
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
  updated$control$intrinsic_dim_semantics <- raw$control$intrinsic_dim_semantics
  updated <- .mpcurve_finalize_automatic_dimension(updated)
  out <- as_mpcurve(updated)
  out$dimension_selection <- object$dimension_selection %||% NULL
  if (!is.null(out$dimension_selection)) {
    out$dimension_selection$continued_after_selection <- TRUE
  }
  out$dimension_initialization <- updated$dimension_initialization
  out$requested_intrinsic_dim <- if (isTRUE(updated$dimension_initialization$automatic)) {
    out$model_intrinsic_dim
  } else object$requested_intrinsic_dim %||% out$model_intrinsic_dim
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
  .mpcurve_check_legacy_public_args(dots, caller = "fit_mpcurve()")
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
    raw$control$intrinsic_dim_semantics <- "model"
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
