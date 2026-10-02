# Finalize automatic dimension estimation without additional CAVI sweeps.
.mpcurve_finalize_automatic_dimension <- function(fit) {
  initialization <- fit$dimension_initialization
  if (!isTRUE(initialization$automatic)) return(fit)

  partition <- inherits(fit, "soft_partition_cavi")
  M <- if (partition) as.integer(fit$M) else 1L
  d <- ncol(fit$data)
  tol <- .mpcurve_effective_count_tol(fit, M, d)
  adaptive_partition <- partition && identical(fit$control$partition_prior, "adaptive")
  if (adaptive_partition && tol > 1e-6) {
    warning(
      "effective_count_tol exceeds 1e-6. Automatic pruning removes orderings ",
      "with estimated prior weight <= effective_count_tol / ncol(X); ",
      "large tolerances may remove supported orderings. ",
      "Use the default 1e-8 for numerical-zero cleanup.",
      call. = FALSE
    )
  }
  counts <- if (partition) colSums(fit$pi_weights) else
    stats::setNames(as.numeric(d),
      names(fit$dimension_estimation$expected_feature_counts) %||%
        .mpcurve_single_ordering_label())
  # Adaptive EB gives omega_m = counts_m / d. Compare before division to
  # preserve the existing numerical boundary at omega_m <= tol / d.
  keep <- if (adaptive_partition) {
    which(counts > tol)
  } else seq_len(M)
  if (!length(keep)) {
    stop("Automatic dimension finalization must retain at least one ordering.",
         call. = FALSE)
  }

  estimation <- fit$dimension_estimation %||% list(
    initial_intrinsic_dim = initialization$selected_M,
    pruning = list()
  )
  estimation$intrinsic_dim <- as.integer(length(keep))
  estimation$effective_count_tol <- tol
  fit$control$intrinsic_dim_semantics <- "model"
  if (length(keep) == M) {
    estimation$expected_feature_counts <- counts
    fit$dimension_estimation <- estimation
    return(fit)
  }

  state <- .structural_partition_state_from_fit(fit)
  event <- list(
    from_intrinsic_dim = M,
    to_intrinsic_dim = as.integer(length(keep)),
    retained_orderings = state$ordering_labels[keep],
    removed_orderings = state$ordering_labels[-keep],
    expected_feature_counts = counts,
    prior_weights = state$assignment_info$omega,
    objective_before = tail(fit$objective_history, 1L),
    refitted = FALSE,
    iterations = fit$iter,
    # These traces belong to the preceding dimension, not the returned model.
    fitting_history = list(
      objective = fit$objective_history,
      temperature = fit$temperature_history,
      weights = fit$weight_history,
      local_blocks = fit$score_history,
      sigma2 = fit$sigma2_trace,
      lambda = fit$lambda_trace
    )
  )
  state$M <- as.integer(length(keep))
  state$ordering_labels <- state$ordering_labels[keep]
  state$gamma <- state$gamma[keep]
  state$position_pi <- state$position_pi[keep]
  state$q_u <- state$q_u[keep]
  state$lambda_mat <- state$lambda_mat[, keep, drop = FALSE]
  state$local_blocks <- state$local_blocks[, keep, drop = FALSE]
  state$pi_weights <- state$pi_weights[, keep, drop = FALSE]
  retained_mass <- rowSums(state$pi_weights)
  if (any(!is.finite(retained_mass)) || any(retained_mass <= 0)) {
    stop(
      "Automatic dimension pruning leaves a feature with no positive ",
      "assignment probability among the retained orderings. ",
      "Lower effective_count_tol; the default 1e-8 removes numerical-zero prior weights.",
      call. = FALSE
    )
  }
  state$pi_weights <- state$pi_weights / retained_mass
  state$init_info <- state$init_info[keep]
  if (is.matrix(state$ordering_similarity) &&
      identical(dim(state$ordering_similarity), c(M, M))) {
    state$ordering_similarity <- state$ordering_similarity[keep, keep, drop = FALSE]
  }
  state$assignment_info <- .structural_partition_assignment_info(
    state$pi_weights, 1, state$control)
  state$objective_terms <- .structural_partition_objective(
    state$gamma, state$position_pi, state$pi_weights, state$local_blocks,
    state$assignment_info)
  state$objective <- state$objective_terms$objective
  state$iteration_offset <- as.integer(fit$iter)
  event$objective_after <- state$objective
  estimation$expected_feature_counts <- colSums(state$pi_weights)
  estimation$pruning <- c(estimation$pruning, list(event))
  state$dimension_estimation <- estimation

  # Start a fresh trace for the reduced model, retaining earlier traces above.
  histories <- .structural_partition_append_history(
    list(objective = numeric(), temperature = numeric(), weights = list(),
         local_blocks = list(), sigma2 = list(), lambda = list()), state, 1)
  reduced <- .structural_partition_finalize(state, histories, fit$converged)
  if (state$M == 1L) {
    # Use the ordinary single-ordering schema, preserving the fitted posterior.
    reduced <- reduced$fits[[1L]]
    reduced$elbo_trace <- histories$objective
    reduced$lambda_trace <- list(as.numeric(state$lambda_mat[, 1L]))
    reduced$sigma2_trace <- list(state$sigma2)
    reduced$iter <- state$iteration_offset
    reduced$iteration_offset <- state$iteration_offset
    reduced$converged <- fit$converged
    reduced$control$compatibility_view <- NULL
    reduced$control$variational_family <- NULL
    reduced$control$method <- state$init_info[[1L]]$method_used %||% "PCA"
    reduced$control$tol <- fit$control$tol_outer
    reduced$init_info <- state$init_info[[1L]]
    reduced$dimension_initialization <- initialization
    reduced$dimension_estimation <- estimation
  }
  reduced
}
