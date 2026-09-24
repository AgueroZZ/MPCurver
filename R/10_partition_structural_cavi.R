# Structural variational inference for fixed-M MPCurve partition models.
#
# The canonical variational family is
#   q(C) prod_j q(Z_j) q(U_j | Z_j).
# Conditional trajectory factors are indexed by feature and ordering, but they
# are mutually exclusive branches under q(Z_j), not independent generative
# trajectories. New fits and continuation use only the top-level state below;
# `$fits` is generated as a read-only compatibility view.

.structural_partition_warn_legacy_controls <- function(
    freeze_unused_ordering = NULL,
    freeze_unused_ordering_threshold = NULL,
    freeze_feature = NULL,
    freeze_feature_weight_threshold = NULL,
    drop_unused_ordering = NULL,
    caller = "soft_partition_cavi()") {
  supplied <- c(
    freeze_unused_ordering = !is.null(freeze_unused_ordering),
    freeze_unused_ordering_threshold = !is.null(freeze_unused_ordering_threshold),
    freeze_feature = !is.null(freeze_feature),
    freeze_feature_weight_threshold = !is.null(freeze_feature_weight_threshold),
    drop_unused_ordering = !is.null(drop_unused_ordering)
  )
  if (any(supplied)) {
    warning(
      caller,
      ": ",
      paste(names(supplied)[supplied], collapse = ", "),
      " are deprecated and ignored. Structural VI does not use active, ",
      "freeze, or drop states.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

.structural_partition_validate_weights <- function(weights, d, M) {
  weights <- as.matrix(weights)
  if (!all(dim(weights) == c(d, M)) || any(!is.finite(weights)) ||
      any(weights < 0)) {
    stop("pi_weights must be a finite nonnegative d x M matrix.", call. = FALSE)
  }
  rs <- rowSums(weights)
  if (any(rs <= 0)) {
    stop("Every row of pi_weights must have positive mass.", call. = FALSE)
  }
  weights / rs
}

.structural_partition_stabilize_gamma <- function(R) {
  # Exact categorical coordinate updates are strictly positive for finite
  # logits. Floating-point softmax can nevertheless underflow to exact zeros,
  # so restore a numerically representable simplex point. This is not a
  # precision ridge and is not applied to feature-assignment weights.
  R <- as.matrix(R)
  floor_prob <- sqrt(.Machine$double.eps) / ncol(R)
  R <- pmax(R, floor_prob)
  R / rowSums(R)
}

.structural_partition_lambda_matrix <- function(lambda_init, fits, d, M,
                                                lambda_min, lambda_max) {
  if (is.null(lambda_init)) {
    out <- vapply(seq_len(M), function(m) {
      value <- as.numeric(fits[[m]]$lambda_vec %||% rep(1, d))
      if (length(value) == 1L) value <- rep(value, d)
      if (length(value) != d) rep(1, d) else value
    }, numeric(d))
  } else if (is.matrix(lambda_init)) {
    out <- as.matrix(lambda_init)
    if (!all(dim(out) == c(d, M))) {
      stop("lambda_init matrix must have dimensions ncol(X) x M.", call. = FALSE)
    }
  } else {
    value <- as.numeric(lambda_init)
    if (length(value) == 1L) {
      out <- matrix(value, nrow = d, ncol = M)
    } else if (length(value) == d) {
      out <- matrix(value, nrow = d, ncol = M)
    } else if (length(value) == d * M) {
      out <- matrix(value, nrow = d, ncol = M)
    } else {
      stop(
        "lambda_init must have length 1, ncol(X), ncol(X) * M, or be a d x M matrix.",
        call. = FALSE
      )
    }
  }
  if (any(!is.finite(out)) || any(out <= 0)) {
    stop("lambda_init values must be finite and positive.", call. = FALSE)
  }
  pmin(pmax(out, lambda_min), lambda_max)
}

.structural_partition_sigma2 <- function(X, sigma2_init, measurement_sd,
                                         sigma_min, sigma_max) {
  d <- ncol(X)
  if (!is.null(measurement_sd)) {
    if (!is.null(sigma2_init)) {
      stop("sigma2_init cannot be supplied when S is non-NULL.", call. = FALSE)
    }
    return(NULL)
  }
  if (is.null(sigma2_init)) {
    sigma2 <- apply(X, 2, stats::var)
  } else {
    sigma2 <- as.numeric(sigma2_init)
    if (length(sigma2) == 1L) sigma2 <- rep(sigma2, d)
    if (length(sigma2) != d) {
      stop("sigma2_init must have length 1 or ncol(X).", call. = FALSE)
    }
  }
  if (any(!is.finite(sigma2)) || any(sigma2 < 0)) {
    stop("sigma2_init must contain finite nonnegative values.", call. = FALSE)
  }
  pmin(pmax(sigma2, sigma_min), sigma_max)
}

.structural_partition_q_u <- function(X, R, sigma2, measurement_sd,
                                      lambda_vec, Q_K,
                                      ordering_label = NULL) {
  X <- as.matrix(X)
  R <- as.matrix(R)
  n <- nrow(X)
  d <- ncol(X)
  K <- ncol(R)
  if (nrow(R) != n || length(lambda_vec) != d) {
    stop("Incompatible X, R, or lambda dimensions in structural q(U).", call. = FALSE)
  }
  noise_info <- .cavi_resolve_measurement_sd(
    S = measurement_sd,
    X = X,
    caller = ".structural_partition_q_u()"
  )
  known_noise <- !is.null(noise_info$measurement_sd)
  if (!known_noise && (length(sigma2) != d || any(!is.finite(sigma2)) || any(sigma2 <= 0))) {
    stop("sigma2 must be a finite positive vector of length ncol(X).", call. = FALSE)
  }

  Nk <- colSums(R)
  GX <- if (known_noise) NULL else t(X) %*% R
  m_mat <- matrix(0, nrow = d, ncol = K)
  sdiag_mat <- matrix(0, nrow = d, ncol = K)
  S_list <- vector("list", d)
  logdetS <- numeric(d)
  eq_quad <- numeric(d)

  for (j in seq_len(d)) {
    if (known_noise) {
      inv_vj <- 1 / noise_info$var_mat[, j]
      diag_j <- as.numeric(crossprod(inv_vj, R))
      rhs_j <- as.numeric(crossprod(inv_vj * X[, j], R))
    } else {
      diag_j <- as.numeric(Nk / sigma2[j])
      rhs_j <- as.numeric(GX[j, ] / sigma2[j])
    }
    A_j <- lambda_vec[j] * Q_K + diag(diag_j, K, K)
    A_j <- 0.5 * (A_j + t(A_j))
    U_j <- tryCatch(chol(A_j), error = function(e) NULL)
    if (is.null(U_j)) {
      label <- ordering_label %||% "unknown"
      stop(
        sprintf(
          paste0(
            "Structural q(U) precision is not positive definite for feature %d ",
            "and ordering %s. The likelihood does not identify the intrinsic-prior ",
            "null space under the current position responsibilities."
          ),
          j,
          label
        ),
        call. = FALSE
      )
    }
    S_j <- chol2inv(U_j)
    m_j <- as.numeric(backsolve(U_j, forwardsolve(t(U_j), rhs_j)))
    S_list[[j]] <- S_j
    m_mat[j, ] <- m_j
    sdiag_mat[j, ] <- diag(S_j)
    logdetS[j] <- -2 * sum(log(diag(U_j)))
    eq_quad[j] <- as.numeric(crossprod(m_j, Q_K %*% m_j) + sum(Q_K * S_j))
  }

  list(
    m_mat = m_mat,
    sdiag_mat = sdiag_mat,
    S_list = S_list,
    logdetS = logdetS,
    eq_quad = eq_quad
  )
}

.structural_partition_update_gamma <- function(X, q_u, pi_vec, sigma2,
                                               measurement_sd, weights) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  K <- length(pi_vec)
  weights <- as.numeric(weights)
  if (length(weights) != d || any(!is.finite(weights)) || any(weights < 0)) {
    stop("Structural q(C) weights must be finite, nonnegative, and length d.", call. = FALSE)
  }
  noise_info <- .cavi_resolve_measurement_sd(
    S = measurement_sd,
    X = X,
    caller = ".structural_partition_update_gamma()"
  )
  known_noise <- !is.null(noise_info$measurement_sd)
  log_r <- matrix(0, nrow = n, ncol = K)
  if (!known_noise) {
    inv_var_weight <- matrix(weights / sigma2, nrow = n, ncol = d, byrow = TRUE)
    const_term <- sum(weights * log(2 * base::pi * sigma2))
  } else {
    inv_var_weight <- matrix(weights, nrow = n, ncol = d, byrow = TRUE) /
      noise_info$var_mat
    const_term <- rowSums(
      matrix(weights, nrow = n, ncol = d, byrow = TRUE) *
        log(2 * base::pi * noise_info$var_mat)
    )
  }
  for (k in seq_len(K)) {
    diff_k <- sweep(X, 2, q_u$m_mat[, k], "-")
    quad_det <- rowSums(inv_var_weight * diff_k^2)
    quad_unc <- if (!known_noise) {
      rep(sum(weights * q_u$sdiag_mat[, k] / sigma2), n)
    } else {
      drop(inv_var_weight %*% q_u$sdiag_mat[, k])
    }
    log_r[, k] <- log(pmax(pi_vec[k], .Machine$double.eps)) -
      0.5 * (const_term + quad_det + quad_unc)
  }
  lse <- matrixStats::rowLogSumExps(log_r)
  .structural_partition_stabilize_gamma(exp(log_r - lse))
}

.structural_partition_expected_residuals <- function(X, R, q_u) {
  X <- as.matrix(X)
  R <- as.matrix(R)
  n <- nrow(X)
  d <- ncol(X)
  K <- ncol(R)
  Nk <- colSums(R)
  out <- numeric(d)
  for (j in seq_len(d)) {
    value <- 0
    for (k in seq_len(K)) {
      value <- value + sum(R[, k] * (X[, j] - q_u$m_mat[j, k])^2) +
        Nk[k] * q_u$sdiag_mat[j, k]
    }
    out[j] <- value
  }
  out
}

.structural_partition_update_sigma2 <- function(X, gamma, q_u, weights,
                                                sigma_min, sigma_max) {
  X <- as.matrix(X)
  d <- ncol(X)
  M <- length(gamma)
  residuals <- vapply(seq_len(M), function(m) {
    .structural_partition_expected_residuals(X, gamma[[m]], q_u[[m]])
  }, numeric(d))
  sigma2 <- rowSums(weights * residuals) / nrow(X)
  list(
    sigma2 = pmin(pmax(sigma2, sigma_min), sigma_max),
    residuals = residuals
  )
}

.structural_partition_local_blocks <- function(X, gamma, q_u, sigma2,
                                               measurement_sd, lambda_mat,
                                               Q_K, rw_q,
                                               lambda_sd_prior_rate) {
  d <- ncol(X)
  M <- length(gamma)
  out <- matrix(NA_real_, nrow = d, ncol = M)
  for (m in seq_len(M)) {
    data_prior <- .cavi_feature_scores_from_state(
      X = X,
      R = gamma[[m]],
      q_u = q_u[[m]],
      sigma2 = sigma2,
      measurement_sd = measurement_sd,
      lambda_vec = lambda_mat[, m],
      Q_K = Q_K,
      rw_q = rw_q,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      include_prior = TRUE
    )
    out[, m] <- data_prior + .cavi_entropy_u_terms(q_u[[m]], K = ncol(gamma[[m]]))
  }
  out
}

.structural_partition_assignment_info <- function(weights, T_now, control) {
  M <- ncol(weights)
  .cavi_partition_assignment_info(
    weights = weights,
    T_now = T_now,
    partition_prior = control$partition_prior,
    partition_prior_init = control$partition_prior_init,
    assignment_prior = control$assignment_prior,
    ordering_alpha = control$ordering_alpha,
    assignment_M = M,
    active_orderings = rep(TRUE, M),
    active_feature_pairs = matrix(TRUE, nrow(weights), M),
    drop_unused_ordering = FALSE,
    partition_prior_missing = FALSE,
    caller = ".structural_partition_assignment_info()"
  )
}

.structural_partition_update_weights <- function(local_blocks, weights, T_now,
                                                 control) {
  if (!is.numeric(T_now) || length(T_now) != 1L || !is.finite(T_now) || T_now <= 0) {
    stop("T_now must be a finite positive number.", call. = FALSE)
  }
  old_info <- .structural_partition_assignment_info(weights, T_now, control)
  logits <- (local_blocks +
               matrix(old_info$e_log_omega, nrow(local_blocks), ncol(local_blocks),
                      byrow = TRUE)) / T_now
  lse <- matrixStats::rowLogSumExps(logits)
  weights_new <- exp(logits - lse)
  info_new <- .structural_partition_assignment_info(weights_new, T_now, control)
  list(weights = weights_new, assignment_info = info_new)
}

.structural_partition_objective <- function(gamma, position_pi, weights,
                                            local_blocks, assignment_info) {
  cell_terms <- sum(vapply(seq_along(gamma), function(m) {
    .cavi_cell_terms_from_state(gamma[[m]], position_pi[[m]])
  }, numeric(1)))
  list(
    objective = as.numeric(
      cell_terms + sum(weights * local_blocks) + assignment_info$objective
    ),
    cell_terms = as.numeric(cell_terms),
    local_blocks = local_blocks,
    assignment_info = assignment_info
  )
}

.structural_partition_sweep <- function(state, T_now = 1) {
  X <- state$data
  M <- state$M
  d <- state$d
  labels <- state$ordering_labels
  q_u <- vector("list", M)
  for (m in seq_len(M)) {
    q_u[[m]] <- .structural_partition_q_u(
      X = X,
      R = state$gamma[[m]],
      sigma2 = state$sigma2,
      measurement_sd = state$measurement_sd,
      lambda_vec = state$lambda_mat[, m],
      Q_K = state$Q_K,
      ordering_label = labels[m]
    )
  }

  gamma <- vector("list", M)
  for (m in seq_len(M)) {
    gamma[[m]] <- .structural_partition_update_gamma(
      X = X,
      q_u = q_u[[m]],
      pi_vec = state$position_pi[[m]],
      sigma2 = state$sigma2,
      measurement_sd = state$measurement_sd,
      weights = state$pi_weights[, m]
    )
  }
  names(gamma) <- labels

  position_pi <- lapply(seq_len(M), function(m) {
    if (identical(state$control$position_prior, "fixed")) {
      state$position_pi[[m]]
    } else {
      value <- pmax(colMeans(gamma[[m]]), .Machine$double.eps)
      value / sum(value)
    }
  })
  names(position_pi) <- labels

  sigma_update <- NULL
  sigma2 <- state$sigma2
  if (is.null(state$measurement_sd)) {
    sigma_update <- .structural_partition_update_sigma2(
      X = X,
      gamma = gamma,
      q_u = q_u,
      weights = state$pi_weights,
      sigma_min = state$control$sigma_min,
      sigma_max = state$control$sigma_max
    )
    sigma2 <- sigma_update$sigma2
  }

  lambda_mat <- state$lambda_mat
  if (!isTRUE(state$control$fix_lambda)) {
    rank_q <- .rw_precision_metadata(state$Q_K, rw_q = state$rw_q)$rank
    for (m in seq_len(M)) {
      lambda_mat[, m] <- .cavi_update_lambda_weighted(
        q_u = q_u[[m]],
        r_rank = rank_q,
        lambda_sd_prior_rate = state$control$lambda_sd_prior_rate,
        lambda_min = state$control$lambda_min,
        lambda_max = state$control$lambda_max,
        feature_weights = rep(1, d),
        feature_active = rep(TRUE, d)
      )
    }
  }

  local_blocks <- .structural_partition_local_blocks(
    X = X,
    gamma = gamma,
    q_u = q_u,
    sigma2 = sigma2,
    measurement_sd = state$measurement_sd,
    lambda_mat = lambda_mat,
    Q_K = state$Q_K,
    rw_q = state$rw_q,
    lambda_sd_prior_rate = state$control$lambda_sd_prior_rate
  )
  weight_update <- .structural_partition_update_weights(
    local_blocks = local_blocks,
    weights = state$pi_weights,
    T_now = T_now,
    control = state$control
  )
  objective <- .structural_partition_objective(
    gamma = gamma,
    position_pi = position_pi,
    weights = weight_update$weights,
    local_blocks = local_blocks,
    assignment_info = weight_update$assignment_info
  )

  state$gamma <- gamma
  state$position_pi <- position_pi
  state$q_u <- q_u
  state$sigma2 <- sigma2
  state$lambda_mat <- lambda_mat
  state$pi_weights <- weight_update$weights
  colnames(state$pi_weights) <- labels
  state$local_blocks <- local_blocks
  colnames(state$local_blocks) <- labels
  state$assignment_info <- weight_update$assignment_info
  state$objective <- objective$objective
  state$objective_terms <- objective
  state$expected_residuals <- sigma_update$residuals %||% NULL
  state
}

.structural_partition_initial_state <- function(
    X,
    S,
    M,
    fits_init,
    init_methods,
    pca_components,
    partition_init,
    similarity_metric,
    smooth_fit_lambda_mode,
    smooth_fit_lambda_value,
    cluster_linkage,
    similarity_min_feature_sd,
    K,
    discretization,
    rw_q,
    ridge,
    lambda_init,
    fix_lambda,
    lambda_sd_prior_rate,
    lambda_min,
    lambda_max,
    sigma2_init,
    sigma_min,
    sigma_max,
    position_prior,
    position_prior_init,
    partition_prior,
    partition_prior_init,
    assignment_prior,
    ordering_alpha,
    verbose) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  labels <- .cavi_partition_order_labels(M)

  init_result <- NULL
  if (is.null(fits_init)) {
    init_result <- init_m_trajectories_cavi(
      X = X,
      S = S,
      M = M,
      methods = init_methods,
      pca_components = pca_components,
      K = K,
      rw_q = rw_q,
      ridge = ridge,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      smooth_fit_lambda_mode = smooth_fit_lambda_mode,
      smooth_fit_lambda_value = smooth_fit_lambda_value,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      discretization = discretization,
      partition_init = partition_init,
      similarity_metric = similarity_metric,
      cluster_linkage = cluster_linkage,
      similarity_min_feature_sd = similarity_min_feature_sd,
      num_iter = 2L,
      verbose = verbose
    )
    fits <- init_result$fits
  } else {
    fits <- fits_init
    if (!is.list(fits) || length(fits) != M ||
        any(!vapply(fits, inherits, logical(1), "cavi"))) {
      stop("fits_init must be a list of M cavi objects.", call. = FALSE)
    }
  }
  K_values <- vapply(fits, function(fit) ncol(fit$gamma), integer(1))
  if (length(unique(K_values)) != 1L) {
    stop("All initialization fits must use a common K.", call. = FALSE)
  }
  K_use <- K_values[1]
  if (!all(vapply(fits, function(fit) {
    isTRUE(all.equal(as.matrix(fit$data), X, check.attributes = FALSE))
  }, logical(1)))) {
    stop("All initialization fits must contain the supplied X.", call. = FALSE)
  }

  stored_S <- fits[[1]]$measurement_sd %||% NULL
  S_use <- if (is.null(S)) stored_S else S
  if (!is.null(stored_S) && !.cavi_same_measurement_sd(stored_S, S_use)) {
    stop("Supplied S does not match measurement_sd in fits_init.", call. = FALSE)
  }
  noise_info <- .cavi_resolve_measurement_sd(
    S = S_use,
    X = X,
    caller = ".structural_partition_initial_state()"
  )

  position_ctl <- .cavi_validate_position_prior(
    position_prior = position_prior,
    position_prior_init = position_prior_init,
    K = K_use,
    caller = "soft_partition_cavi()"
  )
  assignment_ctl <- .validate_partition_assignment_controls(
    partition_prior = partition_prior,
    partition_prior_init = partition_prior_init,
    assignment_prior = assignment_prior,
    ordering_alpha = ordering_alpha,
    M = M,
    partition_prior_missing = FALSE,
    caller = "soft_partition_cavi()"
  )

  gamma <- stats::setNames(lapply(fits, function(fit) {
    .structural_partition_stabilize_gamma(fit$gamma)
  }), labels)
  position_pi <- stats::setNames(lapply(gamma, function(R) {
    .cavi_initialize_position_prior(
      R = R,
      K = K_use,
      position_prior = position_ctl$position_prior,
      position_prior_init = position_ctl$position_prior_init
    )
  }), labels)
  sigma2 <- .structural_partition_sigma2(
    X = X,
    sigma2_init = sigma2_init,
    measurement_sd = noise_info$measurement_sd,
    sigma_min = sigma_min,
    sigma_max = sigma_max
  )
  lambda_mat <- .structural_partition_lambda_matrix(
    lambda_init = lambda_init,
    fits = fits,
    d = d,
    M = M,
    lambda_min = lambda_min,
    lambda_max = lambda_max
  )
  colnames(lambda_mat) <- labels
  Q_K <- make_random_walk_precision(
    K = K_use,
    d = 1,
    q = rw_q,
    lambda = 1,
    ridge = ridge
  )
  q_u <- stats::setNames(lapply(seq_len(M), function(m) {
    .structural_partition_q_u(
      X = X,
      R = gamma[[m]],
      sigma2 = sigma2,
      measurement_sd = noise_info$measurement_sd,
      lambda_vec = lambda_mat[, m],
      Q_K = Q_K,
      ordering_label = labels[m]
    )
  }), labels)
  weights <- matrix(1 / M, nrow = d, ncol = M, dimnames = list(colnames(X), labels))
  control <- list(
    variational_family = "structured",
    position_prior = position_ctl$position_prior,
    partition_prior = assignment_ctl$partition_prior,
    partition_prior_init = assignment_ctl$partition_prior_init,
    assignment_prior = assignment_ctl$legacy_assignment_prior,
    ordering_alpha = assignment_ctl$ordering_alpha,
    fix_lambda = isTRUE(fix_lambda),
    lambda_sd_prior_rate = lambda_sd_prior_rate,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    sigma_min = sigma_min,
    sigma_max = sigma_max,
    ridge = ridge,
    rw_q = rw_q,
    discretization = discretization,
    partition_init = partition_init,
    similarity_metric = similarity_metric,
    smooth_fit_lambda_mode = smooth_fit_lambda_mode,
    smooth_fit_lambda_value = smooth_fit_lambda_value,
    cluster_linkage = cluster_linkage,
    similarity_min_feature_sd = similarity_min_feature_sd,
    noise_model = noise_info$noise_model
  )
  local_blocks <- .structural_partition_local_blocks(
    X = X,
    gamma = gamma,
    q_u = q_u,
    sigma2 = sigma2,
    measurement_sd = noise_info$measurement_sd,
    lambda_mat = lambda_mat,
    Q_K = Q_K,
    rw_q = rw_q,
    lambda_sd_prior_rate = lambda_sd_prior_rate
  )
  colnames(local_blocks) <- labels
  assignment_info <- .structural_partition_assignment_info(weights, 1, control)
  objective <- .structural_partition_objective(
    gamma = gamma,
    position_pi = position_pi,
    weights = weights,
    local_blocks = local_blocks,
    assignment_info = assignment_info
  )

  list(
    data = X,
    measurement_sd = noise_info$measurement_sd,
    M = M,
    K = K_use,
    n = n,
    d = d,
    ordering_labels = labels,
    Q_K = Q_K,
    rw_q = rw_q,
    gamma = gamma,
    position_pi = position_pi,
    q_u = q_u,
    sigma2 = sigma2,
    lambda_mat = lambda_mat,
    pi_weights = weights,
    local_blocks = local_blocks,
    assignment_info = assignment_info,
    objective = objective$objective,
    objective_terms = objective,
    control = control,
    init_info = init_result$init_info %||% NULL,
    ordering_similarity = init_result$ordering_similarity %||% NULL,
    similarity_init = init_result$similarity_init %||% NULL
  )
}

.structural_partition_compatibility_fits <- function(state, objective_history) {
  labels <- state$ordering_labels
  fits <- lapply(seq_len(state$M), function(m) {
    q_u <- state$q_u[[m]]
    out <- structure(
      list(
        params = list(
          pi = state$position_pi[[m]],
          mu = lapply(seq_len(state$K), function(k) q_u$m_mat[, k]),
          sigma2 = state$sigma2
        ),
        gamma = state$gamma[[m]],
        measurement_sd = state$measurement_sd,
        posterior = list(
          mean = q_u$m_mat,
          cov = q_u$S_list,
          var = q_u$sdiag_mat,
          diag = q_u$sdiag_mat
        ),
        lambda_vec = state$lambda_mat[, m],
        Q_K = state$Q_K,
        rw_q = state$rw_q,
        elbo_trace = numeric(0),
        loglik_trace = numeric(0),
        ml_trace = numeric(0),
        lambda_trace = list(state$lambda_mat[, m]),
        sigma2_trace = list(state$sigma2),
        pi_trace = list(state$position_pi[[m]]),
        iter = max(0L, length(objective_history) - 1L),
        converged = FALSE,
        control = utils::modifyList(
          state$control,
          list(
            modelName = "homoskedastic",
            compatibility_view = TRUE,
            ordering = labels[m]
          )
        ),
        data = state$data
      ),
      class = "cavi"
    )
    out$priors <- .mpcurve_priors_from_cavi_fit(out)
    out
  })
  stats::setNames(fits, labels)
}

.structural_partition_finalize <- function(state, histories, converged) {
  labels <- state$ordering_labels
  means <- stats::setNames(lapply(state$q_u, `[[`, "m_mat"), labels)
  covs <- stats::setNames(lapply(state$q_u, `[[`, "S_list"), labels)
  vars <- stats::setNames(lapply(state$q_u, `[[`, "sdiag_mat"), labels)
  fits <- .structural_partition_compatibility_fits(state, histories$objective)
  prior_meta <- .rw_precision_metadata(state$Q_K, rw_q = state$rw_q)
  structure(
    list(
      data = state$data,
      params = list(
        pi = state$position_pi,
        mu = means,
        sigma2 = state$sigma2
      ),
      measurement_sd = state$measurement_sd,
      posterior = list(mean = means, cov = covs, var = vars, diag = vars),
      conditional_posterior = list(mean = means, cov = covs, var = vars, diag = vars),
      gamma = state$gamma,
      position_pi = state$position_pi,
      lambda_mat = state$lambda_mat,
      Q_K = state$Q_K,
      rw_q = state$rw_q,
      pi_weights = state$pi_weights,
      assign = labels[max.col(state$pi_weights, ties.method = "first")],
      local_blocks = state$local_blocks,
      assignment_posterior = .cavi_partition_assignment_state(state$assignment_info),
      objective_terms = state$objective_terms,
      objective_history = histories$objective,
      temperature_history = histories$temperature,
      n_anneal = 1L + state$control$n_outer * state$control$inner_iter,
      weight_history = histories$weights,
      score_history = histories$local_blocks,
      sigma2_trace = histories$sigma2,
      lambda_trace = histories$lambda,
      iter = max(0L, length(histories$objective) - 1L),
      M = state$M,
      K = state$K,
      n = state$n,
      d = state$d,
      ordering_labels = labels,
      variational_family = "structured",
      fits = fits,
      init_info = state$init_info,
      ordering_similarity = state$ordering_similarity,
      similarity_init = state$similarity_init,
      converged = isTRUE(converged),
      convergence_info = list(
        tol_outer = state$control$tol_outer,
        iterations = max(0L, length(histories$objective) - 1L),
        converged = isTRUE(converged)
      ),
      control = utils::modifyList(
        state$control,
        list(
          variational_family = "structured",
          prior_proper = prior_meta$proper,
          prior_rank = prior_meta$rank,
          prior_logdet = prior_meta$logdet
        )
      )
    ),
    class = "soft_partition_cavi"
  )
}

.structural_partition_append_history <- function(histories, state, T_now) {
  histories$objective <- c(histories$objective, state$objective)
  histories$temperature <- c(histories$temperature, T_now)
  histories$weights[length(histories$weights) + 1L] <- list(state$pi_weights)
  histories$local_blocks[length(histories$local_blocks) + 1L] <- list(state$local_blocks)
  histories$sigma2[length(histories$sigma2) + 1L] <- list(state$sigma2)
  histories$lambda[length(histories$lambda) + 1L] <- list(state$lambda_mat)
  histories
}

.structural_partition_cavi <- function(
    X,
    S = NULL,
    M = 2L,
    fits_init = NULL,
    init_methods = NULL,
    pca_components = NULL,
    partition_init = c("similarity", "ordering_methods"),
    similarity_metric = c("spearman", "pearson", "smooth_fit"),
    smooth_fit_lambda_mode = c("optimize", "fixed"),
    smooth_fit_lambda_value = 1,
    cluster_linkage = "single",
    similarity_min_feature_sd = 1e-8,
    K = NULL,
    discretization = c("quantile", "equal", "kmeans"),
    T_start = 5,
    T_end = 1,
    n_outer = 25L,
    inner_iter = 1L,
    max_converge_iter = 100L,
    tol_outer = 1e-5,
    rw_q = 2L,
    ridge = 0,
    lambda_init = 1,
    fix_lambda = FALSE,
    lambda_sd_prior_rate = NULL,
    lambda_min = 1e-10,
    lambda_max = 1e10,
    sigma2_init = NULL,
    sigma_min = 1e-10,
    sigma_max = 1e10,
    position_prior = c("adaptive", "fixed"),
    position_prior_init = NULL,
    partition_prior = c("adaptive", "fixed"),
    partition_prior_init = NULL,
    assignment_prior = NULL,
    ordering_alpha = NULL,
    hard_assign_final = FALSE,
    freeze_unused_ordering = NULL,
    freeze_unused_ordering_threshold = NULL,
    freeze_feature = NULL,
    freeze_feature_weight_threshold = NULL,
    drop_unused_ordering = NULL,
    verbose = TRUE) {
  .structural_partition_warn_legacy_controls(
    freeze_unused_ordering = freeze_unused_ordering,
    freeze_unused_ordering_threshold = freeze_unused_ordering_threshold,
    freeze_feature = freeze_feature,
    freeze_feature_weight_threshold = freeze_feature_weight_threshold,
    drop_unused_ordering = drop_unused_ordering,
    caller = "soft_partition_cavi()"
  )
  X <- as.matrix(X)
  M <- as.integer(M)
  n_outer <- as.integer(n_outer)
  inner_iter <- as.integer(inner_iter)
  max_converge_iter <- as.integer(max_converge_iter)
  if (M < 2L) stop("M must be >= 2 for partition models.", call. = FALSE)
  if (n_outer < 1L || inner_iter < 1L || max_converge_iter < 0L) {
    stop("n_outer and inner_iter must be positive; max_converge_iter must be nonnegative.",
         call. = FALSE)
  }
  if (!is.finite(T_start) || !is.finite(T_end) || T_start <= 0 || T_end <= 0) {
    stop("T_start and T_end must be finite positive values.", call. = FALSE)
  }
  if (!is.finite(tol_outer) || tol_outer < 0) {
    stop("tol_outer must be a finite nonnegative value.", call. = FALSE)
  }
  partition_init <- match.arg(partition_init)
  similarity_metric <- match.arg(similarity_metric)
  smooth_fit_lambda_mode <- match.arg(smooth_fit_lambda_mode)
  discretization <- match.arg(discretization)
  position_prior <- match.arg(position_prior)
  partition_prior <- match.arg(partition_prior)
  lambda_sd_prior_rate <- .normalize_lambda_sd_prior_rate(lambda_sd_prior_rate)
  cluster_linkage <- .cavi_validate_cluster_linkage(cluster_linkage)

  state <- .structural_partition_initial_state(
    X = X,
    S = S,
    M = M,
    fits_init = fits_init,
    init_methods = init_methods,
    pca_components = pca_components,
    partition_init = partition_init,
    similarity_metric = similarity_metric,
    smooth_fit_lambda_mode = smooth_fit_lambda_mode,
    smooth_fit_lambda_value = smooth_fit_lambda_value,
    cluster_linkage = cluster_linkage,
    similarity_min_feature_sd = similarity_min_feature_sd,
    K = K,
    discretization = discretization,
    rw_q = rw_q,
    ridge = ridge,
    lambda_init = lambda_init,
    fix_lambda = fix_lambda,
    lambda_sd_prior_rate = lambda_sd_prior_rate,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    sigma2_init = sigma2_init,
    sigma_min = sigma_min,
    sigma_max = sigma_max,
    position_prior = position_prior,
    position_prior_init = position_prior_init,
    partition_prior = partition_prior,
    partition_prior_init = partition_prior_init,
    assignment_prior = assignment_prior,
    ordering_alpha = ordering_alpha,
    verbose = verbose
  )
  state$control$tol_outer <- tol_outer
  state$control$T_start <- T_start
  state$control$T_end <- T_end
  state$control$n_outer <- n_outer
  state$control$inner_iter <- inner_iter
  state$control$max_converge_iter <- max_converge_iter
  state$assignment_info <- .structural_partition_assignment_info(
    state$pi_weights,
    T_start,
    state$control
  )
  state$objective_terms <- .structural_partition_objective(
    gamma = state$gamma,
    position_pi = state$position_pi,
    weights = state$pi_weights,
    local_blocks = state$local_blocks,
    assignment_info = state$assignment_info
  )
  state$objective <- state$objective_terms$objective
  histories <- list(
    objective = state$objective,
    temperature = T_start,
    weights = list(state$pi_weights),
    local_blocks = list(state$local_blocks),
    sigma2 = list(state$sigma2),
    lambda = list(state$lambda_mat)
  )

  schedule <- if (n_outer == 1L) {
    T_end
  } else {
    exp(seq(log(T_start), log(T_end), length.out = n_outer))
  }
  for (T_now in schedule) {
    for (inner in seq_len(inner_iter)) {
      state <- .structural_partition_sweep(state, T_now = T_now)
      histories <- .structural_partition_append_history(histories, state, T_now)
    }
  }

  if (!isTRUE(all.equal(tail(histories$temperature, 1L), 1))) {
    state$assignment_info <- .structural_partition_assignment_info(
      state$pi_weights,
      1,
      state$control
    )
    state$objective_terms <- .structural_partition_objective(
      gamma = state$gamma,
      position_pi = state$position_pi,
      weights = state$pi_weights,
      local_blocks = state$local_blocks,
      assignment_info = state$assignment_info
    )
    state$objective <- state$objective_terms$objective
    histories <- .structural_partition_append_history(histories, state, 1)
  }

  converged <- FALSE
  previous <- tail(histories$objective, 1L)
  if (max_converge_iter > 0L) {
    for (iter in seq_len(max_converge_iter)) {
      state <- .structural_partition_sweep(state, T_now = 1)
      histories <- .structural_partition_append_history(histories, state, 1)
      current <- state$objective
      delta <- current - previous
      rel_delta <- delta / (abs(previous) + 1e-12)
      if (delta < -1e-7 * (abs(previous) + 1)) {
        warning(
          sprintf("Structural partition ELBO decreased by %.3e at T = 1 iteration %d.",
                  delta, iter),
          call. = FALSE
        )
      }
      if (delta >= -1e-8 * (abs(previous) + 1) && abs(rel_delta) < tol_outer) {
        converged <- TRUE
        break
      }
      previous <- current
    }
  }

  if (isTRUE(hard_assign_final)) {
    hard <- matrix(0, nrow = state$d, ncol = state$M)
    hard[cbind(seq_len(state$d), max.col(state$pi_weights, ties.method = "first"))] <- 1
    colnames(hard) <- state$ordering_labels
    state$pi_weights <- hard
    state$assignment_info <- .structural_partition_assignment_info(hard, 1, state$control)
    state$objective_terms <- .structural_partition_objective(
      gamma = state$gamma,
      position_pi = state$position_pi,
      weights = hard,
      local_blocks = state$local_blocks,
      assignment_info = state$assignment_info
    )
    state$objective <- state$objective_terms$objective
    histories <- .structural_partition_append_history(histories, state, 1)
  }

  .structural_partition_finalize(state, histories, converged)
}

.structural_partition_state_from_fit <- function(fit) {
  if (!inherits(fit, "soft_partition_cavi") ||
      !identical(fit$variational_family %||% fit$control$variational_family,
                 "structured")) {
    stop("fit is not a structural partition-CAVI object.", call. = FALSE)
  }
  labels <- fit$ordering_labels %||% colnames(fit$pi_weights) %||%
    .cavi_partition_order_labels(fit$M)
  q_u <- stats::setNames(lapply(seq_len(fit$M), function(m) {
    list(
      m_mat = fit$conditional_posterior$mean[[m]],
      S_list = fit$conditional_posterior$cov[[m]],
      sdiag_mat = fit$conditional_posterior$var[[m]],
      logdetS = vapply(
        fit$conditional_posterior$cov[[m]],
        function(S) as.numeric(determinant(S, logarithm = TRUE)$modulus),
        numeric(1)
      ),
      eq_quad = vapply(seq_len(fit$d), function(j) {
        mu <- fit$conditional_posterior$mean[[m]][j, ]
        V <- fit$conditional_posterior$cov[[m]][[j]]
        as.numeric(crossprod(mu, fit$Q_K %*% mu) + sum(fit$Q_K * V))
      }, numeric(1))
    )
  }), labels)
  list(
    data = fit$data,
    measurement_sd = fit$measurement_sd,
    M = fit$M,
    K = fit$K,
    n = fit$n,
    d = fit$d,
    ordering_labels = labels,
    Q_K = fit$Q_K,
    rw_q = fit$rw_q,
    gamma = fit$gamma,
    position_pi = fit$position_pi,
    q_u = q_u,
    sigma2 = fit$params$sigma2,
    lambda_mat = fit$lambda_mat,
    pi_weights = fit$pi_weights,
    local_blocks = fit$local_blocks,
    assignment_info = .structural_partition_assignment_info(fit$pi_weights, 1, fit$control),
    objective = tail(fit$objective_history, 1L),
    objective_terms = fit$objective_terms,
    control = fit$control,
    init_info = fit$init_info,
    ordering_similarity = fit$ordering_similarity,
    similarity_init = fit$similarity_init
  )
}

.continue_structural_partition_cavi <- function(
    fit,
    iter = 1L,
    lambda = NULL,
    S = NULL,
    tol_outer = NULL,
    lambda_sd_prior_rate = NULL,
    lambda_min = NULL,
    lambda_max = NULL,
    sigma_min = NULL,
    sigma_max = NULL,
    verbose = FALSE) {
  iter <- as.integer(iter)
  if (iter < 1L) stop("iter must be at least 1.", call. = FALSE)
  state <- .structural_partition_state_from_fit(fit)
  if (!is.null(S) && !.cavi_same_measurement_sd(S, state$measurement_sd)) {
    stop("Supplied S does not match measurement_sd stored on the fit.", call. = FALSE)
  }
  if (!is.null(lambda)) {
    state$lambda_mat <- .structural_partition_lambda_matrix(
      lambda_init = lambda,
      fits = fit$fits,
      d = state$d,
      M = state$M,
      lambda_min = lambda_min %||% state$control$lambda_min,
      lambda_max = lambda_max %||% state$control$lambda_max
    )
  }
  state$control$lambda_sd_prior_rate <- if (is.null(lambda_sd_prior_rate)) {
    state$control$lambda_sd_prior_rate
  } else {
    .normalize_lambda_sd_prior_rate(lambda_sd_prior_rate)
  }
  state$control$lambda_min <- lambda_min %||% state$control$lambda_min
  state$control$lambda_max <- lambda_max %||% state$control$lambda_max
  state$control$sigma_min <- sigma_min %||% state$control$sigma_min
  state$control$sigma_max <- sigma_max %||% state$control$sigma_max
  state$control$tol_outer <- tol_outer %||% state$control$tol_outer

  histories <- list(
    objective = fit$objective_history,
    temperature = fit$temperature_history,
    weights = fit$weight_history,
    local_blocks = fit$score_history,
    sigma2 = fit$sigma2_trace,
    lambda = fit$lambda_trace
  )
  converged <- FALSE
  previous <- tail(histories$objective, 1L)
  for (i in seq_len(iter)) {
    state <- .structural_partition_sweep(state, T_now = 1)
    histories <- .structural_partition_append_history(histories, state, 1)
    delta <- state$objective - previous
    rel_delta <- delta / (abs(previous) + 1e-12)
    if (delta < -1e-7 * (abs(previous) + 1)) {
      warning(
        sprintf("Structural partition ELBO decreased by %.3e during continuation.", delta),
        call. = FALSE
      )
    }
    if (delta >= -1e-8 * (abs(previous) + 1) &&
        abs(rel_delta) < state$control$tol_outer) {
      converged <- TRUE
      break
    }
    previous <- state$objective
    if (isTRUE(verbose)) {
      message(sprintf("[structural partition] continuation %d objective %.6f", i,
                      state$objective))
    }
  }
  .structural_partition_finalize(state, histories, converged)
}

# Internal fixed-M partition entry point used by fit_mpcurve() and the
# two-ordering compatibility wrapper. The legacy augmented implementation has
# a distinct archival name in 10_partition_cavi.R.
soft_partition_cavi <- .structural_partition_cavi
