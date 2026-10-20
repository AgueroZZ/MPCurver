`%||%` <- function(a, b) if (!is.null(a)) a else b


# ============================================================
# Internal helpers for multi-ordering CAVI: feature similarities and
# similarity-based initialization, per-ordering starting fits,
# feature-assignment controls, and objective terms used by structural
# inference.
# ============================================================

.cavi_partition_weight_floor <- function() sqrt(.Machine$double.eps)

.cavi_resolve_K <- function(X, K, rw_q = 2L) {
  rw_q <- .mpcurve_integer_setting(rw_q, "rw_order")
  if (is.null(K)) {
    if (nrow(X) <= rw_q) {
      stop("Automatic num_bins requires at least rw_order + 1 samples.", call. = FALSE)
    }
    return(max(rw_q + 1L, min(50L, as.integer(floor(nrow(X) / 5)))))
  }
  K <- .mpcurve_integer_setting(K, "num_bins", 2L)
  if (K <= rw_q) stop("num_bins must be at least rw_order + 1.", call. = FALSE)
  K
}

.cavi_resolve_feature_weights <- function(feature_weights, d,
                                          floor = .cavi_partition_weight_floor()) {
  if (is.null(feature_weights)) return(rep(1, d))
  w <- pmax(as.numeric(feature_weights), floor)
  if (length(w) != d) stop("feature_weights must have length d.")
  w
}

.cavi_subset_measurement_sd <- function(S, cols) {
  if (is.null(S)) return(NULL)
  if (is.matrix(S)) {
    return(as.matrix(S)[, cols, drop = FALSE])
  }
  as.numeric(S)[cols]
}

.cavi_partition_order_labels <- function(M) {
  if (M <= 26L) LETTERS[seq_len(M)] else paste0("ord", seq_len(M))
}

.cavi_default_subset_pca_components <- function(methods) {
  methods <- as.character(methods)
  out <- rep(NA_integer_, length(methods))
  out[methods == "PCA"] <- 1L
  out
}

.cavi_validate_cluster_linkage <- function(cluster_linkage = "single") {
  valid <- c("ward.D", "ward.D2", "single", "complete",
             "average", "mcquitty", "median", "centroid")
  cluster_linkage <- as.character(cluster_linkage)[1]
  if (is.na(cluster_linkage) || !nzchar(cluster_linkage) || !cluster_linkage %in% valid) {
    stop(
      "cluster_linkage must be one of: ",
      paste(valid, collapse = ", "),
      "."
    )
  }
  cluster_linkage
}

.cavi_constant_null_score <- function(x,
                                      sigma_min = 1e-10,
                                      sigma_max = 1e10) {
  x <- as.numeric(x)
  sigma2 <- mean((x - mean(x))^2)
  sigma2 <- pmin(pmax(sigma2, sigma_min), sigma_max)
  -0.5 * length(x) * (log(2 * base::pi * sigma2) + 1)
}

.cavi_constant_null_sigma2 <- function(x,
                                       sigma_min = 1e-10,
                                       sigma_max = 1e10) {
  x <- as.numeric(x)
  sigma2 <- mean((x - mean(x))^2)
  pmin(pmax(sigma2, sigma_min), sigma_max)
}

.cavi_constant_null_score_known <- function(x, measurement_sd) {
  x <- as.numeric(x)
  sd_vec <- pmax(as.numeric(measurement_sd), sqrt(.Machine$double.eps))
  if (length(sd_vec) == 1L) sd_vec <- rep(sd_vec, length(x))
  if (length(sd_vec) != length(x)) {
    stop("measurement_sd must have length equal to length(x).")
  }
  var_vec <- sd_vec^2
  prec_vec <- 1 / var_vec
  mu_hat <- sum(prec_vec * x) / sum(prec_vec)
  as.numeric(
    -0.5 * (
      sum(log(2 * base::pi * var_vec)) +
        sum((x - mu_hat)^2 * prec_vec)
    )
  )
}

.cavi_constant_null_sigma2_known <- function(measurement_sd) {
  mean(pmax(as.numeric(measurement_sd), sqrt(.Machine$double.eps))^2)
}

.compute_same_ordering_similarity_cor <- function(X,
                                                  metric = c("spearman", "pearson"),
                                                  use = "pairwise.complete.obs",
                                                  abs_value = TRUE,
                                                  min_feature_sd = 1e-8) {
  X <- as.matrix(X)
  metric <- match.arg(metric)
  min_feature_sd <- as.numeric(min_feature_sd)[1]
  if (!is.finite(min_feature_sd) || min_feature_sd < 0) {
    stop("min_feature_sd must be a single finite nonnegative number.")
  }

  d <- ncol(X)
  feature_names <- colnames(X) %||% paste0("V", seq_len(d))
  feature_sd <- apply(X, 2, stats::sd, na.rm = TRUE)
  low_variance <- !is.finite(feature_sd) | feature_sd < min_feature_sd

  if (d == 1L) {
    S <- matrix(1, nrow = 1L, ncol = 1L)
  } else {
    S <- suppressWarnings(stats::cor(X, method = metric, use = use))
    if (!is.matrix(S)) {
      S <- matrix(S, nrow = d, ncol = d)
    }
    S[!is.finite(S)] <- 0
    if (isTRUE(abs_value)) {
      S <- abs(S)
    }
    S <- pmin(pmax(S, 0), 1)
    S <- 0.5 * (S + t(S))
  }

  dimnames(S) <- list(feature_names, feature_names)
  if (any(low_variance)) {
    S[low_variance, ] <- 0
    S[, low_variance] <- 0
  }
  diag(S) <- 1

  list(
    S = S,
    distance = 1 - S,
    metric = metric,
    feature_info = data.frame(
      feature = feature_names,
      sd = as.numeric(feature_sd),
      low_variance = low_variance,
      stringsAsFactors = FALSE
    )
  )
}

.cavi_equal_width_bin_index <- function(x, K) {
  x <- as.numeric(x)
  K <- as.integer(K)[1]
  xr <- range(x, na.rm = TRUE)
  if (!is.finite(xr[1]) || !is.finite(xr[2])) {
    stop("x must contain only finite values.")
  }
  if (K < 2L) stop("K must be >= 2.")
  if (xr[1] == xr[2]) {
    return(list(
      bin_index = rep(as.integer(ceiling(K / 2)), length(x)),
      breaks = c(xr[1], rep(xr[2], K))
    ))
  }
  breaks <- seq(xr[1], xr[2], length.out = K + 1L)
  idx <- findInterval(x, vec = breaks, rightmost.closed = TRUE, all.inside = TRUE)
  idx <- pmin(pmax(as.integer(idx), 1L), K)
  list(bin_index = idx, breaks = breaks)
}

.cavi_binned_directional_stats <- function(x_cov, y, K) {
  bin_info <- .cavi_equal_width_bin_index(x_cov, K = K)
  idx <- bin_info$bin_index
  counts <- tabulate(idx, nbins = K)
  y_sum <- numeric(K)
  grouped <- rowsum(matrix(as.numeric(y), ncol = 1L), group = idx, reorder = FALSE)
  y_sum[as.integer(rownames(grouped))] <- grouped[, 1]
  list(
    bin_index = idx,
    breaks = bin_info$breaks,
    counts = as.numeric(counts),
    y_sum = as.numeric(y_sum),
    y_sq = sum(as.numeric(y)^2),
    n = length(y)
  )
}

.cavi_binned_directional_stats_known <- function(x_cov, y, measurement_sd, K) {
  bin_info <- .cavi_equal_width_bin_index(x_cov, K = K)
  idx <- bin_info$bin_index
  counts <- tabulate(idx, nbins = K)
  sd_vec <- pmax(as.numeric(measurement_sd), sqrt(.Machine$double.eps))
  if (length(sd_vec) == 1L) sd_vec <- rep(sd_vec, length(y))
  if (length(sd_vec) != length(y)) {
    stop("measurement_sd must have length equal to length(y).")
  }
  var_vec <- sd_vec^2
  prec_vec <- 1 / var_vec

  prec_sum <- numeric(K)
  wy_sum <- numeric(K)
  prec_grouped <- rowsum(matrix(prec_vec, ncol = 1L), group = idx, reorder = FALSE)
  wy_grouped <- rowsum(matrix(prec_vec * as.numeric(y), ncol = 1L), group = idx, reorder = FALSE)
  prec_sum[as.integer(rownames(prec_grouped))] <- prec_grouped[, 1]
  wy_sum[as.integer(rownames(wy_grouped))] <- wy_grouped[, 1]

  list(
    bin_index = idx,
    breaks = bin_info$breaks,
    counts = as.numeric(counts),
    prec_sum = as.numeric(prec_sum),
    wy_sum = as.numeric(wy_sum),
    wy_sq = sum((as.numeric(y)^2) * prec_vec),
    log_const = sum(log(2 * base::pi * var_vec)),
    n = length(y)
  )
}

.cavi_smooth_fit_log_evidence <- function(stats,
                                          lambda,
                                          sigma2,
                                          Q_K,
                                          logdet_Q,
                                          rank_Q) {
  K <- nrow(Q_K)
  b <- stats$y_sum / sigma2
  P <- lambda * Q_K + diag(stats$counts / sigma2, K)
  chol_info <- .cavi_chol_with_jitter(P)
  L <- chol_info$chol
  alpha <- backsolve(L, forwardsolve(t(L), b))
  logdet_P <- 2 * sum(log(diag(L)))
  as.numeric(
    -0.5 * stats$n * log(2 * base::pi * sigma2) +
      0.5 * (rank_Q * log(lambda) + logdet_Q) -
      0.5 * logdet_P -
      0.5 * stats$y_sq / sigma2 +
      0.5 * sum(b * alpha)
  )
}

.cavi_smooth_fit_log_evidence_known <- function(stats,
                                                lambda,
                                                Q_K,
                                                logdet_Q,
                                                rank_Q) {
  K <- nrow(Q_K)
  b <- stats$wy_sum
  P <- lambda * Q_K + diag(stats$prec_sum, K)
  chol_info <- .cavi_chol_with_jitter(P)
  L <- chol_info$chol
  alpha <- backsolve(L, forwardsolve(t(L), b))
  logdet_P <- 2 * sum(log(diag(L)))
  as.numeric(
    -0.5 * stats$log_const +
      0.5 * (rank_Q * log(lambda) + logdet_Q) -
      0.5 * logdet_P -
      0.5 * stats$wy_sq +
      0.5 * sum(b * alpha)
  )
}

.cavi_fit_directional_smooth_evidence <- function(x_cov,
                                                  y,
                                                  measurement_sd = NULL,
                                                  K,
                                                  rw_q = 2L,
                                                  ridge = 0,
                                                  lambda_sd_prior_rate = NULL,
                                                  lambda_mode = c("optimize", "fixed"),
                                                  lambda_value = 1,
                                                  lambda_min = 1e-10,
                                                  lambda_max = 1e10,
                                                  sigma_min = 1e-10,
                                                  sigma_max = 1e10) {
  lambda_mode <- match.arg(lambda_mode)
  lambda_value <- as.numeric(lambda_value)[1]
  if (!is.finite(lambda_value) || lambda_value <= 0) {
    stop("smooth_fit_lambda_value must be a single finite positive number.")
  }
  if (!is.finite(sigma_min) || !is.finite(sigma_max) || sigma_min <= 0 || sigma_max <= 0 || sigma_min > sigma_max) {
    stop("sigma_min/sigma_max must be positive finite numbers with sigma_min <= sigma_max.")
  }
  Q_K <- make_random_walk_precision(K = K, d = 1, q = rw_q, lambda = 1, ridge = ridge)
  Q_K <- 0.5 * (as.matrix(Q_K) + t(as.matrix(Q_K)))
  prior_meta <- .rw_precision_metadata(Q_K, rw_q = rw_q)
  has_known_sd <- !is.null(measurement_sd)
  if (has_known_sd) {
    measurement_sd <- pmax(as.numeric(measurement_sd), sqrt(.Machine$double.eps))
    if (length(measurement_sd) == 1L) measurement_sd <- rep(measurement_sd, length(y))
    if (length(measurement_sd) != length(y)) {
      stop("measurement_sd must have length equal to length(y).")
    }
    null_sigma2 <- .cavi_constant_null_sigma2_known(measurement_sd)
    null_score <- .cavi_constant_null_score_known(y, measurement_sd = measurement_sd)
    stats <- .cavi_binned_directional_stats_known(
      x_cov = x_cov,
      y = y,
      measurement_sd = measurement_sd,
      K = K
    )
  } else {
    null_sigma2 <- .cavi_constant_null_sigma2(y, sigma_min = sigma_min, sigma_max = sigma_max)
    null_score <- .cavi_constant_null_score(y, sigma_min = sigma_min, sigma_max = sigma_max)
    stats <- .cavi_binned_directional_stats(x_cov = x_cov, y = y, K = K)
  }
  rate_val <- .lambda_sd_prior_rate_value(lambda_sd_prior_rate)

  score_at <- function(lambda, sigma2) {
    if (has_known_sd) {
      .cavi_smooth_fit_log_evidence_known(
        stats = stats,
        lambda = lambda,
        Q_K = Q_K,
        logdet_Q = prior_meta$logdet,
        rank_Q = prior_meta$rank
      )
    } else {
      .cavi_smooth_fit_log_evidence(
        stats = stats,
        lambda = lambda,
        sigma2 = sigma2,
        Q_K = Q_K,
        logdet_Q = prior_meta$logdet,
        rank_Q = prior_meta$rank
      )
    }
  }

  objective_at <- function(lambda, sigma2) {
    score <- score_at(lambda, sigma2)
    if (rate_val > 0) {
      score <- score + .lambda_sd_prior_terms(
        lambda_vec = lambda,
        rate = rate_val,
        include_constant = TRUE
      )
    }
    score
  }

  lambda_init <- pmin(pmax(lambda_value, lambda_min), lambda_max)
  sigma_init <- null_sigma2

  if (has_known_sd) {
    if (identical(lambda_mode, "fixed")) {
      log_evidence <- score_at(lambda = lambda_init, sigma2 = null_sigma2)
      return(list(
        log_evidence = as.numeric(log_evidence),
        null_score = as.numeric(null_score),
        delta = as.numeric(max(log_evidence - null_score, 0)),
        lambda = as.numeric(lambda_init),
        sigma2 = as.numeric(null_sigma2),
        null_sigma2 = as.numeric(null_sigma2),
        success = TRUE,
        breaks = stats$breaks,
        bin_counts = stats$counts
      ))
    }

    opt_lambda <- tryCatch(
      stats::optimize(
        f = function(log_lambda) {
          -objective_at(lambda = exp(log_lambda), sigma2 = null_sigma2)
        },
        interval = c(log(lambda_min), log(lambda_max))
      ),
      error = function(e) NULL
    )
    if (is.null(opt_lambda)) {
      return(list(
        log_evidence = as.numeric(null_score),
        null_score = as.numeric(null_score),
        delta = 0,
        lambda = as.numeric(lambda_init),
        sigma2 = as.numeric(null_sigma2),
        null_sigma2 = as.numeric(null_sigma2),
        success = FALSE,
        breaks = stats$breaks,
        bin_counts = stats$counts
      ))
    }

    lambda_hat <- exp(opt_lambda$minimum)
    log_evidence <- score_at(lambda = lambda_hat, sigma2 = null_sigma2)
    return(list(
      log_evidence = as.numeric(log_evidence),
      null_score = as.numeric(null_score),
      delta = as.numeric(max(log_evidence - null_score, 0)),
      lambda = as.numeric(lambda_hat),
      sigma2 = as.numeric(null_sigma2),
      null_sigma2 = as.numeric(null_sigma2),
      success = TRUE,
      breaks = stats$breaks,
      bin_counts = stats$counts
    ))
  }

  if (identical(lambda_mode, "fixed")) {
    opt_sigma <- stats::optimize(
      f = function(log_sigma2) {
        -score_at(lambda = lambda_init, sigma2 = exp(log_sigma2))
      },
      interval = c(log(sigma_min), log(sigma_max))
    )
    sigma_hat <- exp(opt_sigma$minimum)
    log_evidence <- score_at(lambda = lambda_init, sigma2 = sigma_hat)
    return(list(
      log_evidence = as.numeric(log_evidence),
      null_score = as.numeric(null_score),
      delta = as.numeric(max(log_evidence - null_score, 0)),
      lambda = as.numeric(lambda_init),
      sigma2 = as.numeric(sigma_hat),
      null_sigma2 = as.numeric(null_sigma2),
      success = TRUE,
      breaks = stats$breaks,
      bin_counts = stats$counts
    ))
  }

  opt_res <- tryCatch(
    stats::optim(
      par = c(log(lambda_init), log(sigma_init)),
      fn = function(theta) {
        lambda <- exp(theta[1])
        sigma2 <- exp(theta[2])
        -objective_at(lambda = lambda, sigma2 = sigma2)
      },
      method = "L-BFGS-B",
      lower = c(log(lambda_min), log(sigma_min)),
      upper = c(log(lambda_max), log(sigma_max))
    ),
    error = function(e) NULL
  )

  if (is.null(opt_res)) {
    return(list(
      log_evidence = as.numeric(null_score),
      null_score = as.numeric(null_score),
      delta = 0,
      lambda = as.numeric(lambda_init),
      sigma2 = as.numeric(null_sigma2),
      null_sigma2 = as.numeric(null_sigma2),
      success = FALSE,
      breaks = stats$breaks,
      bin_counts = stats$counts
    ))
  }

  lambda_hat <- exp(opt_res$par[1])
  sigma_hat <- exp(opt_res$par[2])
  log_evidence <- score_at(lambda = lambda_hat, sigma2 = sigma_hat)

  list(
    log_evidence = as.numeric(log_evidence),
    null_score = as.numeric(null_score),
    delta = as.numeric(max(log_evidence - null_score, 0)),
    lambda = as.numeric(lambda_hat),
    sigma2 = as.numeric(sigma_hat),
    null_sigma2 = as.numeric(null_sigma2),
    success = TRUE,
    breaks = stats$breaks,
    bin_counts = stats$counts
  )
}

.compute_same_ordering_similarity_smooth_fit <- function(X,
                                                         S = NULL,
                                                         K,
                                                         rw_q = 2L,
                                                         ridge = 0,
                                                         lambda_sd_prior_rate = NULL,
                                                         smooth_fit_lambda_mode = c("optimize", "fixed"),
                                                         smooth_fit_lambda_value = 1,
                                                         lambda_min = 1e-10,
                                                         lambda_max = 1e10,
                                                         min_feature_sd = 1e-8,
                                                         sigma_min = 1e-10,
                                                         sigma_max = 1e10) {
  X <- as.matrix(X)
  smooth_fit_lambda_mode <- match.arg(smooth_fit_lambda_mode)
  min_feature_sd <- as.numeric(min_feature_sd)[1]
  if (!is.finite(min_feature_sd) || min_feature_sd < 0) {
    stop("min_feature_sd must be a single finite nonnegative number.")
  }
  if (!is.finite(as.numeric(smooth_fit_lambda_value)[1]) ||
      as.numeric(smooth_fit_lambda_value)[1] <= 0) {
    stop("smooth_fit_lambda_value must be a single finite positive number.")
  }

  n <- nrow(X)
  d <- ncol(X)
  noise_info <- .cavi_resolve_measurement_sd(S = S, X = X,
                                             caller = ".compute_same_ordering_similarity_smooth_fit()")
  feature_names <- colnames(X) %||% paste0("V", seq_len(d))
  feature_sd <- apply(X, 2, stats::sd, na.rm = TRUE)
  low_variance <- !is.finite(feature_sd) | feature_sd < min_feature_sd
  observed <- !is.na(X)
  if (is.null(noise_info$measurement_sd)) {
    null_sigma2 <- vapply(seq_len(d), function(k) {
      .cavi_constant_null_sigma2(X[observed[, k], k], sigma_min = sigma_min, sigma_max = sigma_max)
    }, numeric(1))
    null_scores <- vapply(seq_len(d), function(k) {
      .cavi_constant_null_score(X[observed[, k], k], sigma_min = sigma_min, sigma_max = sigma_max)
    }, numeric(1))
  } else {
    null_sigma2 <- vapply(seq_len(d), function(k) {
      .cavi_constant_null_sigma2_known(noise_info$sd_mat[observed[, k], k])
    }, numeric(1))
    null_scores <- vapply(seq_len(d), function(k) {
      .cavi_constant_null_score_known(X[observed[, k], k],
                                      measurement_sd = noise_info$sd_mat[observed[, k], k])
    }, numeric(1))
  }
  lambda_default <- pmin(pmax(as.numeric(smooth_fit_lambda_value)[1], lambda_min), lambda_max)

  if (d == 1L) {
    S <- matrix(1, nrow = 1L, ncol = 1L, dimnames = list(feature_names, feature_names))
    distance <- matrix(0, nrow = 1L, ncol = 1L, dimnames = list(feature_names, feature_names))
    directional_score <- matrix(null_scores, nrow = 1L, ncol = 1L, dimnames = list(feature_names, feature_names))
    directional_delta <- matrix(0, nrow = 1L, ncol = 1L, dimnames = list(feature_names, feature_names))
    directional_success <- matrix(FALSE, nrow = 1L, ncol = 1L, dimnames = list(feature_names, feature_names))
    directional_lambda <- matrix(lambda_default, nrow = 1L, ncol = 1L, dimnames = list(feature_names, feature_names))
    directional_sigma2 <- matrix(null_sigma2, nrow = 1L, ncol = 1L, dimnames = list(feature_names, feature_names))
    feature_info <- data.frame(
      feature = feature_names,
      sd = as.numeric(feature_sd),
      low_variance = low_variance,
      null_score = as.numeric(null_scores),
      null_sigma2 = as.numeric(null_sigma2),
      stringsAsFactors = FALSE
    )
    return(list(
      S = S,
      distance = distance,
      metric = "smooth_fit",
      feature_info = feature_info,
      directional_score = directional_score,
      directional_delta = directional_delta,
      directional_success = directional_success,
      directional_lambda = directional_lambda,
      directional_sigma2 = directional_sigma2
    ))
  }

  directional_score <- matrix(
    rep(null_scores, each = d),
    nrow = d,
    ncol = d,
    dimnames = list(feature_names, feature_names)
  )
  directional_delta <- matrix(
    0,
    nrow = d,
    ncol = d,
    dimnames = list(feature_names, feature_names)
  )
  directional_success <- matrix(
    FALSE,
    nrow = d,
    ncol = d,
    dimnames = list(feature_names, feature_names)
  )
  directional_lambda <- matrix(
    lambda_default,
    nrow = d,
    ncol = d,
    dimnames = list(feature_names, feature_names)
  )
  directional_sigma2 <- matrix(
    rep(null_sigma2, each = d),
    nrow = d,
    ncol = d,
    dimnames = list(feature_names, feature_names)
  )

  for (j in seq_len(d)) {
    if (low_variance[j]) next
    for (k in seq_len(d)) {
      if (j == k || low_variance[k]) next
      # With missing data each pair uses its jointly observed samples; the
      # evidence gain is rescaled by n / m so pairs remain comparable per
      # observation (complete data: m = n, unchanged).
      rows <- which(observed[, j] & observed[, k])
      fit_jk <- tryCatch(
        .cavi_fit_directional_smooth_evidence(
          x_cov = X[rows, j],
          y = X[rows, k],
          measurement_sd = if (is.null(noise_info$measurement_sd)) NULL else noise_info$sd_mat[rows, k],
          K = K,
          rw_q = rw_q,
          ridge = ridge,
          lambda_sd_prior_rate = lambda_sd_prior_rate,
          lambda_mode = smooth_fit_lambda_mode,
          lambda_value = smooth_fit_lambda_value,
          lambda_min = lambda_min,
          lambda_max = lambda_max,
          sigma_min = sigma_min,
          sigma_max = sigma_max
        ),
        error = function(e) NULL
      )
      if (is.null(fit_jk)) {
        directional_success[j, k] <- FALSE
        next
      }
      directional_score[j, k] <- fit_jk$log_evidence
      directional_delta[j, k] <- fit_jk$delta * n / length(rows)
      directional_lambda[j, k] <- fit_jk$lambda
      directional_sigma2[j, k] <- fit_jk$sigma2
      directional_success[j, k] <- isTRUE(fit_jk$success)
    }
  }

  S_raw <- pmax(directional_delta, t(directional_delta))
  diag(S_raw) <- 0
  off_diag <- S_raw[row(S_raw) != col(S_raw)]
  max_off_diag <- if (length(off_diag)) max(off_diag, na.rm = TRUE) else 0
  if (!is.finite(max_off_diag) || max_off_diag <= 0) {
    S <- diag(d)
  } else {
    S <- S_raw / max_off_diag
    S <- pmin(pmax(S, 0), 1)
    diag(S) <- 1
  }

  dimnames(S) <- list(feature_names, feature_names)
  distance <- 1 - S
  dimnames(distance) <- dimnames(S)

  feature_info <- data.frame(
    feature = feature_names,
    sd = as.numeric(feature_sd),
    low_variance = low_variance,
    null_score = as.numeric(null_scores),
    null_sigma2 = as.numeric(null_sigma2),
    stringsAsFactors = FALSE
  )

  list(
    S = S,
    distance = distance,
    metric = "smooth_fit",
    feature_info = feature_info,
    directional_score = directional_score,
    directional_delta = directional_delta,
    directional_success = directional_success,
    directional_lambda = directional_lambda,
    directional_sigma2 = directional_sigma2
  )
}

# Center and standardize only the copy used to compute dimensionless similarity.
# Rescale centered columns before squaring so tiny/large finite units do not
# underflow/overflow the variance calculation. Constant columns remain zero.
# Missing entries stay NA; moments use each feature's observed values.
.cavi_standardize_similarity_features <- function(X) {
  centered <- sweep(X, 2L, colMeans(X, na.rm = TRUE), FUN = "-")
  column_scale <- apply(abs(centered), 2L, max, na.rm = TRUE)
  varying <- is.finite(column_scale) & column_scale > 0
  standardized <- matrix(0, nrow = nrow(X), ncol = ncol(X), dimnames = dimnames(X))
  rescaled_sd <- rep(0, ncol(X))
  if (any(varying)) {
    standardized[, varying] <- sweep(centered[, varying, drop = FALSE],
                                     2L, column_scale[varying], FUN = "/")
    rescaled_sd[varying] <- apply(standardized[, varying, drop = FALSE], 2L, stats::sd,
                                  na.rm = TRUE)
  }
  usable <- varying & is.finite(rescaled_sd) & rescaled_sd > 0
  if (any(usable)) {
    standardized[, usable] <- sweep(standardized[, usable, drop = FALSE],
                                    2L, rescaled_sd[usable], FUN = "/")
  }
  standardized[, !usable] <- 0
  feature_sd <- column_scale * rescaled_sd
  list(X = standardized, sd = feature_sd,
       standardized_sd = apply(standardized, 2L, stats::sd, na.rm = TRUE), usable = usable)
}

# Directional spline R^2 on jointly observed samples. For each predictor p and
# response j, samples observed in both are ordered by x_p and j is regressed on
# a natural spline of their rank position, as in the complete-data metric with
# n replaced by the pair's count m. Pairs with m < df + 3 receive 0.
.cavi_pairwise_spline_r2 <- function(X, low_variance, spline_r2_df) {
  observed <- !is.na(X)
  d <- ncol(X)
  bases <- list()
  basis_Q <- function(m) {
    key <- as.character(m)
    if (is.null(bases[[key]])) {
      position <- (seq_len(m) - 0.5) / m
      basis_qr <- qr(cbind(1, splines::ns(position, df = spline_r2_df, intercept = FALSE)))
      bases[[key]] <<- if (basis_qr$rank == spline_r2_df + 1L) qr.Q(basis_qr) else NA
    }
    bases[[key]]
  }
  directional <- matrix(0, d, d)
  for (predictor in which(!low_variance)) for (response in which(!low_variance)) {
    if (predictor == response) next
    rows <- which(observed[, predictor] & observed[, response])
    if (length(rows) < spline_r2_df + 3L) next
    Q <- basis_Q(length(rows))
    if (!is.matrix(Q)) next
    y <- X[rows[order(X[rows, predictor], method = "radix")], response]
    y <- y - mean(y)
    total <- sum(y^2)
    if (total <= 0) next
    directional[predictor, response] <- 1 - sum((y - Q %*% crossprod(Q, y))^2) / total
  }
  directional
}

.compute_same_ordering_similarity_spline_r2 <- function(X,
                                                        spline_r2_df = 5L,
                                                        min_feature_sd = 1e-8) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  spline_r2_df_raw <- spline_r2_df
  spline_r2_df <- as.integer(spline_r2_df)
  min_feature_sd <- as.numeric(min_feature_sd)[1]
  if (length(spline_r2_df_raw) != 1L || !is.finite(spline_r2_df_raw) ||
      is.na(spline_r2_df) || spline_r2_df_raw != spline_r2_df ||
      spline_r2_df < 1L || spline_r2_df > n - 2L) {
    stop("spline_r2_df must be a single integer between 1 and nrow(X) - 2.",
         call. = FALSE)
  }
  if (!is.finite(min_feature_sd) || min_feature_sd < 0) {
    stop("min_feature_sd must be a single finite nonnegative number.",
         call. = FALSE)
  }

  feature_names <- colnames(X) %||% paste0("V", seq_len(d))
  normalized <- .cavi_standardize_similarity_features(X)
  # The numerical SD guard is now applied in standardized units, so changing
  # the original feature units cannot discard a nonconstant column.
  low_variance <- !normalized$usable | normalized$standardized_sd < min_feature_sd
  rank_position <- (seq_len(n) - 0.5) / n
  spline_basis <- cbind(
    "(Intercept)" = 1,
    splines::ns(rank_position, df = spline_r2_df, intercept = FALSE)
  )
  basis_qr <- qr(spline_basis)
  if (basis_qr$rank != ncol(spline_basis)) {
    stop("The fixed-df spline basis is rank deficient.", call. = FALSE)
  }
  Q <- qr.Q(basis_qr)
  total_sum_squares <- colSums(normalized$X^2)
  # Each usable column has SST = n - 1; constant columns are excluded above.
  valid_response <- !low_variance
  directional_r2 <- matrix(
    0, nrow = d, ncol = d,
    dimnames = list(feature_names, feature_names)
  )

  if (anyNA(X)) {
    directional_r2[, ] <- .cavi_pairwise_spline_r2(X, low_variance, spline_r2_df)
  } else for (predictor in seq_len(d)) {
    if (low_variance[predictor]) next
    ordered_X <- normalized$X[order(X[, predictor], method = "radix"), , drop = FALSE]
    projected_coordinates <- crossprod(Q, ordered_X)
    residuals <- ordered_X - Q %*% projected_coordinates
    residual_sum_squares <- colSums(residuals^2)
    directional_r2[predictor, valid_response] <-
      1 - residual_sum_squares[valid_response] /
        total_sum_squares[valid_response]
  }
  directional_r2 <- pmin(pmax(directional_r2, 0), 1)
  similarity <- pmax(directional_r2, t(directional_r2))
  if (any(low_variance)) {
    similarity[low_variance, ] <- 0
    similarity[, low_variance] <- 0
  }
  diag(similarity) <- 1
  dimnames(similarity) <- list(feature_names, feature_names)
  distance <- 1 - similarity
  dimnames(distance) <- dimnames(similarity)

  list(
    S = similarity,
    distance = distance,
    metric = "spline_r2",
    normalization = "centered_unit_variance",
    feature_info = data.frame(
      feature = feature_names,
      sd = as.numeric(normalized$sd),
      standardized_sd = as.numeric(normalized$standardized_sd),
      low_variance = low_variance,
      stringsAsFactors = FALSE
    ),
    directional_r2 = directional_r2,
    spline_r2_df = spline_r2_df
  )
}

.compute_same_ordering_similarity <- function(X,
                                              S = NULL,
                                              metric = c("spearman", "pearson", "smooth_fit", "spline_r2"),
                                              use = "pairwise.complete.obs",
                                              abs_value = TRUE,
                                              min_feature_sd = 1e-8,
                                              spline_r2_df = 5L,
                                              K = NULL,
                                              rw_q = 2L,
                                              ridge = 0,
                                              lambda_sd_prior_rate = NULL,
                                              smooth_fit_lambda_mode = c("optimize", "fixed"),
                                              smooth_fit_lambda_value = 1,
                                              lambda_min = 1e-10,
                                              lambda_max = 1e10,
                                              sigma_min = 1e-10,
                                              sigma_max = 1e10,
                                              discretization = c("quantile", "equal", "kmeans")) {
  # With missing data, every metric uses each pair's jointly observed samples
  # (pairwise-complete); missing values are never filled for similarity.
  metric <- match.arg(metric)
  if (metric %in% c("spearman", "pearson")) {
    return(.compute_same_ordering_similarity_cor(
      X = X,
      metric = metric,
      use = use,
      abs_value = abs_value,
      min_feature_sd = min_feature_sd
    ))
  }

  if (identical(metric, "spline_r2")) {
    return(.compute_same_ordering_similarity_spline_r2(
      X = X,
      spline_r2_df = spline_r2_df,
      min_feature_sd = min_feature_sd
    ))
  }

  if (is.null(K)) {
    stop("K must be supplied when similarity_metric = 'smooth_fit'.")
  }
  .compute_same_ordering_similarity_smooth_fit(
    X = X,
    S = S,
    K = as.integer(K)[1],
    rw_q = rw_q,
    ridge = ridge,
    lambda_sd_prior_rate = lambda_sd_prior_rate,
    smooth_fit_lambda_mode = smooth_fit_lambda_mode,
    smooth_fit_lambda_value = smooth_fit_lambda_value,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    min_feature_sd = min_feature_sd,
    sigma_min = sigma_min,
    sigma_max = sigma_max
  )
}

.cavi_canonicalize_feature_clusters <- function(cluster_assign) {
  cluster_assign <- as.integer(cluster_assign)
  raw_ids <- sort(unique(cluster_assign))
  min_feature_idx <- vapply(raw_ids, function(id) {
    min(which(cluster_assign == id))
  }, integer(1))
  cluster_order <- raw_ids[order(min_feature_idx)]
  canonical_assign <- match(cluster_assign, cluster_order)
  cluster_members <- split(seq_along(canonical_assign), canonical_assign)

  list(
    feature_cluster = canonical_assign,
    cluster_order = cluster_order,
    cluster_members = cluster_members,
    cluster_sizes = vapply(cluster_members, length, integer(1))
  )
}

.cavi_validate_positive_integer <- function(x, name) {
  x_raw <- x
  x <- suppressWarnings(as.integer(x))
  if (!is.numeric(x_raw) || length(x_raw) != 1L || !is.finite(x_raw) ||
      is.na(x) || x_raw != x || x < 1L) {
    stop(name, " must be a single positive integer.", call. = FALSE)
  }
  x
}

.cavi_validate_similarity_distance <- function(distance) {
  distance <- as.matrix(distance)
  if (!is.numeric(distance) || length(dim(distance)) != 2L ||
      nrow(distance) != ncol(distance) || nrow(distance) < 1L) {
    stop("distance must be a nonempty square numeric matrix.", call. = FALSE)
  }
  if (any(!is.finite(distance))) {
    stop("distance must contain only finite values.", call. = FALSE)
  }
  tolerance <- 1e-10
  if (any(distance < -tolerance) ||
      max(abs(distance - t(distance))) > tolerance ||
      max(abs(diag(distance))) > tolerance) {
    stop(
      "distance must be symmetric, nonnegative, and have a zero diagonal.",
      call. = FALSE
    )
  }
  distance <- pmax(distance, 0)
  diag(distance) <- 0
  distance
}

.cavi_mean_silhouette <- function(distance, cluster_assign) {
  distance <- .cavi_validate_similarity_distance(distance)
  cluster_assign <- as.integer(cluster_assign)
  d <- nrow(distance)
  if (length(cluster_assign) != d || anyNA(cluster_assign)) {
    stop("cluster_assign must contain one finite cluster label per feature.",
         call. = FALSE)
  }
  cluster_ids <- sort(unique(cluster_assign))
  if (length(cluster_ids) < 2L) {
    return(NA_real_)
  }

  silhouette <- numeric(d)
  for (j in seq_len(d)) {
    own_members <- which(cluster_assign == cluster_assign[j])
    if (length(own_members) == 1L) {
      silhouette[j] <- 0
      next
    }
    within_distance <- mean(distance[j, own_members[own_members != j]])
    other_clusters <- cluster_ids[cluster_ids != cluster_assign[j]]
    between_distance <- min(vapply(
      other_clusters,
      function(cluster_id) {
        mean(distance[j, cluster_assign == cluster_id])
      },
      numeric(1)
    ))
    denominator <- max(within_distance, between_distance)
    silhouette[j] <- if (denominator > 0) {
      (between_distance - within_distance) / denominator
    } else {
      0
    }
  }
  mean(silhouette)
}

.cavi_select_similarity_dimension <- function(
    distance,
    cluster_linkage = "single",
    max_intrinsic_dim = 8L,
    min_cluster_size = 2L) {
  distance <- .cavi_validate_similarity_distance(distance)
  cluster_linkage <- .cavi_validate_cluster_linkage(cluster_linkage)
  max_intrinsic_dim <- .cavi_validate_positive_integer(
    max_intrinsic_dim, "max_intrinsic_dim"
  )
  min_cluster_size <- .cavi_validate_positive_integer(
    min_cluster_size, "min_cluster_size"
  )
  d <- nrow(distance)
  feature_names <- rownames(distance) %||% paste0("V", seq_len(d))
  candidate_max <- min(max_intrinsic_dim, d - 1L)
  candidate_M <- if (candidate_max >= 2L) {
    seq.int(2L, candidate_max)
  } else {
    integer(0)
  }
  hc <- if (d >= 2L) {
    stats::hclust(stats::as.dist(distance), method = cluster_linkage)
  } else {
    NULL
  }

  cluster_assignments <- vector("list", length(candidate_M))
  cluster_sizes <- vector("list", length(candidate_M))
  mean_silhouette <- rep(NA_real_, length(candidate_M))
  minimum_cluster_size <- rep(NA_integer_, length(candidate_M))
  maximum_cluster_size <- rep(NA_integer_, length(candidate_M))
  feasible <- rep(FALSE, length(candidate_M))

  for (idx in seq_along(candidate_M)) {
    raw_cluster <- stats::cutree(hc, k = candidate_M[idx])
    cluster_info <- .cavi_canonicalize_feature_clusters(raw_cluster)
    cluster_assignments[[idx]] <- cluster_info$feature_cluster
    cluster_sizes[[idx]] <- cluster_info$cluster_sizes
    mean_silhouette[idx] <- .cavi_mean_silhouette(
      distance, cluster_info$feature_cluster
    )
    minimum_cluster_size[idx] <- min(cluster_info$cluster_sizes)
    maximum_cluster_size[idx] <- max(cluster_info$cluster_sizes)
    feasible[idx] <- minimum_cluster_size[idx] >= min_cluster_size
  }

  diagnostics <- data.frame(
    M = candidate_M,
    mean_silhouette = mean_silhouette,
    minimum_cluster_size = minimum_cluster_size,
    maximum_cluster_size = maximum_cluster_size,
    feasible = feasible,
    stringsAsFactors = FALSE
  )
  diagnostics$cluster_sizes <- I(cluster_sizes)

  eligible_idx <- which(feasible & is.finite(mean_silhouette))
  if (length(eligible_idx)) {
    best_idx <- eligible_idx[which.max(mean_silhouette[eligible_idx])]
    selected_M <- candidate_M[best_idx]
    selected_cluster <- cluster_assignments[[best_idx]]
    selected_sizes <- cluster_sizes[[best_idx]]
    fallback_reason <- NULL
  } else {
    selected_M <- 1L
    selected_cluster <- rep(1L, d)
    selected_sizes <- stats::setNames(d, "ordering_1")
    fallback_reason <- paste0(
      "No cut with M >= 2 had every cluster size at least ",
      min_cluster_size, "; used one ordering."
    )
  }

  selected_info <- .cavi_canonicalize_feature_clusters(selected_cluster)
  selected_cluster <- selected_info$feature_cluster
  selected_sizes <- selected_info$cluster_sizes
  ordering_labels <- .cavi_partition_order_labels(selected_M)
  names(selected_sizes) <- ordering_labels
  names(selected_cluster) <- feature_names
  cluster_proportions <- selected_sizes / d

  list(
    selected_M = as.integer(selected_M),
    feature_cluster = selected_cluster,
    cluster_sizes = selected_sizes,
    cluster_proportions = cluster_proportions,
    criterion = "maximum_mean_silhouette_among_feasible_cuts",
    fallback_reason = fallback_reason,
    max_intrinsic_dim = as.integer(candidate_max),
    requested_max_intrinsic_dim = max_intrinsic_dim,
    min_cluster_size = min_cluster_size,
    cluster_linkage = cluster_linkage,
    diagnostics = diagnostics,
    hclust = hc
  )
}

.cavi_resolve_similarity_methods <- function(methods = NULL,
                                             pca_components = NULL,
                                             M) {
  valid_methods <- c("PCA", "fiedler", "pcurve", "tSNE", "random", "isomap")

  if (is.null(methods)) {
    methods <- rep("PCA", M)
  } else {
    methods <- as.character(methods)
    if (length(methods) == 1L) {
      methods <- rep(methods, M)
    } else if (length(methods) != M) {
      stop(sprintf("init_methods must have length 1 or M=%d under similarity initialization.", M))
    }
    methods <- vapply(methods, match.arg, character(1), choices = valid_methods)
  }

  if (is.null(pca_components)) {
    # Similarity initialization fits each ordering on a different feature
    # subset. A repeated PCA method therefore means PC1 within every subset,
    # rather than successive PCs across otherwise unrelated feature blocks.
    pca_components <- .cavi_default_subset_pca_components(methods)
  } else {
    if (!is.numeric(pca_components) || !length(pca_components)) {
      stop("pca_components must contain positive integers.", call. = FALSE)
    }
    pca_components <- vapply(pca_components, .mpcurve_integer_setting,
      integer(1), name = "pca_components")
    if (length(pca_components) == 1L) {
      pca_components <- rep(pca_components, M)
    } else if (length(pca_components) != M) {
      stop(sprintf("pca_components must have length 1 or M=%d under similarity initialization.", M))
    }
  }

  bad_pca <- methods == "PCA" & (is.na(pca_components) | pca_components < 1L)
  if (any(bad_pca)) {
    stop("pca_components for PCA-based methods must contain positive integers.")
  }

  list(
    methods = methods,
    pca_components = pca_components
  )
}

.cavi_similarity_subset_fit <- function(X_sub,
                                        S = NULL,
                                        K,
                                        method,
                                        pca_component = 1L,
                                        rw_q = 2L,
                                        ridge = 0,
                                        lambda_sd_prior_rate = NULL,
                                        lambda_min = 1e-10,
                                        lambda_max = 1e10,
                                        sigma_min = 1e-10,
                                        sigma_max = 1e10,
                                        max_iter = 5L,
                                        tol = 1e-6,
                                        discretization = c("quantile", "equal", "kmeans"),
                                        cluster_label = NULL,
                                        method_args = list(),
                                        init_on_failure = "pca",
                                        verbose = FALSE) {
  X_sub <- as.matrix(X_sub)
  method <- match.arg(method, c("PCA", "fiedler", "pcurve", "tSNE", "random", "isomap"))
  fit <- .cavi_fit_from_method(
    X = X_sub, S = S, K = K, method = method,
    pca_component = pca_component, method_args = method_args,
    init_on_failure = init_on_failure, ordering_label = cluster_label,
    rw_q = rw_q, ridge = ridge, lambda_sd_prior_rate = lambda_sd_prior_rate,
    lambda_min = lambda_min, lambda_max = lambda_max,
    sigma_min = sigma_min, sigma_max = sigma_max,
    max_iter = max_iter, tol = tol, discretization = discretization,
    strict_K = TRUE, verbose = verbose
  )
  info <- fit$init_info %||% list(method_requested = method, method_used = method,
    pca_component_requested = NA_integer_, pca_component_used = NA_integer_,
    fallback = FALSE, fallback_reason = NULL)
  c(list(fit = fit), info)
}

.cavi_init_m_trajectories_similarity <- function(X,
                                                 S = NULL,
                                                 M = 2L,
                                                 methods = NULL,
                                                 method_args = list(),
                                                 init_on_failure = "pca",
                                                 pca_components = NULL,
                                                 K = NULL,
                                                 rw_q = 2L,
                                                 ridge = 0,
                                                 lambda_sd_prior_rate = NULL,
                                                 smooth_fit_lambda_mode = c("optimize", "fixed"),
                                                 smooth_fit_lambda_value = 1,
                                                 lambda_min = 1e-10,
                                                 lambda_max = 1e10,
                                                 sigma_min = 1e-10,
                                                 sigma_max = 1e10,
                                                 discretization = c("quantile", "equal", "kmeans"),
                                                 num_iter = 0L,
                                                 similarity_metric = c("spearman", "pearson", "smooth_fit", "spline_r2"),
                                                 cluster_linkage = "single",
                                                 similarity_min_feature_sd = 1e-8,
                                                 spline_r2_df = 5L,
                                                 similarity_precomputed = NULL,
                                                 verbose = FALSE) {
  X <- as.matrix(X)
  M <- as.integer(M)
  d <- ncol(X)
  if (M > d) {
    stop(sprintf("similarity initialization requires M=%d to be <= ncol(X)=%d.", M, d))
  }

  discretization <- match.arg(discretization)
  similarity_metric <- match.arg(similarity_metric)
  smooth_fit_lambda_mode <- match.arg(smooth_fit_lambda_mode)
  cluster_linkage <- .cavi_validate_cluster_linkage(cluster_linkage)
  K_use <- .cavi_resolve_K(X, K, rw_q)
  if (identical(similarity_metric, "smooth_fit") && isTRUE(M > 1L) &&
      is.finite(ridge) && ridge <= 0) {
    warning(
      paste0(
        "similarity_metric = 'smooth_fit' with ridge = 0 uses an intrinsic RW prior. ",
        "This is fine for pairwise initialization, but the directional scores are pseudo-evidence ",
        "rather than fully proper marginal likelihoods; use a small positive ridge if you want ",
        "the smoother evidence to be theoretically proper."
      ),
      call. = FALSE
    )
  }
  method_info <- .cavi_resolve_similarity_methods(
    methods = methods,
    pca_components = pca_components,
    M = M
  )

  if (is.null(similarity_precomputed)) {
    similarity <- .compute_same_ordering_similarity(
      X = X,
      S = S,
      metric = similarity_metric,
      min_feature_sd = similarity_min_feature_sd,
      spline_r2_df = spline_r2_df,
      K = K_use,
      rw_q = rw_q,
      ridge = ridge,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      smooth_fit_lambda_mode = smooth_fit_lambda_mode,
      smooth_fit_lambda_value = smooth_fit_lambda_value,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      discretization = discretization
    )
  } else {
    similarity <- similarity_precomputed
    if (!is.list(similarity) || is.null(similarity$S) ||
        is.null(similarity$distance) ||
        !identical(as.character(similarity$metric)[1], similarity_metric)) {
      stop(
        "similarity_precomputed must contain S and distance for similarity_metric.",
        call. = FALSE
      )
    }
    similarity$distance <- .cavi_validate_similarity_distance(
      similarity$distance
    )
    similarity$S <- as.matrix(similarity$S)
    if (!all(dim(similarity$S) == c(d, d)) ||
        !all(dim(similarity$distance) == c(d, d)) ||
        any(!is.finite(similarity$S))) {
      stop("similarity_precomputed matrices must be finite d x d matrices.",
           call. = FALSE)
    }
  }
  hc <- stats::hclust(stats::as.dist(similarity$distance), method = cluster_linkage)
  raw_cluster <- stats::cutree(hc, k = M)
  cluster_info <- .cavi_canonicalize_feature_clusters(raw_cluster)
  ord_labels <- .cavi_partition_order_labels(M)
  # A one-feature cluster initializes an ordering from a single feature. Such
  # orderings are weakly identified, and with missing values their positions
  # can collapse; report it rather than change the requested cut.
  singletons <- which(cluster_info$cluster_sizes == 1L)
  if (length(singletons)) {
    names_X <- colnames(X) %||% paste0("V", seq_len(ncol(X)))
    message(sprintf(paste0(
      "Similarity initialization placed %s alone in %s (cluster sizes %s). A ",
      "one-feature ordering is weakly identified, especially with missing values; ",
      "consider a smaller intrinsic_dim or another similarity_metric."),
      paste(sprintf("feature '%s'", names_X[unlist(cluster_info$cluster_members[singletons])]),
            collapse = ", "),
      if (length(singletons) == 1L) sprintf("ordering %s", ord_labels[singletons]) else
        sprintf("orderings %s", paste(ord_labels[singletons], collapse = ", ")),
      paste(cluster_info$cluster_sizes, collapse = "+")))
  }

  fits <- vector("list", M)
  init_info <- vector("list", M)

  for (m in seq_len(M)) {
    cols_m <- cluster_info$cluster_members[[m]]
    subset_res <- .cavi_similarity_subset_fit(
      X_sub = X[, cols_m, drop = FALSE],
      S = .cavi_subset_measurement_sd(S, cols_m),
      K = K_use,
      method = method_info$methods[m],
      method_args = method_args,
      init_on_failure = init_on_failure,
      pca_component = method_info$pca_components[m],
      rw_q = rw_q,
      ridge = ridge,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      max_iter = num_iter,
      tol = 1e-6,
      discretization = discretization,
      cluster_label = sprintf("cluster %s", ord_labels[m]),
      verbose = FALSE
    )

    # Use max_iter = 0 to preserve the subset-specialized cell ordering R.
    # Iterating on all d features without feature weights would corrupt the
    # ordering (especially with known S where all features contribute equally).
    # The partition CAVI loop (with feature weights) handles refinement.
    fits[[m]] <- cavi(
      X = X,
      convergence = "relative",
      K = length(subset_res$fit$params$pi),
      responsibilities_init = subset_res$fit$gamma,
      position_prior_init = colMeans(subset_res$fit$gamma),
      S = S,
      rw_q = rw_q,
      ridge = ridge,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      max_iter = 0L,
      tol = 1e-6,
      verbose = FALSE
    )

    init_info[[m]] <- list(
      label = ord_labels[m],
      method_requested = subset_res$method_requested,
      method_used = subset_res$method_used,
      pca_component_requested = subset_res$pca_component_requested,
      pca_component_used = subset_res$pca_component_used,
      K = K_use,
      discretization = discretization,
      feature_idx = cols_m,
      cluster_size = length(cols_m),
      cluster_min_feature = min(cols_m),
      fallback = subset_res$fallback,
      fallback_reason = subset_res$fallback_reason
    )
    if (!is.null(subset_res$k_used)) init_info[[m]]$k_used <- subset_res$k_used
  }

  names(init_info) <- ord_labels
  names(cluster_info$cluster_sizes) <- ord_labels
  names(cluster_info$cluster_members) <- ord_labels
  feature_cluster_named <- cluster_info$feature_cluster
  names(feature_cluster_named) <- colnames(X) %||% paste0("V", seq_len(ncol(X)))

  similarity_init <- list(
    partition_init = "similarity",
    similarity_metric = similarity_metric,
    similarity_normalization = similarity$normalization %||% NULL,
    cluster_linkage = cluster_linkage,
    similarity_min_feature_sd = as.numeric(similarity_min_feature_sd)[1],
    spline_r2_df = if (identical(similarity_metric, "spline_r2")) {
      similarity$spline_r2_df %||% as.integer(spline_r2_df)[1]
    } else NULL,
    smooth_fit_lambda_mode = smooth_fit_lambda_mode,
    smooth_fit_lambda_value = as.numeric(smooth_fit_lambda_value)[1],
    S = similarity$S,
    distance = similarity$distance,
    feature_cluster = feature_cluster_named,
    cluster_order = cluster_info$cluster_order,
    cluster_members = cluster_info$cluster_members,
    cluster_sizes = cluster_info$cluster_sizes,
    feature_info = similarity$feature_info,
    init_methods = method_info$methods,
    pca_components = method_info$pca_components,
    hclust_order = hc$order,
    directional_r2 = similarity$directional_r2 %||% NULL,
    directional_score = similarity$directional_score %||% NULL,
    directional_delta = similarity$directional_delta %||% NULL,
    directional_success = similarity$directional_success %||% NULL,
    directional_lambda = similarity$directional_lambda %||% NULL,
    directional_sigma2 = similarity$directional_sigma2 %||% NULL
  )

  if (isTRUE(verbose)) {
    cat(sprintf(
      "[init_m_cavi similarity] M=%d metric=%s linkage=%s sizes=%s\n",
      M,
      similarity_metric,
      cluster_linkage,
      paste(cluster_info$cluster_sizes, collapse = ",")
    ))
  }

  list(
    fits = fits,
    init_info = init_info,
    ordering_similarity = NULL,
    similarity_init = similarity_init
  )
}


.validate_partition_assignment_controls <- function(partition_prior = c("adaptive", "fixed"),
                                                    partition_prior_init = NULL,
                                                    M,
                                                    caller = "soft_partition_cavi()") {
  partition_prior <- match.arg(partition_prior)

  if (!is.null(partition_prior_init)) {
    init <- as.numeric(partition_prior_init)
    if (length(init) != M || any(!is.finite(init)) || any(init < 0) || sum(init) <= 0) {
      stop(
        caller,
        ": `partition_prior_init` must be a nonnegative vector of length M with positive sum.",
        call. = FALSE
      )
    }
    partition_prior_init <- init / sum(init)
  }

  list(
    partition_prior = partition_prior,
    partition_prior_init = partition_prior_init
  )
}

# Ordering-assignment terms of the objective: the partition prior omega (fixed
# or estimated by empirical Bayes from the soft feature assignments), the
# expected log assignment prior, and the assignment entropy at temperature T_now.
.cavi_partition_assignment_info <- function(weights,
                                            T_now = 1,
                                            partition_prior = c("adaptive", "fixed"),
                                            partition_prior_init = NULL,
                                            caller = ".cavi_partition_assignment_info()") {
  weights <- as.matrix(weights)
  M <- ncol(weights)
  d <- nrow(weights)
  ctl <- .validate_partition_assignment_controls(
    partition_prior = partition_prior,
    partition_prior_init = partition_prior_init,
    M = M,
    caller = caller
  )
  z_entropy <- -sum(weights * log(weights + 1e-300))

  if (identical(ctl$partition_prior, "fixed")) {
    omega <- if (is.null(ctl$partition_prior_init)) {
      rep(1 / M, M)
    } else {
      ctl$partition_prior_init / sum(ctl$partition_prior_init)
    }
  } else {
    omega <- colSums(weights)
    omega <- unname(omega / sum(omega))
  }
  e_log_omega <- log(pmax(omega, 1e-300))
  expected_log_assign <- sum(weights * rep(e_log_omega, each = d))
  list(
    partition_prior = ctl$partition_prior,
    omega = omega,
    e_log_omega = e_log_omega,
    expected_log_assign = as.numeric(expected_log_assign),
    z_entropy = as.numeric(z_entropy),
    objective = as.numeric(expected_log_assign + T_now * z_entropy)
  )
}

.cavi_hard_gamma_from_cluster_rank <- function(cluster_rank, n, K) {
  gamma <- matrix(0, nrow = n, ncol = K)
  gamma[cbind(seq_len(length(cluster_rank)), cluster_rank)] <- 1
  gamma
}

.cavi_exact_k_init_from_ordering <- function(X,
                                             ordering_vec,
                                             K,
                                             discretization = c("quantile", "equal", "kmeans"),
                                             strict_K = TRUE) {
  X <- .mpcurve_initialization_data(X)
  discretization <- match.arg(discretization)

  capture_init <- function(disc) {
    warnings <- character(0)
    init <- withCallingHandlers(
      .cavi_bin_ordering(
        X = X,
        ordering_vec = ordering_vec,
        K = K,
        discretization = disc,
        na_action = "drop",
        eps = 1e-12
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(init = init, warnings = unique(warnings))
  }

  primary <- capture_init(discretization)
  K_primary <- length(primary$init$pi)
  if (!strict_K || K_primary == K) {
    attr(primary$init, "discretization_used") <- discretization
    attr(primary$init, "suppressed_warnings") <- primary$warnings
    return(primary$init)
  }

  if (discretization != "equal") {
    fallback <- capture_init("equal")
    if (length(fallback$init$pi) == K) {
      attr(fallback$init, "discretization_used") <- "equal"
      attr(fallback$init, "suppressed_warnings") <- c(primary$warnings, fallback$warnings)
      attr(fallback$init, "fallback_from") <- discretization
      attr(fallback$init, "fallback_reason") <- sprintf(
        "Requested discretization '%s' produced K=%d instead of K=%d.",
        discretization, K_primary, K
      )
      return(fallback$init)
    }
  }

  stop(
    "Unable to construct a partition initialization with the requested K=",
    K, ". For intrinsic_dim > 1, all orderings must use the same K."
  )
}

.cavi_get_ordering_result <- function(X, method, pca_component = 1L,
                                      method_args = list(), init_on_failure = "pca",
                                      ordering_label = NULL) {
  X <- .mpcurve_initialization_data(X)
  method <- match.arg(method, c("PCA", "fiedler", "pcurve", "tSNE", "random", "isomap"))
  init_on_failure <- match.arg(init_on_failure, c("pca", "error"))
  if (method == "PCA" && !"component" %in% names(method_args)) {
    method_args$component <- pca_component
  }
  .mpcurve_validate_ordering_args(X, method, method_args)
  helper <- switch(method, PCA = PCA_ordering, fiedler = fiedler_ordering,
    pcurve = pcurve_ordering, tSNE = tSNE_ordering, isomap = isomap_ordering,
    random = .random_ordering)
  # Catch only ordering computation failures; invalid settings and downstream
  # model-fitting errors must not change the requested initializer.
  result <- tryCatch(do.call(helper, c(list(X = X), method_args)), error = identity)
  info <- list(method_requested = method, method_used = method,
    pca_component_requested = if (method == "PCA") method_args$component else NA_integer_,
    pca_component_used = if (method == "PCA") method_args$component else NA_integer_,
    fallback = FALSE, fallback_reason = NULL)
  if (inherits(result, "error")) {
    reason <- conditionMessage(result)
    prefix <- if (is.null(ordering_label)) "" else paste0(ordering_label, ": ")
    if (init_on_failure == "error" || method == "PCA") {
      stop(prefix, method, " initialization failed: ", reason, call. = FALSE)
    }
    warning(prefix, method, " initialization failed: ", reason,
      ". Falling back to PCA component 1. Set init_control = ",
      "mpcurve_init_control(on_failure = 'error') to stop instead.", call. = FALSE)
    result <- tryCatch(PCA_ordering(X, component = 1L), error = identity)
    if (inherits(result, "error")) {
      stop(prefix, "PCA fallback failed: ", conditionMessage(result),
        ". Original ", method, " initialization failure: ", reason, call. = FALSE)
    }
    info$method_used <- "PCA"
    info$pca_component_used <- 1L
    info$fallback <- TRUE
    info$fallback_reason <- reason
  }
  if (!is.null(result$k_used)) info$k_used <- result$k_used
  attr(result, "init_info") <- info
  result
}

.cavi_build_from_ordering <- function(X,
                                      ordering_vec,
                                      S = NULL,
                                      K,
                                      rw_q = 2L,
                                      ridge = 0,
                                      lambda_init = 1,
                                      lambda_sd_prior_rate = NULL,
                                      lambda_min = 1e-10,
                                      lambda_max = 1e10,
                                      sigma_min = 1e-10,
                                      sigma_max = 1e10,
                                      max_iter = 5L,
                                      tol = 1e-6,
                                      discretization = c("quantile", "equal", "kmeans"),
                                      strict_K = FALSE,
                                      ordering_label = NULL,
                                      verbose = FALSE) {
  discretization <- match.arg(discretization)
  init <- .cavi_exact_k_init_from_ordering(
    X = X,
    ordering_vec = ordering_vec,
    K = K,
    discretization = discretization,
    strict_K = strict_K
  )
  disc_used <- attr(init, "discretization_used") %||% discretization
  fallback_from <- attr(init, "fallback_from")
  if (!is.null(fallback_from)) {
    prefix <- if (is.null(ordering_label) || !nzchar(ordering_label)) "" else paste0(ordering_label, ": ")
    warning(
      prefix,
      attr(init, "fallback_reason"),
      " Falling back to equal-width bins so all orderings use the same K.",
      call. = FALSE
    )
  }
  gamma0 <- .cavi_hard_gamma_from_cluster_rank(init$cluster_rank, nrow(X), length(init$pi))
  if (anyNA(X)) gamma0 <- .mpcurve_missing_initial_responsibilities(gamma0)
  cavi_args <- list(
    X = X,
    K = length(init$pi),
    responsibilities_init = gamma0,
    position_prior_init = init$pi,
    S = S,
    lambda_init = rep(lambda_init, ncol(X)),
    lambda_sd_prior_rate = lambda_sd_prior_rate,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    sigma_min = sigma_min,
    sigma_max = sigma_max,
    rw_q = rw_q,
    ridge = ridge,
    discretization = disc_used,
    max_iter = max_iter,
    tol = tol,
    # Preserve initialization across changes to the outer stopping rule.
    convergence = "relative",
    verbose = verbose
  )
  if (is.null(S)) {
    cavi_args$sigma2_init <- if (anyNA(X)) {
      pmin(pmax(.mpcurve_observed_variance(X), sigma_min), sigma_max)
    } else init$sigma2
  }
  do.call(cavi, cavi_args)
}


.cavi_chol_with_jitter <- function(A,
                                   base_jitter = 1e-8,
                                   max_tries = 6L) {
  A <- 0.5 * (A + t(A))
  L <- tryCatch(chol(A), error = function(e) NULL)
  if (!is.null(L)) {
    return(list(chol = L, matrix = A, jitter = 0))
  }

  jitter0 <- max(base_jitter, 1e-10 * max(1, mean(diag(A))))
  for (i in seq_len(max_tries)) {
    jitter_i <- jitter0 * (10 ^ (i - 1L))
    A_try <- A + diag(jitter_i, nrow(A))
    L <- tryCatch(chol(A_try), error = function(e) NULL)
    if (!is.null(L)) {
      return(list(chol = L, matrix = A_try, jitter = jitter_i))
    }
  }

  stop("Unable to stabilize the weighted q(U) precision matrix for Cholesky factorization.")
}


.cavi_update_lambda_weighted <- function(q_u, r_rank,
                                         lambda_sd_prior_rate = NULL,
                                         lambda_min = 1e-10,
                                         lambda_max = 1e10,
                                         lambda_current = NULL,
                                         feature_weights = NULL,
                                         feature_active = NULL) {
  d <- nrow(q_u$m_mat)
  w <- .cavi_resolve_feature_weights(feature_weights, d)
  if (is.null(feature_active)) {
    feature_active <- rep(TRUE, d)
  }
  feature_active <- as.logical(feature_active)
  if (length(feature_active) != d) {
    stop("feature_active must have length d.")
  }
  denom <- pmax(q_u$eq_quad, 1e-12)
  lambda_sd_prior_rate_val <- .lambda_sd_prior_rate_value(lambda_sd_prior_rate)
  if (lambda_sd_prior_rate_val <= 0) {
    lambda_new <- pmin(pmax(r_rank / denom, lambda_min), lambda_max)
  } else {
    lambda_new <- .optimize_lambda_induced_exp(
      eq_quad = denom,
      r_rank = r_rank,
      rate = lambda_sd_prior_rate_val,
      lambda_min = lambda_min,
      lambda_max = lambda_max
    )
  }
  if (!is.null(lambda_current)) {
    hold <- !feature_active | w <= .cavi_partition_weight_floor()
    lambda_new[hold] <- lambda_current[hold]
  }
  lambda_new
}

.cavi_entropy_u <- function(q_u, K) {
  sum(0.5 * (K * log(2 * base::pi * exp(1)) + q_u$logdetS))
}

.cavi_entropy_u_terms <- function(q_u, K) {
  0.5 * (K * log(2 * base::pi * exp(1)) + q_u$logdetS)
}

.cavi_prior_terms_from_state <- function(q_u, lambda_vec, Q_K, rw_q = 2L,
                                         lambda_sd_prior_rate = NULL) {
  prior_meta <- .rw_precision_metadata(Q_K, rw_q = rw_q)
  r_rank <- prior_meta$rank
  logdet_Q <- prior_meta$logdet
  0.5 * (r_rank * (log(lambda_vec) - log(2 * base::pi)) + logdet_Q) -
    0.5 * lambda_vec * q_u$eq_quad +
    .lambda_sd_prior_terms(lambda_vec, rate = lambda_sd_prior_rate, include_constant = TRUE)
}

.cavi_cell_terms_from_state <- function(R, pi_vec) {
  log_pi_term <- sum(R * rep(log(pmax(pi_vec, .Machine$double.eps)), each = nrow(R)))
  entropy_c <- -sum(R * log(R + 1e-300))
  as.numeric(log_pi_term + entropy_c)
}

.cavi_posterior_state_from_fit <- function(fit) {
  q_u <- list(
    m_mat = fit$posterior$mean,
    sdiag_mat = fit$posterior$var,
    S_list = fit$posterior$cov,
    logdetS = vapply(
      fit$posterior$cov,
      function(S) as.numeric(determinant(S, logarithm = TRUE)$modulus),
      numeric(1)
    )
  )
  q_u$eq_quad <- vapply(seq_len(nrow(q_u$m_mat)), function(j) {
    mj <- q_u$m_mat[j, ]
    Sj <- q_u$S_list[[j]]
    as.numeric(crossprod(mj, fit$Q_K %*% mj) + sum(fit$Q_K * Sj))
  }, numeric(1))
  q_u
}


.cavi_feature_scores_from_state <- function(X, R, q_u, sigma2, lambda_vec,
                                            measurement_sd = NULL,
                                            Q_K, rw_q = 2L,
                                            lambda_sd_prior_rate = NULL,
                                            include_prior = TRUE,
                                            feature_active = NULL) {
  X <- as.matrix(X)
  R <- as.matrix(R)
  n <- nrow(X)
  d <- ncol(X)
  K <- ncol(R)
  if (is.null(feature_active)) {
    feature_active <- rep(TRUE, d)
  }
  feature_active <- as.logical(feature_active)
  if (length(feature_active) != d) {
    stop("feature_active must have length d.")
  }
  prior_meta <- .rw_precision_metadata(Q_K, rw_q = rw_q)
  r_rank <- prior_meta$rank
  logdet_Q <- prior_meta$logdet
  noise_info <- .cavi_resolve_measurement_sd(
    S = measurement_sd,
    X = X,
    caller = ".cavi_feature_scores_from_state()"
  )
  if (anyNA(X)) {
    scores <- .cavi_observed_feature_scores(X, R, q_u, sigma2, noise_info$var_mat)
    if (isTRUE(include_prior)) scores <- scores + .cavi_prior_terms_from_state(
      q_u, lambda_vec, Q_K, rw_q, lambda_sd_prior_rate)
    scores[!feature_active] <- 0
    return(scores)
  }

  scores <- numeric(d)
  m_mat <- q_u$m_mat
  sdiag_mat <- q_u$sdiag_mat

  for (j in seq_len(d)) {
    if (!feature_active[j]) {
      scores[j] <- 0
      next
    }
    diff_j <- matrix(X[, j], nrow = n, ncol = K) -
      matrix(m_mat[j, ], nrow = n, ncol = K, byrow = TRUE)
    quad_jk <- colSums(R * (diff_j^2))
    unc_jk <- colSums(R) * sdiag_mat[j, ]
    if (is.null(noise_info$measurement_sd)) {
      data_j <- -0.5 * sum(
        colSums(R) * log(2 * base::pi * sigma2[j]) +
          (quad_jk + unc_jk) / sigma2[j]
      )
    } else {
      var_j <- noise_info$var_mat[, j]
      const_j <- sum(log(2 * base::pi * var_j))
      unc_row_j <- rowSums(R * rep(sdiag_mat[j, ], each = n))
      data_j <- -0.5 * (
        const_j +
          sum(rowSums(R * (diff_j^2) / var_j)) +
          sum(unc_row_j / var_j)
      )
    }

    prior_j <- 0
    if (isTRUE(include_prior)) {
      prior_j <- 0.5 * (r_rank * (log(lambda_vec[j]) - log(2 * base::pi)) + logdet_Q) -
        0.5 * lambda_vec[j] * q_u$eq_quad[j] +
        .lambda_sd_prior_terms(
          lambda_vec = lambda_vec[j],
          rate = lambda_sd_prior_rate,
          include_constant = TRUE
        )
    }
    scores[j] <- data_j + prior_j
  }

  scores
}


# ============================================================
# Per-ordering starting fits
# ============================================================

.cavi_fit_from_method <- function(X,
                                  S = NULL,
                                  method,
                                  K,
                                  rw_q = 2L,
                                  ridge = 0,
                                  lambda_sd_prior_rate = NULL,
                                  lambda_min = 1e-10,
                                  lambda_max = 1e10,
                                  sigma_min = 1e-10,
                                  sigma_max = 1e10,
                                  max_iter = 5L,
                                  tol = 1e-6,
                                  discretization = c("quantile", "equal", "kmeans"),
                                  pca_component = 1L,
                                  strict_K = FALSE,
                                  method_args = list(),
                                  init_on_failure = "pca",
                                  ordering_label = NULL,
                                  verbose = FALSE) {
  method <- match.arg(method, c("PCA", "fiedler", "pcurve", "tSNE", "random", "isomap"))
  discretization <- match.arg(discretization)

  ordering_result <- .cavi_get_ordering_result(X, method = method, pca_component = pca_component, method_args = method_args,
    init_on_failure = init_on_failure, ordering_label = ordering_label)
  fit <- .cavi_build_from_ordering(
    X = X,
    ordering_vec = ordering_result$t,
    S = S,
    K = K,
    rw_q = rw_q,
    ridge = ridge,
    lambda_init = 1,
    lambda_sd_prior_rate = lambda_sd_prior_rate,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    sigma_min = sigma_min,
    sigma_max = sigma_max,
    max_iter = max_iter,
    tol = tol,
    discretization = discretization,
    strict_K = strict_K,
    ordering_label = ordering_label %||% if (method == "PCA") {
      sprintf("%s(PC%d)", method, pca_component)
    } else {
      method
    },
    verbose = verbose
  )
  fit$init_info <- attr(ordering_result, "init_info")
  fit
}


# ============================================================
# M-ordering initialization
# ============================================================

# Initialize feature groups and one sample ordering per group.
init_m_trajectories_cavi <- function(X,
                                     S = NULL,
                                     M = 2L,
                                     methods = NULL,
                                     pca_components = NULL,
                                     method_args = list(),
                                     init_on_failure = "pca",
                                     K = NULL,
                                     rw_q = 2L,
                                     ridge = 0,
                                     lambda_sd_prior_rate = NULL,
                                     smooth_fit_lambda_mode = c("optimize", "fixed"),
                                     smooth_fit_lambda_value = 1,
                                     lambda_min = 1e-10,
                                     lambda_max = 1e10,
                                     sigma_min = 1e-10,
                                     sigma_max = 1e10,
                                     discretization = c("quantile", "equal", "kmeans"),
                                     similarity_metric = c("spearman", "pearson", "smooth_fit", "spline_r2"),
                                     cluster_linkage = "single",
                                     similarity_min_feature_sd = 1e-8,
                                     spline_r2_df = 5L,
                                     similarity_precomputed = NULL,
                                     num_iter = 5L,
                                     verbose = FALSE) {
  X <- as.matrix(X)
  M <- as.integer(M)
  K_use <- .cavi_resolve_K(X, K, rw_q)
  discretization <- match.arg(discretization)
  similarity_metric <- match.arg(similarity_metric)
  smooth_fit_lambda_mode <- match.arg(smooth_fit_lambda_mode)

  .cavi_init_m_trajectories_similarity(
    X = X,
    S = S,
    M = M,
    methods = methods,
    method_args = method_args,
    init_on_failure = init_on_failure,
    pca_components = pca_components,
    K = K_use,
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
    num_iter = num_iter,
    similarity_metric = similarity_metric,
    cluster_linkage = cluster_linkage,
    similarity_min_feature_sd = similarity_min_feature_sd,
    spline_r2_df = spline_r2_df,
    similarity_precomputed = similarity_precomputed,
    verbose = verbose
  )
}
