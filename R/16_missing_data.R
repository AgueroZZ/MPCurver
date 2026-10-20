# Observed-entry calculations shared by single and structural CAVI. Missing
# entries are integrated out; they never contribute residuals or precision.
.mpcurve_initialization_data <- function(X) {
  X <- as.matrix(X)
  if (!anyNA(X)) return(X)
  means <- colMeans(X, na.rm = TRUE)
  if (any(!is.finite(means))) {
    stop("Every feature must have at least one observed value.", call. = FALSE)
  }
  for (j in seq_len(ncol(X))) X[is.na(X[, j]), j] <- means[j]
  X
}

.mpcurve_observed_variance <- function(X) {
  vapply(seq_len(ncol(X)), function(j) {
    x <- X[!is.na(X[, j]), j]
    # A single observation has no empirical variance; one is only a starting
    # value. Subsequent updates use its observed expected squared residual.
    if (length(x) > 1L) stats::var(x) else 1
  }, numeric(1))
}

# Exported ordering helpers need complete data. fit_mpcurve() fills missing
# entries before calling them; direct calls are warned and filled the same way.
.mpcurve_ordering_input <- function(X, caller) {
  X <- as.matrix(X)
  if (anyNA(X)) {
    warning(sprintf(paste0("%s() does not accept missing values; %d missing ",
                           "value(s) were filled with feature means of the observed ",
                           "values for this ordering."), caller, sum(is.na(X))),
            call. = FALSE)
    X <- .mpcurve_initialization_data(X)
  }
  X
}

.mpcurve_missing_initial_responsibilities <- function(R) {
  # Soften generated starts, not the likelihood precision or user-supplied
  # responsibilities. Sparse columns can otherwise occupy one RW2 bin only.
  0.99 * R + 0.01 / ncol(R)
}

# Sums over observed entries are written as masked matrix products: a missing
# cell has zero precision and zero filled value, so it adds nothing to any sum
# over samples or features. Squared residuals are expanded after centering each
# feature at its observed mean. The expansion is exact algebra, and centering
# keeps it from cancelling digits when a feature has a large common offset;
# residuals are unchanged by shifting data and trajectory means together.
.cavi_observed_terms <- function(X, sigma2 = NULL, var_mat = NULL,
                                 weights = rep(1, ncol(X))) {
  observed <- !is.na(X)
  center <- colMeans(X, na.rm = TRUE)
  centered <- sweep(X, 2L, center, "-")
  centered[!observed] <- 0
  variance <- if (is.null(var_mat)) {
    matrix(sigma2, nrow(X), ncol(X), byrow = TRUE)
  } else var_mat
  precision <- sweep(1 / variance, 2L, weights, "*")
  precision[!observed] <- 0
  log_normalizer <- sweep(log(2 * base::pi * variance), 2L, weights, "*")
  log_normalizer[!observed] <- 0
  list(observed = observed * 1, centered = centered, center = center,
       precision = precision, log_normalizer = log_normalizer)
}

# Expected squared residuals E[(x_ij - U_jk)^2] summed over observed entries,
# per feature and position: A2 - 2 m A1 + (m^2 + s) A0 in centered units.
.cavi_observed_expected_squares <- function(centered, mass, q_u, center, R) {
  mean_centered <- q_u$m_mat - center
  crossprod(centered^2 * mass, R) - 2 * mean_centered * crossprod(centered * mass, R) +
    (mean_centered^2 + q_u$sdiag_mat) * crossprod(mass, R)
}

# Sample-by-position expected log density, summed over each sample's observed
# features with optional feature weights.
.cavi_observed_logdensity <- function(X, q_u, sigma2 = NULL, var_mat = NULL,
                                     weights = rep(1, ncol(X)),
                                     uncertainty = TRUE) {
  terms <- .cavi_observed_terms(X, sigma2, var_mat, weights)
  mean_centered <- q_u$m_mat - terms$center
  second_moment <- mean_centered^2 + if (uncertainty) q_u$sdiag_mat else 0
  weighted <- terms$precision * terms$centered
  quadratic <- rowSums(weighted * terms$centered) - 2 * (weighted %*% mean_centered) +
    terms$precision %*% second_moment
  unname(-0.5 * (rowSums(terms$log_normalizer) + quadratic))
}

# Feature-level expected observed log likelihood under responsibilities R.
.cavi_observed_feature_scores <- function(X, R, q_u, sigma2 = NULL,
                                          var_mat = NULL) {
  terms <- .cavi_observed_terms(X, sigma2, var_mat)
  squares <- .cavi_observed_expected_squares(
    terms$centered, terms$precision, q_u, terms$center, R)
  unname(-0.5 * (colSums(terms$log_normalizer * rowSums(R)) + rowSums(squares)))
}

.cavi_observed_q_u <- function(X, R, sigma2, lambda_vec, Q_K,
                               var_mat = NULL, ordering_label = NULL) {
  d <- ncol(X)
  K <- ncol(R)
  observed <- !is.na(X)
  n_observed <- colSums(observed)
  precision <- if (is.null(var_mat)) {
    sweep(observed * 1, 2L, 1 / sigma2, "*")
  } else {
    observed / var_mat
  }
  precision[!observed] <- 0
  filled <- X
  filled[!observed] <- 0
  # Data precision D[j, ] and right-hand side b[j, ] for every feature at once.
  diagonal_mat <- crossprod(precision, R)
  rhs_mat <- crossprod(precision * filled, R)
  occupied <- .cavi_occupied_bins(crossprod(observed * 1, R), n_observed)
  out <- list(m_mat = matrix(0, d, K), sdiag_mat = matrix(0, d, K),
              S_list = vector("list", d), logdetS = numeric(d),
              eq_quad = numeric(d))
  for (j in seq_len(d)) {
    diagonal <- diagonal_mat[j, ]
    fail_precision <- function() {
      stop(sprintf(paste0("Observed-data q(U) precision is not positive definite ",
                          "for feature %d (%d observed values). The likelihood ",
                          "must identify the intrinsic prior null space."),
                   j, n_observed[j]), call. = FALSE)
    }
    # Cholesky alone can accept a mathematically singular intrinsic system
    # after roundoff. Test null(Q) against the bins carrying data first.
    .cavi_check_precision_support(occupied[j, ], Q_K, j, ordering_label, n_observed[j])
    P <- lambda_vec[j] * Q_K + diag(diagonal, K)
    P <- 0.5 * (P + t(P))
    ch <- tryCatch(chol(P), error = function(e) fail_precision())
    V <- chol2inv(ch)
    # A successful factorization does not guarantee a usable inverse. Reject
    # systems with condition number at least 1e12, where fewer than roughly
    # four significant digits may remain; do not alter the model with a ridge.
    reciprocal_condition <- 1 / (norm(P, "1") * norm(V, "1"))
    if (!is.finite(reciprocal_condition) ||
        reciprocal_condition <= 1e-12) {
      stop(sprintf(paste0("Observed-data q(U) precision is too ill-conditioned ",
                          "for reliable inference for feature %d (%d observed ",
                          "values; reciprocal condition %.3g)."),
                   j, n_observed[j], reciprocal_condition), call. = FALSE)
    }
    mu <- as.numeric(backsolve(ch, forwardsolve(t(ch), rhs_mat[j, ])))
    out$m_mat[j, ] <- mu
    out$sdiag_mat[j, ] <- diag(V)
    out$S_list[[j]] <- V
    out$logdetS[j] <- -2 * sum(log(diag(ch)))
    out$eq_quad[j] <- as.numeric(crossprod(mu, Q_K %*% mu) + sum(Q_K * V))
  }
  out
}

.cavi_observed_residuals <- function(X, R, q_u) {
  observed <- !is.na(X)
  center <- colMeans(X, na.rm = TRUE)
  centered <- sweep(X, 2L, center, "-")
  centered[!observed] <- 0
  unname(rowSums(.cavi_observed_expected_squares(centered, observed * 1, q_u, center, R)))
}

.cavi_observed_elbo <- function(X, R, pi_vec, q_u, sigma2, measurement_sd,
                                lambda_vec, Q_K, rw_q, lambda_sd_prior_rate) {
  scores <- .cavi_feature_scores_from_state(
    X, R, q_u, sigma2, lambda_vec, measurement_sd, Q_K, rw_q,
    lambda_sd_prior_rate)
  .cavi_cell_terms_from_state(R, pi_vec) + sum(scores) +
    .cavi_entropy_u(q_u, ncol(R))
}

.mpcurve_observation_summary <- function(X) {
  observed <- !is.na(X)
  list(n_observed = sum(observed), n_missing = sum(!observed),
       per_sample = rowSums(observed), per_feature = colSums(observed),
       initialization = if (all(observed)) "complete_data" else "feature_mean_fill",
       likelihood = "observed_entries")
}
