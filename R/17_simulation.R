# Internal simulators behind simulate_intrinsic_trajectories() and
# simulate_dual_trajectory(). Each ordering has its own latent positions; every
# signal feature follows one ordering and the remaining features are pure noise.

.simulate_dual_signal_block <- function(t,
                                        d,
                                        family = c("sinusoidal", "linear", "monotone", "quadratic"),
                                        sigma = 0.3,
                                        signal_range = c(0.5, 2),
                                        linear_slope_range = c(0.8, 2),
                                        monotone_power_range = c(0.7, 2.5),
                                        quadratic_curvature_range = c(0.8, 2),
                                        quadratic_center_range = c(0.25, 0.75),
                                        intercept_sd = 0,
                                        sinusoid_freq = 1:4) {
  family <- match.arg(family)
  n <- length(t)
  d <- as.integer(d)
  if (d <= 0L) return(matrix(0, nrow = n, ncol = 0L))

  amps <- stats::runif(d, min(signal_range), max(signal_range))
  intercepts <- stats::rnorm(d, mean = 0, sd = intercept_sd)
  out <- matrix(0, nrow = n, ncol = d)

  for (j in seq_len(d)) {
    base <- switch(
      family,
      sinusoidal = {
        freq_j <- sinusoid_freq[sample.int(length(sinusoid_freq), 1L)]
        phase_j <- stats::runif(1L, 0, 2 * pi)
        sin(freq_j * pi * t + phase_j)
      },
      linear = {
        slope_j <- stats::runif(1L, min(linear_slope_range), max(linear_slope_range))
        slope_j <- slope_j * sample(c(-1, 1), 1L)
        slope_j * (t - 0.5)
      },
      monotone = {
        power_left <- stats::runif(1L, min(monotone_power_range), max(monotone_power_range))
        power_right <- stats::runif(1L, min(monotone_power_range), max(monotone_power_range))
        logistic_mid <- stats::runif(1L, 0.2, 0.8)
        logistic_steep <- stats::runif(1L, 3, 12)
        exp_rate <- stats::runif(1L, 1, 6)
        hill_left <- stats::runif(1L, min(monotone_power_range), max(monotone_power_range))
        hill_right <- stats::runif(1L, min(monotone_power_range), max(monotone_power_range))

        logistic_basis <- stats::plogis((t - logistic_mid) * logistic_steep)
        logistic_basis <- (logistic_basis - min(logistic_basis)) /
          (max(logistic_basis) - min(logistic_basis) + 1e-12)

        exp_basis <- (exp(exp_rate * t) - 1) / (exp(exp_rate) - 1)

        hill_num <- t^hill_left
        hill_den <- hill_num + (1 - t)^hill_right + 1e-12
        hill_basis <- hill_num / hill_den

        monotone_basis <- cbind(
          t^power_left,
          1 - (1 - t)^power_right,
          logistic_basis,
          exp_basis,
          hill_basis
        )
        mix_w <- stats::rexp(ncol(monotone_basis))
        mix_w <- mix_w / sum(mix_w)
        raw <- as.numeric(monotone_basis %*% mix_w)
        if (sample(c(FALSE, TRUE), 1L)) raw <- 1 - raw
        raw - mean(raw)
      },
      quadratic = {
        center_j <- stats::runif(1L, min(quadratic_center_range), max(quadratic_center_range))
        curv_j <- stats::runif(1L, min(quadratic_curvature_range), max(quadratic_curvature_range))
        curv_j <- curv_j * sample(c(-1, 1), 1L)
        raw <- curv_j * (t - center_j)^2
        raw - mean(raw)
      }
    )

    out[, j] <- intercepts[j] + amps[j] * base + stats::rnorm(n, 0, sigma)
  }

  out
}

# Internal numerical implementation of simulate_intrinsic_trajectories.
.simulate_intrinsic_trajectories <- function(n = 200,
                                            d_signal = c(30, 30),
                                            d_noise = 20,
                                            sigma = 0.3,
                                            seed = 42L,
                                            trajectory_family = "sinusoidal",
                                            signal_range = c(0.5, 2),
                                            linear_slope_range = c(0.8, 2),
                                            monotone_power_range = c(0.7, 2.5),
                                            quadratic_curvature_range = c(0.8, 2),
                                            quadratic_center_range = c(0.25, 0.75),
                                            intercept_sd = 0,
                                            sinusoid_freq = 1:4,
                                            latent_positions = NULL) {
  n <- as.integer(n)
  d_signal <- as.integer(d_signal)
  d_noise <- as.integer(d_noise)
  if (length(n) != 1L || is.na(n) || n < 2L) stop("n must be a single integer >= 2.")
  if (!length(d_signal) || any(is.na(d_signal)) || any(d_signal < 0L)) {
    stop("d_signal must be a nonnegative integer vector of positive length.")
  }
  if (length(d_noise) != 1L || is.na(d_noise) || d_noise < 0L) {
    stop("d_noise must be a single nonnegative integer.")
  }
  if ((sum(d_signal) + d_noise) < 1L) {
    stop("At least one of the signal or noise feature counts must be positive.")
  }
  if (!is.finite(sigma) || sigma < 0) stop("sigma must be a nonnegative number.")
  if (!is.null(seed)) set.seed(seed)

  M <- length(d_signal)
  trajectory_family <- rep(trajectory_family, length.out = M)
  valid_family <- c("sinusoidal", "linear", "monotone", "quadratic")
  if (!all(trajectory_family %in% valid_family)) {
    stop("trajectory_family must contain only: ",
         paste(valid_family, collapse = ", "), ".")
  }
  if (!is.numeric(signal_range) || length(signal_range) != 2L ||
      any(!is.finite(signal_range)) || any(signal_range <= 0)) {
    stop("signal_range must be a positive numeric vector of length 2.")
  }
  if (!is.numeric(linear_slope_range) || length(linear_slope_range) != 2L ||
      any(!is.finite(linear_slope_range)) || any(linear_slope_range <= 0)) {
    stop("linear_slope_range must be a positive numeric vector of length 2.")
  }
  if (!is.numeric(monotone_power_range) || length(monotone_power_range) != 2L ||
      any(!is.finite(monotone_power_range)) || any(monotone_power_range <= 0)) {
    stop("monotone_power_range must be a positive numeric vector of length 2.")
  }
  if (!is.numeric(quadratic_curvature_range) || length(quadratic_curvature_range) != 2L ||
      any(!is.finite(quadratic_curvature_range)) || any(quadratic_curvature_range <= 0)) {
    stop("quadratic_curvature_range must be a positive numeric vector of length 2.")
  }
  if (!is.numeric(quadratic_center_range) || length(quadratic_center_range) != 2L ||
      any(!is.finite(quadratic_center_range)) ||
      quadratic_center_range[1] < 0 || quadratic_center_range[2] > 1 ||
      quadratic_center_range[1] > quadratic_center_range[2]) {
    stop("quadratic_center_range must be a finite numeric vector of length 2 inside [0, 1].")
  }
  if (!is.finite(intercept_sd) || intercept_sd < 0) {
    stop("intercept_sd must be nonnegative.")
  }
  if (!length(sinusoid_freq) || any(!is.finite(sinusoid_freq))) {
    stop("sinusoid_freq must contain at least one finite value.")
  }

  if (is.null(latent_positions)) {
    latent_positions <- matrix(stats::runif(n * M), nrow = n, ncol = M)
  } else {
    latent_positions <- as.matrix(latent_positions)
    if (!all(dim(latent_positions) == c(n, M))) {
      stop("latent_positions must be an n x length(d_signal) matrix.")
    }
    if (any(!is.finite(latent_positions)) ||
        any(latent_positions < 0) || any(latent_positions > 1)) {
      stop("latent_positions must contain finite values inside [0, 1].")
    }
  }

  ordering_labels <- if (M <= 26L) LETTERS[seq_len(M)] else paste0("ord", seq_len(M))
  signal_blocks <- lapply(seq_len(M), function(m) {
    .simulate_dual_signal_block(
      t = latent_positions[, m],
      d = d_signal[m],
      family = trajectory_family[m],
      sigma = sigma,
      signal_range = signal_range,
      linear_slope_range = linear_slope_range,
      monotone_power_range = monotone_power_range,
      quadratic_curvature_range = quadratic_curvature_range,
      quadratic_center_range = quadratic_center_range,
      intercept_sd = intercept_sd,
      sinusoid_freq = sinusoid_freq
    )
  })
  X_noise <- if (d_noise > 0L) {
    matrix(stats::rnorm(n * d_noise), nrow = n, ncol = d_noise)
  } else {
    matrix(0, nrow = n, ncol = 0L)
  }

  X_full <- do.call(cbind, c(signal_blocks, list(X_noise)))
  assign_full <- unlist(Map(function(lbl, d_m) rep(lbl, d_m), ordering_labels, d_signal))
  assign_full <- c(assign_full, rep("noise", d_noise))
  perm <- sample.int(sum(d_signal) + d_noise)
  colnames(X_full) <- paste0("V", seq_len(ncol(X_full)))

  list(
    X = X_full[, perm, drop = FALSE],
    true_assign = assign_full[perm],
    latent_positions = latent_positions,
    ordering_labels = ordering_labels,
    original_order = perm,
    trajectory_family = trajectory_family
  )
}

# Internal numerical implementation of simulate_dual_trajectory.
.simulate_dual_trajectory <- function(n = 200,
                                     d1 = 30,
                                     d2 = 30,
                                     d_noise = 20,
                                     sigma = 0.3,
                                     crossing = FALSE,
                                     seed = 42L,
                                     trajectory_family = c("sinusoidal", "sinusoidal"),
                                     signal_range = c(0.5, 2),
                                     linear_slope_range = c(0.8, 2),
                                     monotone_power_range = c(0.7, 2.5),
                                     quadratic_curvature_range = c(0.8, 2),
                                     quadratic_center_range = c(0.25, 0.75),
                                     intercept_sd = 0,
                                     sinusoid_freq = 1:4) {
  n <- as.integer(n)
  d1 <- as.integer(d1)
  d2 <- as.integer(d2)
  d_noise <- as.integer(d_noise)
  if (length(n) != 1L || is.na(n) || n < 2L) stop("n must be a single integer >= 2.")
  if (any(c(d1, d2, d_noise) < 0L) || any(is.na(c(d1, d2, d_noise)))) {
    stop("d1, d2, and d_noise must be nonnegative integers.")
  }
  if ((d1 + d2 + d_noise) < 1L) {
    stop("At least one of d1, d2, or d_noise must be positive.")
  }
  if (!is.null(seed)) set.seed(seed)

  t1 <- stats::runif(n)
  if (crossing) {
    t2_raw <- 1 - t1 + stats::rnorm(n, 0, 0.3)
    t2 <- (t2_raw - min(t2_raw)) / (max(t2_raw) - min(t2_raw))
  } else {
    t2 <- stats::runif(n)
  }

  sim <- .simulate_intrinsic_trajectories(
    n = n,
    d_signal = c(d1, d2),
    d_noise = d_noise,
    sigma = sigma,
    seed = NULL,
    trajectory_family = trajectory_family,
    signal_range = signal_range,
    linear_slope_range = linear_slope_range,
    monotone_power_range = monotone_power_range,
    quadratic_curvature_range = quadratic_curvature_range,
    quadratic_center_range = quadratic_center_range,
    intercept_sd = intercept_sd,
    sinusoid_freq = sinusoid_freq,
    latent_positions = cbind(t1, t2)
  )

  sim$t1 <- t1
  sim$t2 <- t2
  sim
}
