# Run from the package root with argument baseline, current, or compare.
# All installations are isolated. CSV evidence is retained; complete fit
# snapshots are written under the ignored package tmp/ directory.
mode <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(mode) || !mode %in% c("baseline", "current", "compare")) {
  stop("Choose baseline, current, or compare.")
}
root <- "experiments/missing_data_v042"
dir.create(root, recursive = TRUE, showWarnings = FALSE)
if (mode == "baseline") .libPaths(c("tmp/baseline-library", .libPaths()))
if (mode == "current") pkgload::load_all(".", quiet = TRUE)
suppressPackageStartupMessages(library(MPCurver))

trace_of <- function(fit) {
  if (fit$intrinsic_dim == 1L) fit$fit$elbo_trace else fit$objective_history
}

make_data <- function(seed) {
  set.seed(seed)
  t <- seq(-1, 1, length.out = 40)
  cbind(t, t^2, sin(3 * t), cos(3 * t), -t, sin(2 * t)) +
    matrix(rnorm(240, sd = 0.25), 40, 6)
}

snapshot <- function(fit) list(
  positions = fitted_positions(fit, "probability"),
  trajectories = fitted_trajectories(fit),
  covariance = fitted_trajectories(fit, "covariance"),
  assignments = fitted_assignments(fit),
  sigma2 = fit$fit$params$sigma2,
  lambda = if (fit$intrinsic_dim == 1L) fit$fit$lambda_vec else fit$fit$lambda_mat,
  objective = trace_of(fit))

# Independent scalar expected densities plus prior/entropy terms.
scalar_objective <- function(fit) {
  raw <- fit$fit
  X <- fit$data
  n <- nrow(X)
  d <- ncol(X)
  K <- fit$K
  M <- fit$intrinsic_dim
  R <- fitted_positions(fit, "probability")
  mu <- fitted_trajectories(fit)
  cov <- fitted_trajectories(fit, "covariance")
  w <- fitted_assignments(fit)
  lambda <- if (M == 1L) matrix(raw$lambda_vec, d, 1) else raw$lambda_mat
  position_pi <- if (M == 1L) list(raw$params$pi) else raw$position_pi
  Q <- raw$Q_K
  eig <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
  positive <- eig[eig > length(eig) * max(abs(eig)) * .Machine$double.eps]
  variance <- if (is.null(raw$measurement_sd)) {
    matrix(rep(raw$params$sigma2, each = n), n, d)
  } else if (is.matrix(raw$measurement_sd)) raw$measurement_sd^2 else
    matrix(rep(raw$measurement_sd^2, each = n), n, d)
  objective <- 0
  for (m in seq_len(M)) {
    for (i in seq_len(n)) for (k in seq_len(K)) if (R[i, k, m] > 0) {
      objective <- objective + R[i, k, m] *
        (log(position_pi[[m]][k]) - log(R[i, k, m]))
    }
    for (j in seq_len(d)) {
      local <- 0
      for (i in seq_len(n)) if (!is.na(X[i, j])) for (k in seq_len(K)) {
        local <- local + R[i, k, m] *
          (dnorm(X[i, j], mu[j, k, m], sqrt(variance[i, j]), log = TRUE) -
             cov[k, k, j, m] / (2 * variance[i, j]))
      }
      V <- cov[, , j, m]
      u <- mu[j, , m]
      local <- local + 0.5 * (length(positive) * log(lambda[j, m] / (2 * pi)) +
        sum(log(positive)) - lambda[j, m] * (sum(u * (Q %*% u)) + sum(diag(Q %*% V)))) +
        0.5 * (K * log(2 * pi * exp(1)) + as.numeric(determinant(V, logarithm = TRUE)$modulus))
      objective <- objective + w[j, m] * local
    }
  }
  if (M > 1L) {
    omega <- fitted_prior(fit, type = "partition")$omega
    for (j in seq_len(d)) for (m in seq_len(M)) if (w[j, m] > 0) {
      objective <- objective + w[j, m] * (log(omega[m]) - log(w[j, m]))
    }
  }
  objective
}

grid <- expand.grid(M = 1:2, noise = c("estimated", "feature", "entry"),
                    prior = c("fixed", "adaptive"), ridge = c(0, 0.1),
                    stringsAsFactors = FALSE)
if (mode %in% c("baseline", "current")) {
  X <- make_data(4218)
  results <- lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    S <- switch(g$noise, estimated = NULL, feature = rep(0.3, ncol(X)),
                entry = matrix(seq(0.2, 0.5, length.out = length(X)), nrow(X)))
    fit <- fit_mpcurve(X, S = S, num_bins = 6, intrinsic_dim = g$M,
      position_prior = g$prior, partition_prior = g$prior, max_iter = 8, tol = 0,
      control = mpcurve_control(ridge = g$ridge))
    snapshot(fit)
  })
  saveRDS(results, paste0("tmp/v042-", mode, "-complete.rds"))
}

if (mode == "current") {
  stress <- expand.grid(seed = 4218:4219, pattern = c("random", "unequal", "block"),
                        M = 1:2, noise = c("estimated", "feature", "entry"),
                        ridge = c(0, 0.1), stringsAsFactors = FALSE)
  records <- vector("list", nrow(stress))
  for (i in seq_len(nrow(stress))) {
    g <- stress[i, ]
    X <- make_data(g$seed)
    set.seed(g$seed + 100)
    if (g$pattern == "random") X[matrix(runif(length(X)) < 0.35, nrow(X))] <- NA_real_
    if (g$pattern == "unequal") {
      X[matrix(runif(length(X)) < rep(seq(0.05, 0.75, length.out = 6), each = 40), 40)] <- NA_real_
    }
    if (g$pattern == "block") {
      X[1:15, 1:3] <- NA_real_
      X[26:40, 4:6] <- NA_real_
    }
    X[1, ] <- make_data(g$seed)[1, ]
    empty <- which(rowSums(!is.na(X)) == 0)
    if (length(empty)) X[empty, 1] <- make_data(g$seed)[empty, 1]
    S <- switch(g$noise, estimated = NULL, feature = rep(0.3, ncol(X)),
                entry = matrix(seq(0.2, 0.5, length.out = length(X)), nrow(X)))
    if (is.matrix(S)) S[is.na(X)] <- NA_real_
    warnings <- character()
    fit <- withCallingHandlers(fit_mpcurve(X, S = S, num_bins = 6,
      intrinsic_dim = g$M, max_iter = 25, tol = 0,
      position_prior = if (g$seed %% 2) "adaptive" else "fixed",
      partition_prior = if (g$seed %% 2) "adaptive" else "fixed",
      control = mpcurve_control(ridge = g$ridge)), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    objective <- trace_of(fit)
    error <- abs(tail(objective, 1) - scalar_objective(fit))
    # Tied mean-filled initialization coordinates can trigger the existing
    # documented discretization fallback. Numerical warnings remain failures.
    initialization_warnings <- grepl("Requested discretization.*Falling back", warnings)
    records[[i]] <- cbind(g, n_observed = sum(!is.na(X)),
      minimum_increment = min(diff(objective)), objective_error = error,
      initialization_warning_count = sum(initialization_warnings),
      other_warning_count = sum(!initialization_warnings))
    stopifnot(identical(fit$data, X), all(is.finite(objective)),
              error < 1e-6, min(diff(objective)) > -1e-6, all(initialization_warnings))
    cat("Verified missing-data fit", i, "of", nrow(stress), "\n")
  }
  write.csv(do.call(rbind, records), file.path(root, "missing_pattern_checks.csv"), row.names = FALSE)
  # A single observed value can drive an intrinsic fit toward an unidentified
  # state: its sample's positions occupy too few bins (singular precision) or
  # leave the system numerically unreliable. Verify that the fit stops.
  sparse <- make_data(4218)
  sparse[-1, 6] <- NA_real_
  failure <- tryCatch({
    fit_mpcurve(sparse, num_bins = 6, max_iter = 50, tol = 0)
    "No error"
  }, error = conditionMessage)
  stopifnot(grepl("singular|ill-conditioned", failure))
  proper <- fit_mpcurve(sparse, num_bins = 6, max_iter = 50, tol = 0,
                       control = mpcurve_control(ridge = 0.1))
  proper_error <- abs(tail(trace_of(proper), 1) - scalar_objective(proper))
  stopifnot(proper_error < 1e-6, min(diff(trace_of(proper))) > -1e-6)
  write.csv(data.frame(case = c("intrinsic_single_observation", "proper_single_observation"),
                       result = c(failure, "Finite monotone fit; scalar objective verified"),
                       objective_error = c(NA_real_, proper_error)),
            file.path(root, "sparse_feature_checks.csv"), row.names = FALSE)
  capture.output(sessionInfo(), file = file.path(root, "session-info.txt"))
}

if (mode == "compare") {
  old <- readRDS("tmp/v042-baseline-complete.rds")
  current <- readRDS("tmp/v042-current-complete.rds")
  error <- vapply(seq_along(old), function(i) {
    max(vapply(names(old[[i]]), function(field) {
      a <- old[[i]][[field]]
      b <- current[[i]][[field]]
      stopifnot(identical(dim(a), dim(b)), length(a) == length(b))
      if (!length(a)) 0 else max(abs(as.numeric(a) - as.numeric(b)))
    }, numeric(1)))
  }, numeric(1))
  exact <- vapply(seq_along(old), function(i) identical(old[[i]], current[[i]]), logical(1))
  stopifnot(all(error < 1e-10))
  write.csv(cbind(grid, max_absolute_difference = error, identical = exact),
            file.path(root, "complete_data_regression.csv"), row.names = FALSE)
  print(c(max_absolute_difference = max(error), identical_cases = sum(exact)))
}
