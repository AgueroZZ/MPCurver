# Select an intrinsic dimension by comparing fixed-M structural fits.

.mpcurve_dimension_single_score_check <- function(fit, score) {
  raw <- fit$fit
  q_u <- .cavi_posterior_state_from_fit(raw)
  feature_scores <- .cavi_feature_scores_from_state(
    X = raw$data,
    R = raw$gamma,
    q_u = q_u,
    sigma2 = raw$params$sigma2,
    lambda_vec = raw$lambda_vec,
    measurement_sd = raw$measurement_sd,
    Q_K = raw$Q_K,
    rw_q = raw$rw_q,
    lambda_sd_prior_rate = raw$control$lambda_sd_prior_rate,
    include_prior = TRUE
  )
  entropy_u <- .cavi_entropy_u_terms(q_u, ncol(raw$gamma))
  cell_terms <- .cavi_cell_terms_from_state(raw$gamma, raw$params$pi)
  reconstructed <- cell_terms + sum(feature_scores + entropy_u)
  difference <- reconstructed - score
  if (!is.finite(difference) ||
      abs(difference) > 1e-6 * (1 + abs(score))) {
    stop(
      sprintf(
        "The M=1 objective does not match the structural score convention (difference %.6g).",
        difference
      ),
      call. = FALSE
    )
  }
  invisible(difference)
}

.mpcurve_dimension_empty_candidates <- function() {
  data.frame(
    M = integer(),
    initialization = character(),
    status = character(),
    score = numeric(),
    converged = logical(),
    K = integer(),
    final_temperature = numeric(),
    selected_for_M = logical(),
    warnings = character(),
    error = character(),
    stringsAsFactors = FALSE
  )
}

.mpcurve_dimension_empty_history <- function() {
  data.frame(
    step = integer(),
    direction = character(),
    current_M = integer(),
    candidate_M = integer(),
    current_score = numeric(),
    candidate_score = numeric(),
    delta = numeric(),
    accepted = logical(),
    stringsAsFactors = FALSE
  )
}

.mpcurve_dimension_fit_candidate <- function(X, M, initialization, fit_options) {
  warnings <- character()
  error <- NA_character_
  fit <- tryCatch(
    withCallingHandlers(
      {
        do.call(fit_mpcurve, c(list(X = X, intrinsic_dim = M,
                                   partition_prior = "fixed"), fit_options))
      },
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
      }
    ),
    error = function(e) {
      error <<- conditionMessage(e)
      NULL
    }
  )

  score <- NA_real_
  actual_K <- NA_integer_
  converged <- FALSE
  final_temperature <- NA_real_
  if (!is.null(fit)) {
    raw <- fit$fit
    is_partition <- M >= 2L
    trace <- if (is_partition) raw$objective_history else raw$elbo_trace
    final_temperature <- if (is_partition) {
      as.numeric(tail(raw$temperature_history, 1L))
    } else {
      1
    }
    actual_K <- as.integer(fit$K)
    converged <- isTRUE(raw$converged)
    score <- if (length(trace)) as.numeric(tail(trace, 1L)) else NA_real_
    validation_error <- tryCatch({
      if (length(score) != 1L || !is.finite(score)) {
        stop("The final objective is missing or non-finite.")
      }
      if (length(actual_K) != 1L || is.na(actual_K) || actual_K < 2L) {
        stop("The fitted K is missing or invalid.")
      }
      if (length(final_temperature) != 1L ||
          !isTRUE(all.equal(final_temperature, 1, tolerance = 1e-12))) {
        stop("The final structural objective is not evaluated at T = 1.")
      }
      if (is_partition) {
        if (!identical(raw$variational_family, "structured")) {
          stop("The candidate is not a structural partition fit.")
        }
        if (!identical(raw$control$partition_prior, "fixed") ||
            !is.null(raw$control$partition_prior_init)) {
          stop("The candidate does not use a fixed uniform assignment prior.")
        }
        if (!isTRUE(all.equal(rowSums(raw$pi_weights),
                              rep(1, nrow(raw$pi_weights)),
                              tolerance = 1e-8))) {
          stop("The candidate has invalid soft assignment weights.")
        }
      } else {
        .mpcurve_dimension_single_score_check(fit, score)
      }
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(validation_error)) {
      error <- validation_error
      fit <- NULL
    }
  }

  row <- data.frame(
    M = as.integer(M),
    initialization = initialization,
    status = if (is.null(fit)) "failed" else "success",
    score = if (is.null(fit)) NA_real_ else score,
    converged = converged,
    K = actual_K,
    final_temperature = final_temperature,
    selected_for_M = FALSE,
    warnings = paste(unique(warnings), collapse = " | "),
    error = error,
    stringsAsFactors = FALSE
  )
  list(fit = fit, row = row)
}

#' Select the number of intrinsic orderings under a uniform prior
#'
#' Fit adjacent candidate dimensions with a fixed uniform prior over feature
#' orderings. Forward selection starts at one ordering; backward selection
#' starts at `max_intrinsic_dim`. At each step, the candidate is accepted only
#' if its final soft-assignment objective at `T = 1` is strictly larger.
#' Each candidate dimension is fitted once, grouping features by similarity
#' and applying the selected `initial_method` independently within each group.
#'
#' This is a greedy comparison of variational solutions. The selected dimension
#' can depend on initialization and convergence; inspect the returned candidate
#' table and comparison history before interpreting it.
#'
#' @param X Sample-by-feature data matrix, as for [fit_mpcurve()].
#' @param max_intrinsic_dim Positive integer upper bound for the search.
#' @param direction Either `"forward"` or `"backward"`.
#' @inheritParams fit_mpcurve
#' @param init_control Initialization settings from [mpcurve_init_control()].
#'   Supplied fits, responsibilities, and PCA components cannot be shared across
#'   candidate counts.
#' @param control Model and optimization settings from [mpcurve_control()].
#'   Comparisons require a fixed uniform partition prior and soft final
#'   assignments. Custom partition weights and hard assignments are rejected.
#'
#' @return The selected [fit_mpcurve()] `mpcurve` object. Its
#'   `$dimension_selection` record contains `direction`,
#'   `max_intrinsic_dim`, `selected_M`, `stop_reason`, a
#'   `continued_after_selection` flag, a `$candidates` table
#'   with every attempted fit and its score, convergence state, warnings, and
#'   error, and a `$history` table including the rejected stopping step.
#' @md
#' @export
select_mpcurve_dimension <- function(
    X, max_intrinsic_dim, direction = c("forward", "backward"),
    S = NULL, num_bins = NULL, initial_method = "PCA",
    position_prior = c("adaptive", "fixed"), max_iter = 100L, tol = 1e-6,
    verbose = FALSE, init_control = NULL, control = NULL) {
  X <- .mpcurve_validate_data(X, S)
  direction <- match.arg(direction)
  max_intrinsic_dim <- .mpcurve_integer_setting(max_intrinsic_dim, "max_intrinsic_dim")
  position_prior <- match.arg(position_prior)
  fit_options <- list(S = S, num_bins = num_bins, initial_method = initial_method,
    position_prior = position_prior, max_iter = max_iter, tol = tol,
    verbose = verbose, init_control = init_control, control = control)
  init <- .mpcurve_interface_options(fit_options$init_control,
                                     mpcurve_init_control, "init_control")
  if (!is.null(init$fits_init) || !is.null(init$responsibilities_init) ||
      !is.null(init$pca_components)) {
    stop("The selector cannot share supplied fits, responsibilities, or PCA components across dimensions.",
         call. = FALSE)
  }
  ctrl <- .mpcurve_interface_options(fit_options$control, mpcurve_control, "control")
  if (!is.null(ctrl$partition_prior_weights)) {
    stop("Dimension selection requires a uniform partition prior; partition_prior_weights must be NULL.",
         call. = FALSE)
  }
  if (!is.null(fit_options$initial_method) && length(fit_options$initial_method) != 1L) {
    stop("initial_method must select exactly one initialization method.", call. = FALSE)
  }
  fit_options$partition_prior <- NULL
  fit_options$control <- ctrl
  fit_options$init_control <- init

  candidates <- .mpcurve_dimension_empty_candidates()
  history <- .mpcurve_dimension_empty_history()
  fitted_K <- NULL

  fit_dimension <- function(M) {
    initialization <- if (M == 1L) "single" else "similarity"
    result <- .mpcurve_dimension_fit_candidate(X, M, initialization, fit_options)
    if (identical(result$row$status, "failed")) {
      stop(sprintf("Fit failed for M=%d using %s initialization: %s",
                   M, initialization, result$row$error), call. = FALSE)
    }
    result$row$selected_for_M <- TRUE
    list(
      fit = result$fit,
      score = result$row$score,
      K = result$row$K,
      rows = result$row
    )
  }

  fit_and_record <- function(M) {
    result <- fit_dimension(M)
    if (is.null(fitted_K)) {
      fitted_K <<- result$K
    } else if (result$K != fitted_K) {
      stop(
        sprintf("Cannot compare M=%d and other candidates: fitted K differs (%d versus %d).",
                M, result$K, fitted_K),
        call. = FALSE
      )
    }
    candidates <<- rbind(candidates, result$rows)
    result
  }

  current_M <- if (identical(direction, "forward")) 1L else max_intrinsic_dim
  current <- fit_and_record(current_M)
  stop_reason <- if (identical(direction, "forward"))
    "reached_upper_bound" else "reached_dimension_1"

  repeat {
    candidate_M <- if (identical(direction, "forward")) current_M + 1L else current_M - 1L
    if (candidate_M < 1L || candidate_M > max_intrinsic_dim) break
    candidate <- fit_and_record(candidate_M)
    delta <- candidate$score - current$score
    accepted <- is.finite(delta) && delta > 0
    history <- rbind(history, data.frame(
      step = nrow(history) + 1L,
      direction = direction,
      current_M = current_M,
      candidate_M = candidate_M,
      current_score = current$score,
      candidate_score = candidate$score,
      delta = delta,
      accepted = accepted,
      stringsAsFactors = FALSE
    ))
    if (!accepted) {
      stop_reason <- "no_improvement"
      break
    }
    current_M <- candidate_M
    current <- candidate
  }

  selected <- current$fit
  selected$dimension_selection <- list(
    direction = direction,
    max_intrinsic_dim = max_intrinsic_dim,
    selected_M = current_M,
    criterion = "final_soft_T1_ELBO",
    stop_reason = stop_reason,
    continued_after_selection = FALSE,
    candidates = candidates,
    history = history
  )
  selected
}
