# Run from the MPCurver repository root. Compare the current automatic default
# with the frozen k=2 result from the adjacent internal Isomap experiment.
# Known B-group membership defines the input; truth only evaluates recovery.
pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
input_path <- "../InferOrder/experiments/estimate_intrinsic_m_smooth_v032/data/main_M5_S4_r001.rds"
reference_path <- "../InferOrder/experiments/isomap_elbo_screen_v040/narrow_range.rds"
dataset <- readRDS(input_path)
reference <- readRDS(reference_path)
features <- which(dataset$true_assign == "B")
X <- dataset$X[, features, drop = FALSE]
truth <- dataset$latent_positions[, match("B", dataset$ordering_labels)]
stopifnot(identical(digest::digest(X, algo = "sha256"),
                    reference$provenance$group_input_sha256))
seed <- as.integer(dataset$row$seed + 1000000L)
init <- mpcurve_init_control(discretization = "quantile", on_failure = "error",
  method_args = list(seed = seed, control = list(num_landmarks = nrow(X), keep = "all")))
control <- mpcurve_control(rw_order = 2L, lambda_init = 1, ridge = 0,
                           convergence = "normalized")
set.seed(seed)
early <- fit_mpcurve(X, num_bins = 50, initial_method = "isomap",
  position_prior = "adaptive", max_iter = 1, tol = 0,
  init_control = init, control = control)
final <- do_mpcurve(early, max_iter = 2000, tol = 1e-6)
expected <- reference$candidates[[match(2L, reference$candidate_neighbors)]]
checks <- c(
  early_elbo = max(abs(early$elbo_trace - expected$early$elbo_trace)),
  final_elbo_trace = max(abs(final$elbo_trace - expected$final$elbo_trace)),
  final_positions = max(abs(as.numeric(fitted_positions(final)) - expected$final$position)),
  noise_variances = max(abs(final$params$sigma2 - expected$final$sigma2)),
  precisions = max(abs(final$fit$lambda_vec - expected$final$lambda)),
  position_prior = max(abs(final$params$pi - expected$final$position_prior)))
stopifnot(early$fit$init_info$k_used == 2L, final$converged,
  final$fit$iter == expected$final$iter, all(checks < 1e-8),
  identical(early$fit$init_info, final$fit$init_info))
print(data.frame(k_used = final$fit$init_info$k_used,
  one_sweep_elbo = tail(early$elbo_trace, 1), final_elbo = tail(final$elbo_trace, 1),
  final_recovery = abs(cor(truth, as.numeric(fitted_positions(final)), method = "spearman")),
  sweeps = final$fit$iter, converged = final$converged), digits = 12)
print(checks)
cat("Input SHA256:", digest::digest(X, algo = "sha256"), "\n")
cat("Package version:", as.character(packageVersion("MPCurver")), "\n")
