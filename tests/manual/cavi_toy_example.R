pkgload::load_all(".", quiet = TRUE)

set.seed(1)

sim <- simulate_mpcurve(
  n = 150,
  d = 20,
  num_bins = 8,
  seed = 1,
  control = list(rw_order = 2, lambda_range = c(0.8, 3), noise_sd_range = c(0.08, 0.18))
)

fit <- cavi(
  sim$X,
  K = 8,
  method = "PCA",
  rw_q = 2,
  max_iter = 50,
  tol = 1e-6,
  verbose = TRUE
)

z_hat <- max.col(fit$gamma, ties.method = "first")

cat("\n=== MPCurve simulation and CAVI fitting example ===\n")
cat(sprintf("iter: %d\n", fit$iter))
cat(sprintf("converged: %s\n", fit$converged))
cat(sprintf("ELBO start/end: %.6f -> %.6f\n", fit$elbo_trace[1], tail(fit$elbo_trace, 1)))
cat(sprintf("ELBO min delta: %.3e\n", min(diff(fit$elbo_trace))))
cat(sprintf("direct MAP assignment accuracy: %.3f\n", mean(z_hat == sim$z)))
cat(sprintf("sigma2 range: [%.3g, %.3g]\n", min(fit$params$sigma2), max(fit$params$sigma2)))
cat(sprintf("lambda range: [%.3g, %.3g]\n", min(fit$lambda_vec), max(fit$lambda_vec)))
