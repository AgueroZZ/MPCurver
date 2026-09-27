# Run from the MPCurver repository root with the installed release candidate.
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1,
           VECLIB_MAXIMUM_THREADS = 1, MKL_NUM_THREADS = 1)
library(MPCurver)
stopifnot(as.character(packageVersion('MPCurver')) == '0.3.2')
root <- 'experiments/v032_convergence_regression'
inputs <- readRDS(file.path(root, 'inputs.rds'))
dir.create(file.path(root, 'fits'), showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(root, 'session_info.txt'))
rows <- list()
record <- function(row) {
  rows[[length(rows) + 1L]] <<- row
  write.csv(do.call(rbind, rows), file.path(root, 'results.csv'), row.names = FALSE)
  print(row); flush.console()
}
fit_once <- function(key, seed, args) {
  path <- file.path(root, 'fits', paste0(key, '.rds'))
  if (file.exists(path)) return(readRDS(path))
  set.seed(seed)
  warnings <- character()
  cat(format(Sys.time()), 'START', key, '\n'); flush.console()
  timing <- system.time(fit <- withCallingHandlers(do.call(fit_mpcurve, args),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning') }))
  stopifnot(fit$fit$control$convergence == 'normalized')
  save <- list(fit = fit, seconds = unname(timing['elapsed']), warnings = warnings)
  saveRDS(save, path)
  save
}
validate <- function(fit) {
  gamma <- if (is.matrix(fit$gamma)) list(fit$gamma) else fit$gamma
  stopifnot(all(vapply(gamma, function(g) all(is.finite(g)) &&
    min(g) >= 0 && max(abs(rowSums(g) - 1)) < 1e-8, logical(1))))
  trace <- fit$elbo_trace
  if (!is.null(fit$temperature_history)) trace <- trace[fit$temperature_history == 1]
  stopifnot(all(is.finite(trace)), min(diff(trace)) > -1e-5,
            all(is.finite(fit$params$sigma2)), all(fit$params$sigma2 > 0))
  if (fit$intrinsic_dim > 1) {
    W <- fit$partition$pi_weights
    stopifnot(all(is.finite(W)), min(W) >= 0, max(abs(rowSums(W) - 1)) < 1e-8)
  }
}
base_row <- function(id, kind, result, reference_M = NA_integer_, reference_recovery = NA_real_,
                     recovered = NA_real_, ARI = NA_real_, reference_ARI = NA_real_,
                     prefix_error = NA_real_, passed = TRUE) {
  fit <- result$fit
  data.frame(id = id, kind = kind, fitted_M = fit$intrinsic_dim,
    effective_M = fit$effective_intrinsic_dim, reference_M = reference_M,
    iterations = fit$fit$iter, converged = fit$fit$converged, seconds = result$seconds,
    recovery = recovered, reference_recovery = reference_recovery,
    ARI = ARI, reference_ARI = reference_ARI, prefix_error = prefix_error,
    warnings = length(result$warnings), passed = passed)
}
# The original two website examples at their original iteration budgets.
for (M in 1:2) {
  old <- inputs$site[[M]]
  X <- if (M == 1) as.matrix(old$simulation$obs) else old$simulation$X
  args <- list(X = X, intrinsic_dim = M, method = 'fiedler', K = if (M == 1) 60 else 50,
               iter = if (M == 1) 120 else 100)
  result <- fit_once(paste0('site_M', M), 42, args)
  fit <- result$fit; validate(fit)
  if (M == 1) {
    rho <- abs(cor(as.numeric(fit$gamma %*% seq_len(fit$K)), old$simulation$t, method = 'spearman'))
    previous <- old$summary$abs_spearman
    ari <- previous_ari <- NA_real_
  } else {
    ari <- mclust::adjustedRandIndex(fit$partition$assign, old$simulation$true_assign)
    previous_ari <- 1
    cor_mat <- abs(cor(old$simulation$latent_positions,
      vapply(fit$locations, function(z) z$mean$pseudotime, numeric(nrow(X))), method='spearman'))
    assignment <- as.integer(clue::solve_LSAP(cor_mat, maximum=TRUE))
    rho <- mean(cor_mat[cbind(1:2, assignment)])
    previous <- mean(unlist(old$summary[c('abs_spearman_A','abs_spearman_B')]))
  }
  passed <- rho >= previous - 0.005 && (is.na(ari) || ari >= 0.98)
  record(base_row(paste0('site_M', M), 'website', result, M, previous, rho,
                  ari, previous_ari, passed = passed))
}
# Confirm the released default reproduces all earlier experimental ND fits.
for (old in inputs$nd) {
  stopifnot(identical(digest::digest(old$X,algo='sha256'), old$input_hash))
  args <- list(X=old$X, intrinsic_dim=1L, method='isomap', K=50L,
    rw_q=2L, ridge=0, lambda=1, fix_lambda=FALSE, S=NULL,
    position_prior='adaptive', discretization='quantile', iter=10000L,
    num_cores=1L, verbose=FALSE)
  result <- fit_once(old$row$id, old$row$seed+100L, args)
  fit <- result$fit; validate(fit)
  same_length <- length(fit$elbo_trace) == length(old$objective)
  error <- if (same_length) max(abs(fit$elbo_trace-old$objective)) else Inf
  passed <- fit$fit$converged && same_length && error < 1e-7 &&
    isTRUE(all.equal(fit$params, old$params, tolerance=1e-7)) &&
    isTRUE(all.equal(fit$gamma, old$gamma, tolerance=1e-7)) &&
    isTRUE(all.equal(fit$fit$lambda_vec, old$lambda, tolerance=1e-7))
  record(base_row(old$row$id, 'nd_reproduction', result, prefix_error=error, passed=passed))
}
# Original pilot protocol: one Isomap/similarity initialization per candidate.
evaluate <- function(fit, dataset, adaptive) {
  M <- fit$intrinsic_dim
  W <- if (M == 1) matrix(1,ncol(dataset$X),1) else fit$partition$pi_weights
  active <- if (adaptive) which(colMeans(W)>1e-12) else seq_len(M)
  positions <- if (M == 1) matrix(fit$locations$mean$pseudotime,ncol=1) else
    vapply(fit$locations,function(z) z$mean$pseudotime,numeric(nrow(dataset$X)))
  correlations <- abs(cor(dataset$latent_positions,positions[,active,drop=FALSE],method='spearman'))
  correlations[!is.finite(correlations)] <- 0
  size <- max(dim(correlations)); square <- matrix(0,size,size)
  square[seq_len(nrow(correlations)),seq_len(ncol(correlations))] <- correlations
  assignment <- as.integer(clue::solve_LSAP(square,maximum=TRUE))
  rho <- mean(square[cbind(seq_len(ncol(dataset$latent_positions)),assignment[seq_len(ncol(dataset$latent_positions))])])
  list(M=length(active), ARI=mclust::adjustedRandIndex(dataset$true_assign,max.col(W)), recovery=rho)
}
for (old in inputs$pilot) {
  dat <- old$dataset
  stopifnot(identical(digest::digest(dat$X,algo='sha256'), dat$input_hash))
  for (method in c('adaptive','forward')) {
    selected <- NULL
    for (M in if (method == 'adaptive') 8L else 1:8) {
      key <- sprintf('%s_M%02d',method,M)
      result <- fit_once(paste(dat$row$id,key,sep='_'),dat$row$seed+1000000L,
        list(X=dat$X,intrinsic_dim=M,method='isomap',partition_init='similarity',
          partition_prior=if(method=='adaptive') 'adaptive' else 'fixed',
          effective_weight_tol=1e-12,K=50L,rw_q=2L,ridge=0,lambda=1,fix_lambda=FALSE,
          S=NULL,position_prior='adaptive',similarity_metric='spearman',cluster_linkage='single',
          discretization='quantile',num_cores=1L,iter=3000L,max_converge_iter=3000L,
          T_start=5,T_end=1,n_outer=25L,inner_iter=1L,verbose=FALSE))
      fit <- result$fit; validate(fit)
      baseline <- old$candidates[[key]]
      if (is.null(baseline)) stop('New search moved beyond the original baseline; review before release.')
      reference <- baseline$fit$objective
      trace <- fit$elbo_trace
      count <- min(length(trace),length(reference))
      prefix_error <- max(abs(head(trace,count)-head(reference,count)))
      passed <- fit$fit$converged && prefix_error < 1e-6 && length(result$warnings)==0
      row <- base_row(paste(dat$row$id,key,sep='_'),'pilot_candidate',result,
        prefix_error=prefix_error,passed=passed)
      record(row)
      if (!passed) stop('Pilot candidate regression gate failed; inspect saved results.')
      score <- tail(trace,1)
      if (method=='adaptive' || M==1 || score > tail(selected$fit$elbo_trace,1)) {
        selected <- result
      } else break
    }
    metric <- evaluate(selected$fit,dat,method=='adaptive')
    baseline <- old$results[[method]]$row
    passed <- metric$M==baseline$estimated_M && metric$ARI>=baseline$ARI-0.02 &&
      metric$recovery>=baseline$ordering_recovery-0.005
    record(base_row(paste(dat$row$id,method,sep='_'),'pilot_result',selected,
      baseline$estimated_M,baseline$ordering_recovery,metric$recovery,
      metric$ARI,baseline$ARI,passed=passed))
  }
}
results <- do.call(rbind,rows)
stopifnot(all(results$passed),all(results$warnings==0))
writeLines('PASS: all predeclared simulation regression checks passed.',file.path(root,'validation.txt'))
