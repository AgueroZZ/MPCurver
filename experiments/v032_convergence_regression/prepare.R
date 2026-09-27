# Freeze earlier simulation inputs and reference summaries without changing them.
origin <- normalizePath('../InferOrder/experiments')
out <- 'experiments/v032_convergence_regression'
site <- lapply(1:2, function(M) {
  z <- readRDS(file.path(origin, 'mpcurve_v030_site/results/simulations',
                         sprintf('simulation_m%d.rds', M)))
  list(simulation = z$simulation, summary = z$summary,
       objective = z$fit$elbo_trace)
})
nd_root <- file.path(origin, 'elbo_nd_stopping_v031')
manifest <- read.csv(file.path(nd_root, 'manifest.csv'))
nd <- lapply(seq_len(nrow(manifest)), function(i) {
  row <- manifest[i, ]
  z <- readRDS(file.path(nd_root, 'fits', paste0(row$id, '_nd.rds')))
  list(row = row, X = z$fit$data, input_hash = z$input_hash,
       objective = z$fit$elbo_trace, params = z$fit$params,
       gamma = z$fit$gamma, lambda = z$fit$fit$lambda_vec)
})
pilot_root <- file.path(origin, 'estimate_intrinsic_m_v031')
pilot <- lapply(3:5, function(M) {
  id <- sprintf('pilot_M%d_S4', M)
  candidates <- list()
  for (method in c('adaptive', 'forward')) {
    for (m in if (method == 'adaptive') 8L else seq_len(M + 1L)) {
      key <- sprintf('%s_M%02d', method, m)
      candidates[[key]] <- readRDS(file.path(pilot_root, 'candidates',
        paste0(id, '_', key, '.rds')))
    }
  }
  results <- lapply(c('adaptive', 'forward'), function(method)
    readRDS(file.path(pilot_root, 'results', paste0(id, '_', method, '.rds'))))
  names(results) <- c('adaptive', 'forward')
  list(dataset = readRDS(file.path(pilot_root, 'data', paste0(id, '.rds'))),
       candidates = candidates, results = results)
})
saveRDS(list(site = site, nd = nd, pilot = pilot), file.path(out, 'inputs.rds'),
        compress = 'xz')
writeLines(c('Reference sources: InferOrder experiments/mpcurve_v030_site,',
  'elbo_nd_stopping_v031 and estimate_intrinsic_m_v031.',
  'ND and intrinsic-M references: MPCurver 0.3.1 commit',
  'e38a69a563e338a549382b4b3ece74d655c4654f.',
  paste('Input bundle SHA256:', digest::digest(file=file.path(out,'inputs.rds'),algo='sha256'))),
  file.path(out, 'provenance.txt'))
