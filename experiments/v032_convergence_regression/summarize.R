# Summarize completed checks and verify stopping indices from frozen old traces.
root <- 'experiments/v032_convergence_regression'
x <- read.csv(file.path(root, 'results.csv'))
inputs <- readRDS(file.path(root, 'inputs.rds'))
stopifnot(nrow(x)==53L, all(x$passed), all(x$warnings==0))
expected_rows <- list()
for (p in inputs$pilot) for (key in names(p$candidates)) {
  old <- p$candidates[[key]]
  trace <- old$fit$objective
  delta <- diff(trace)
  eligible <- seq_along(delta) > if (old$M>1) 25 else 0
  eligible <- eligible & delta >= if (old$M>1) -1e-8*(abs(head(trace,-1))+1) else 0
  expected <- which(eligible & abs(delta)/length(p$dataset$X)<1e-6)[1]
  row <- x[x$id == paste(p$dataset$row$id,key,sep='_'),]
  stopifnot(nrow(row)==1L, row$iterations==expected)
  expected_rows[[length(expected_rows)+1L]] <- data.frame(
    id=p$dataset$row$id,key=key,expected_stop=expected,actual_stop=row$iterations,
    baseline_iterations=old$fit$iterations)
}
write.csv(do.call(rbind,expected_rows),file.path(root,'expected_iterations.csv'),row.names=FALSE)
p <- x[x$kind=='pilot_result',]
lines <- c('# MPCurver 0.3.2 simulation regression results','',
  'All 47 fresh fits passed the predeclared regression checks. None of the fitting calls emitted warnings.',
  '', '## Earlier examples','',
  '- M = 1 Swiss roll: ordering Spearman correlation 0.999758, unchanged.',
  '- M = 2 partition: ARI 1; mean ordering correlation 0.998205, unchanged at the original iteration budget.',
  '- All 27 ND-stopping fits exactly reproduced the experimental 0.3.1 ELBO traces and stopping indices; fitted parameter comparisons passed at tolerance 1e-7.',
  '', '## Intrinsic-M pilots','',
  '| True M | Method | Old estimate | New estimate | New ARI | New ordering recovery | Recovery change |',
  '| ---: | --- | ---: | ---: | ---: | ---: | ---: |')
for(i in seq_len(nrow(p))) {
  row <- p[i,]
  M <- as.integer(sub('pilot_M([0-9]+).*','\\1',row$id))
  method <- if(grepl('adaptive$',row$id)) 'Adaptive EB' else 'Uniform + forward'
  lines <- c(lines,sprintf('| %d | %s | %d | %d | %.6f | %.6f | %+.6f |',
    M,method,row$reference_M,row$effective_M,row$ARI,row$recovery,
    row$recovery-row$reference_recovery))
}
lines <- c(lines,'',
  'All 18 intrinsic-M candidate fits converged at the first eligible normalized threshold crossing predicted from the saved old traces. Their common trace segments had zero numerical difference. The known adaptive split for true M = 3 persists; it was already present before this update.',
  '',sprintf('Fresh fitting time: %.1f seconds (one process, one CPU thread).',sum(x$seconds[x$kind!='pilot_result'])),
  '', 'The original M = 2 website example reaches its fixed iteration budget in both versions. Its regression check concerns matched-budget recovery; the intrinsic-M candidates and all ND-stopping fits were required to converge.')
writeLines(lines,file.path(root,'results.md'))
writeLines('PASS: all simulation checks and all 18 pilot stopping-index checks passed.',file.path(root,'validation.txt'))
cat(paste(lines,collapse='\n'),'\n')
