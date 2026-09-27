# MPCurver 0.3.2: normalized ELBO convergence

- Agent: Codex
- Date: 2026-09-26
- Request: record the stopping-rule experiment and adopt ND-normalized ELBO
  change as the package default, with updated version and DESCRIPTION.
- Plan: `plan/2026-09-26-normalized-elbo-convergence.md`.

## Experiment recorded

The reproducible study is in the sibling InferOrder repository:
`experiments/elbo_nd_stopping_v031/` (`run.R`, `report.Rmd`, `report.html`,
`results.csv`, `summary.csv`, `validation.json`, saved fits and figures).
It pins MPCurver 0.3.1 at commit
`e38a69a563e338a549382b4b3ece74d655c4654f`.

There were 27 M = 1 datasets: three seeds crossed with N = 100, 200, 400
and D = 10, 20, 40; K = 50, signal standard deviation 1, noise standard
deviation 0.5, Isomap initialization, RW2 with ridge 0, adaptive position
prior, and estimated noise and smoothness. Data were nested subsets of
each seed's master matrix. Each dataset received two new fits, with identical
data and initialization and counterbalanced rule order, plus an existing
stricter reference fit.

| Rule | Tolerance | Fits | Median iterations | Total seconds |
| --- | ---: | ---: | ---: | ---: |
| Relative reference | 1e-8 | 27 | 3884 | 334.024 |
| Relative | 1e-6 | 27 | 633 | 50.753 |
| ND-normalized | 1e-6 | 27 | 591 | 46.632 |

At the same numerical tolerance, normalized stopping took 8.1% less total
time, used fewer iterations in 27/27 cases, and took less time in 23/27 cases.
All 54 new fits converged; their traces matched prefixes of the reference
traces, and their stopping indices matched the first eligible threshold
crossing. Input hashes were also verified.

Relative to the stricter reference, normalized fits had minimum ordering
Spearman correlation 0.9999625, median prediction RMSE difference 0.0090272,
and maximum prediction RMSE difference 0.0190425. Mean remaining ELBO
difference per entry was 0.00049015 (maximum 0.00129517). Prediction RMSE
against truth averaged 0.132617, versus 0.132495 for the reference.

The normalized rule was modestly faster at equal tolerance; the much larger
gain over the reference mostly reflects relaxing 1e-8 to 1e-6. These simple
single-ordering examples do not establish universal runtime gains or
multi-ordering model-selection accuracy. A small increment is not a bound
on the remaining optimization gap.

## Package changes

- Version 0.3.2; DESCRIPTION, NEWS, public help, README and introductory
  vignette updated.
- New fits use `convergence = "normalized"`: `abs(delta_ELBO) / (N * D)`.
  Both `tol` and `tol_outer` default to 1e-6; zero disables early stopping.
- Single-ordering known-noise and estimated-noise paths and structural
  multi-ordering fits use the same normalization. Existing numerical decrease
  guards and post-annealing T = 1 stopping are retained.
- `convergence = "relative"` preserves the old denominators. For the previous
  defaults also use `tol = 1e-6, tol_outer = 1e-5`.
- Continuation inherits the stored rule. Objects without a rule use relative
  stopping unless explicitly overridden. ELBO traces keep their original scale.
- The rule passes through public fits, multiple initialization methods,
  dimension selection and the two-ordering compatibility wrapper. Preliminary
  partition initialization fits retain their existing relative stopping so
  changing the final stopping rule does not also change the starting state.
- Frozen 0.3.1 experiment libraries and results are preserved.

## Validation

- Full test suite: 578 passed assertions, 0 failures, 0 errors, 0 warnings;
  four existing skips for retired greedy/compaction behavior. The new
  convergence tests account for 80 passed assertions.
- Tests cover first-threshold crossing against independently calculated trace
  increments; known and estimated noise; final T = 1 partition stopping;
  continuation and legacy-object fallback; and all dimension-selection
  candidates receiving the requested rule.
- Raw test results: `log/2026-09-26-convergence-test-results.rds`.
- README rendered and roxygen help regenerated.
- `devtools::check(document = FALSE, args = "--no-manual")`: 0 errors,
  0 warnings, one pre-existing top-level-file NOTE for `README.Rmd` and
  `elife-61271-fig2-data1-v2.csv`. Examples, tests, vignette build and vignette
  rebuild passed. Initial sandbox execution could not start processx; the
  check completed successfully with approved execution outside the sandbox.
- Check evidence: `log/2026-09-26-convergence-package-check.txt` and
  `log/2026-09-26-convergence-check-results.rds`.
- Installed the current checkout using `R CMD INSTALL` into the normal R user
  library, replacing installed 0.3.0. A fresh R process verified version 0.3.2,
  both 1e-6 defaults, and normalized controls in actual M = 1 and M = 2 fits.
  Evidence: `log/2026-09-26-convergence-installed-verification.txt`.
- No GitHub push or website deployment was performed in this update.
