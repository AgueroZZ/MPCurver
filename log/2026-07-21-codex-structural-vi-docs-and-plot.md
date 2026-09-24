# Structural VI documentation and trajectory plotting

- **Agent:** codex
- **Update title:** structural-vi-docs-and-plot
- **Update date:** 2026-07-21

## Key updates

- Implemented a genuine fixed-`M` structural `plot_type = "mu"` path that
  reads the canonical `conditional_posterior$mean` state, supports one- and
  two-feature trajectories, annotates `q(Z_j = m)`, and uses an adaptive panel
  layout for larger `M`.
- Corrected the single-ordering plot delegation so the public `dims` argument
  is honored for `plot_type = "mu"`.
- Rewrote all three current mathematical derivations to match the implemented
  single-ordering and structural variational families, including intrinsic
  pseudo-ELBO conventions, known `S`, shared estimated `sigma2`, lambda
  updates, null-space identifiability, annealing, and the fixed-`M` scope.
- Expanded and regenerated public help for fitting, continuation, plotting,
  summaries, fitted priors, and the complete `mpcurve` object schema. The
  public object page is now included in the curated pkgdown reference index.
- Updated the README, partition vignette, fitness vignette, and public site to
  describe fixed-`M` structural VI and to demonstrate canonical structural
  trajectory plots for estimated and known measurement noise.
- Updated `scripts/build_public_site.R` to load the current checkout before
  executing vignette code. This prevents a stale installed package from
  silently generating outdated figures while the site prose uses new sources.
- Versioned editable `important_derivations/*.Rmd` sources while continuing to
  ignore rendered/internal assets. Added a filename-level ignore rule for the
  temporary Chinese derivation PDF and removed both repository-local copies.

## Additional relevant correction

- The previous README implied that a vector of initialization methods would
  all run with `num_cores = 1`, but the implementation runs only the first in
  that case. The example now uses an explicit serial loop, and the help page
  documents both serial and parallel behavior, including `NULL` entries for
  failed parallel fits.

## Main files

- `R/07_mpcurve.R`
- `R/MPCurver-package.R`
- `important_derivations/cavi_math_details.Rmd`
- `important_derivations/known_noise_cavi_derivation.Rmd`
- `important_derivations/partition_coherent_elbo_derivation.Rmd`
- `README.Rmd`, `README.md`
- `vignettes/partition.Rmd`, `vignettes/fitness.rmd`
- `scripts/build_public_site.R`, `_pkgdown.yml`, `docs/`
- `tests/testthat/test-partition-m-cavi.R`

## Verification

- Focused structural plot/API/analytic tests passed with no failures or
  warnings.
- Full testthat suite: 413 passed, 0 failed, 0 warnings, and 4 intentional
  skips for retired greedy/active-compaction semantics.
- All 37 generated Rd files passed `tools::checkRd()`.
- All three mathematical Rmd files rendered successfully; focused
  single-ordering and structural analytic suites passed 88/88 and 42/42.
- The curated public site built successfully with
  `Rscript scripts/build_public_site.R`; only `mpcurve_intro`, `partition`, and
  `fitness` are public articles, and the internal-material leak scan was clean.
- Structural 1D, 2D, and known-`S` trajectory figures were visually inspected
  after rendering against the current checkout.
- `R CMD build .` succeeded and
  `R CMD check --no-manual MPCurver_0.3.0.tar.gz` returned `Status: OK`.
  Repository-index access warnings were environmental and did not change the
  final check status.

## Follow-up

- A separate predictive or explicitly complexity-penalized criterion is still
  required before reintroducing automatic cross-`M` selection.
- Deprecated freeze/drop formals remain warning-only compatibility arguments
  for the 0.3.x transition and can be removed in a later breaking release.
