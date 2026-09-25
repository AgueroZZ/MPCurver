# Effective M and uniform-prior greedy selection

- **Agent:** codex
- **Update title:** effective-m-uniform-greedy
- **Update date:** 2026-09-25

## Key updates

- Every `mpcurve` fit reports `effective_intrinsic_dim`. For an adaptive partition prior, it counts fitted ordering weights strictly above the configurable `effective_weight_tol` (default `1e-12`). For a fixed partition prior or a single-ordering fit, it equals the fitted `intrinsic_dim`.
- Added `select_mpcurve_dimension()` as a separate forward/backward selector under a fixed uniform partition prior. It compares final soft `T = 1` structural objectives, tries similarity and ordering-method initialization for each `M >= 2`, and returns candidate scores, convergence states, decisions, and the selected fit.
- `fit_mpcurve(greedy = ...)` remains disabled; its message now points to the separate selector. Continuation preserves the selection record and marks that fitting occurred after selection.
- Updated the package help pages, README, partition vignette, reference index, public API test, and focused tests.

## Verification

- `devtools::test()`: 498 passed, 4 existing legacy tests skipped, no failures or test warnings.
- `devtools::check(document = FALSE, args = "--no-manual", error_on = "never")`: 0 errors, 0 warnings, 1 pre-existing top-level-file note for `README.Rmd` and `elife-61271-fig2-data1-v2.csv`.
- `git diff --check`: clean.

The greedy result is a comparison of fitted variational solutions; convergence and initialization can affect the selected dimension. The adaptive effective count is a thresholded occupancy diagnostic, not a posterior over `M`.
