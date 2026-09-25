# Effective intrinsic dimension and uniform-prior greedy selection

**Goal:** Add a transparent effective-ordering count to fitted `mpcurve` objects and a separate public selector for the uniform feature-assignment prior.

**Confirmed API choices:** Keep the current default `partition_prior = "adaptive"`. Whenever the fitted partition prior is adaptive, count ordering weights above a configurable default threshold of `1e-12`. A fixed prior retains the fitted `M` as its effective dimension. Keep `intrinsic_dim`, `active_intrinsic_dim`, and the fitted trajectories unchanged.

## Implementation

- [x] Add `effective_intrinsic_dim` to `mpcurve` results and summaries. Store the threshold in fit control and retain the raw estimated `omega` through the existing fitted-prior interface. Ensure `do_mpcurve()` continuation recalculates the count.
- [x] Add exported `select_mpcurve_dimension(X, max_intrinsic_dim, direction, ...)` that orchestrates fixed-`M` fits under a uniform feature-assignment prior. Compare the best finite final soft `T = 1` objective across similarity and ordering-method initializations at each candidate `M`; accept an adjacent step only for a strictly higher objective. Return the selected `mpcurve` plus a search record with candidate scores, convergence status, score differences, and stop reason.
- [x] Keep the old `fit_mpcurve(greedy = ...)` path disabled; its historical helpers depend on obsolete active-ordering and compaction behavior. Point the error message to the new selector.
- [x] Update roxygen-generated documentation, curated public export test, and relevant user-facing package guidance. Record the change in `log/` using the required update-log fields.

## Verification

- [x] Test fixed/adaptive effective counts, numerical-tail behavior, threshold validation, and continuation behavior.
- [x] Test both greedy directions, output path accounting, fixed-uniform prior validation, and score matching against individually fitted candidates on small data.
- [x] Run targeted testthat files, then the applicable package checks and documentation generation. Inspect changed files and generated docs for unintended edits.

The selector is exploratory: a greedy stop can depend on local optima and fit convergence. Candidate diagnostics must remain visible so that the selected dimension is reviewable.
