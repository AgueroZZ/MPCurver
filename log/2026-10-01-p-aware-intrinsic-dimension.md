# Feature-count-aware intrinsic dimension — 2026-10-01

The user selected a numerical occupancy rule based on total soft feature
assignments: count ordering m when N_m = sum_j Pr(Z_j = m) > 1e-8. P denotes
ncol(X), the number of feature columns. This replaces a fixed threshold on
ordering proportions. The strict inequality excludes a count exactly at the
threshold; no minimum of one feature and no hard-assignment rule is imposed.

Under the adaptive EB update, omega_m = N_m / P. Therefore the equivalent
prior-weight cutoff is effective_count_tol / P. The implementation sums the
posterior feature-assignment probabilities directly. It does not rely on
cached prior weights or allocate additional data matrices. The advanced control
is renamed from effective_weight_tol to effective_count_tol, with default 1e-8
and units of expected feature count. The tolerance must be nonnegative and below
P / num_orderings so at least one ordering exceeds it in a row-stochastic fit.

The fixed-prior convention remains unchanged: intrinsic_dim reports the fitted
ordering count. num_orderings continues to record retained slots in both modes.
No slots, trajectories, or assignments are removed by the reporting rule.
Continuation recalculates the estimate from its updated assignments. Explicit
fraction thresholds in older saved fits are translated into count units by
multiplying by P, preserving those saved settings; new fits use the count rule.

README, the partition tutorial, object documentation, control reference, and
release notes describe the formula and units. The old reporting log documents
the earlier definition; this change supersedes its threshold rule.

## Verification

- The full local suite has 824 passing expectations, zero failures, errors,
  warnings, or skips. Coverage includes scaling the same tiny ordering
  proportion by different P, strict-boundary behavior, invalid counts, the EB
  relation between prior proportions and expected counts, legacy setting
  conversion, continuation, and plotting when reported dimension is smaller
  than retained slots.
- Before/after one- and three-ordering fits have exactly identical underlying
  numerical state, parameters, responsibilities, and locations after removing
  the renamed reporting control from the comparison.
- Updated README and partition tutorial execute through the documented
  renderers; roxygen and affected references are regenerated. Site audit passes
  for 30 pages and 996 local links, assets, and anchors. `git diff --check`
  passes. Actual Firefox rendering of the revised formula, units, and
  example expected counts is inspected; review images remain under `/tmp/`.
  The earlier full package check predates this threshold change and
  is historical evidence; it is not reported as a current check.

All changes remain local. No commit, push, or release is performed.
