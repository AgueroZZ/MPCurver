# MPCurver 0.3.2 simulation regression and release

- Agent: Codex
- Update date: 2026-09-26
- Update title: v032-simulation-release
- Authorization: push the package after the earlier simulations and basic
  regression checks show no update-related problems.
- Plan: `plan/2026-09-26-v032-simulation-regression-release.md`.

## Scope and provenance

The new validation bundle is `experiments/v032_convergence_regression/`.
It freezes the earlier website examples, all 27 ND-stopping examples, and
three intrinsic-M pilots (true M = 3, 4, 5; SNR = 4). Original InferOrder
artifacts and unrelated work are preserved. The intrinsic-M protocol uses
one Isomap/similarity initialization per candidate, matching those earlier
pilots rather than changing the search protocol during a regression check.

The previous implementation validation passed 578 assertions, with four
existing retired-feature skips; package check had 0 errors, 0 warnings and
one existing top-level-file NOTE. No fitting implementation has changed
since those checks.

## Publication preparation

- Public package metadata and docs report version 0.3.2.
- The curated site's Changelog item points to the tracked GitHub NEWS.md,
  avoiding an automatic link to an unbuilt pkgdown news page.
- Rebuild with `scripts/build_public_site.R`, then audit with
  `experiments/v032_convergence_regression/audit_site.py`.

## Status

All 47 fresh fits passed: two website examples, 27 ND reproductions, and
18 intrinsic-M candidate fits. All six pilot method outcomes retain their
previous estimated M (adaptive: 4, 4, 5; uniform forward: 3, 4, 5 for true
M = 3, 4, 5). ARI is unchanged, and the largest ordering-recovery decrease
is 0.0000123. Every ND trace and common pilot trace segment matches its
reference exactly. All 18 pilot stopping indices match the first normalized
threshold crossing in the frozen reference traces. Total fresh fitting time
was 1027.7 seconds, serially on one CPU thread.

No fitting call emitted warnings. The initial validation script also computed
unused correlation diagnostics for underfitted candidates, causing three
zero-standard-deviation warnings. That unused computation was removed and
all validation was rerun against the saved new fits, with no further fitting.
The original run log is retained alongside the clean validation recheck log.

The original M = 2 website example still reaches its original iteration cap;
its matched-budget assignments and recovery are unchanged. All 27 ND fits
and all 18 intrinsic-M candidate fits converged. The known adaptive M = 3
split remains a baseline limitation, not a new regression.

See `experiments/v032_convergence_regression/results.md` for the comparison
and `validation.txt` for the passing gate.

## Release gate

The curated site rebuild completed successfully. The audit passed for 27 HTML
pages, 900 local links/assets/anchors and 139 search entries, with no missing
targets or internal-page links. The homepage, fitting/continuation reference
and introductory vignette contain the new version and stopping-rule guidance.
All simulation and documentation gates passed; the release is ready for the
user-authorized push to origin/main.
