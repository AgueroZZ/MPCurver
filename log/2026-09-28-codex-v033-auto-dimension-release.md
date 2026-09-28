# MPCurver 0.3.3 automatic-dimension release

- Agent: Codex
- Date: 2026-09-28 UTC
- Request: validate automatic ordering-count initialization on the previous
  ordering-count simulation, publish the research update, and release the
  software if the implementation and results are sound.
- Plan: `plan/2026-09-28-v033-auto-dimension-release.md`.

## Software behavior

`fit_mpcurve(intrinsic_dim = "auto")` now computes the feature similarity
once, uses single-linkage clustering by default, and chooses the eligible cut
with the largest mean silhouette. The default similarity is the fast
fixed-df natural-cubic-spline variance-explained score with
`spline_r2_df = 5`. Every selected cluster must contain at least
`similarity_min_cluster_size = 2` features by default. If no multi-ordering
cut is eligible, the fit falls back to one ordering.

For an adaptive partition prior, every feature starts with the same global
ordering-probability vector given by the selected cluster sizes divided by the
number of features. The tree cut is not imposed as a hard feature assignment.
Explicit numeric dimensions preserve their previous behavior, including the
Spearman default, uniform adaptive initialization, and possible singleton
initialization clusters. The current model does not classify singleton or
unstructured features as garbage.

## Simulation pilot

Before the version-only and documentation changes, the current implementation
was installed in an isolated library and run on replicate 1 from all nine
true-M by SNR conditions of the one-monotone-anchor study. All nine automatic
cuts selected the true M, all nine adaptive fits retained the true effective
M, every feature partition had ARI 1, every fit converged, and no warning was
emitted. The prior failure `main_M5_S1_r001` changed from effective M = 1 under
the original eight-slot adaptive fit to selected and effective M = 5, with ARI
1. Its ordering score remained 0.510, reflecting the low-SNR position-recovery
limit rather than a grouping failure. The formal 90-dataset extension uses the
frozen release archive and is recorded in the InferOrder repository.

## Package validation

- Full test suite: 638 successful expectations in 85 test blocks, zero
  failures/errors/warnings, and four existing skips for retired greedy or
  compaction semantics. Saved results are in
  `experiments/v033_auto_dimension_release/test-results.rds`.
- `R CMD check --no-manual`: zero errors, zero warnings, one environment NOTE
  because the optional Suggests `devtools` and `fields` were unavailable.
  Examples, packaged tests, vignette execution, and vignette rebuild all
  passed. The check log is stored beside the test results.
- Source archive: `MPCurver_0.3.3.tar.gz`, SHA-256
  `cb57e1f6e8d859b7ff9aeda17c97d538fae0686eb0ceb616f126023a26d645f9`.
- Curated pkgdown site: 27 pages, 903 local links/assets/anchors, and 138
  search entries checked with no errors. The build used pkgdown 2.2.0 and
  fontawesome 0.5.3, matching the published site's dependency family.
- `git diff --check` passed after normalizing one generated trailing blank
  line. Rendered README and partition figures were visually inspected.

The release was prepared with R 4.3.3, testthat 3.3.2, roxygen2 8.1.0,
rmarkdown 2.26, and knitr 1.52. Publication verification is appended after the
authorized push.

## Publication verification

The release commit is `c905901424e43eab78b58bdcc0d1de367ec8fd73` on
`origin/main`. GitHub Pages deployment
[36487823465](https://github.com/AgueroZZ/MPCurver/actions/runs/36487823465)
completed successfully for that exact commit. The live `fit_mpcurve()`
reference page exposes `intrinsic_dim = "auto"`, `spline_r2_df = 5`, and
`similarity_min_cluster_size = 2`.
