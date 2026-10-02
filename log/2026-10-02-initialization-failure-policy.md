# Explicit initialization failure handling

Review item R3 is fixed according to the user's decision: invalid initializer
settings must error, numerical failures may fall back to PCA by default with
a warning, and fallback must be explicit.

## Behavior

- `mpcurve_init_control(on_failure = "pca")` is the default. An ordering
  computation failure warns with the requested method, original cause, group
  label where applicable, and the option for disabling fallback. It retries
  with PCA component 1.
- `on_failure = "error"` stops on the original failure. A failed PCA attempt
  always stops. The previous last-resort first-feature-rank fallback is removed.
- Shared validation checks helper settings before the fallback handler,
  including sample-count-dependent bounds, embedding-coordinate bounds, and
  integer PCA indices. Fractional indices are no longer truncated.
- Only ordering computation errors are caught. Errors in subsequent model
  fitting propagate instead of being reinterpreted as initializer failures.
- Single-feature groups run the requested helper, honoring its options,
  instead of silently replacing it with a first-feature ranking.
- Single, grouped, and automatic fits use the same failure policy. Initialization
  metadata records the requested/actual method, PCA component, and fallback
  cause. Continuation preserves this metadata without rerunning initialization.

The implementation updates `R/12_interface.R`, `R/13_helpers.R`,
`R/10_partition_cavi.R`, `R/10_partition_structural_cavi.R`, `R/09_cavi.R`, and
the internal initializer in `R/06_cSmoothEM.R`. The public fitting signature
remains unchanged; the new option belongs to `init_control`.

## Current verification

- Full source suite: **1,121 passing expectations**, zero failures, errors,
  unexpected warnings, or skips. This includes 71 new regression expectations
  in `tests/testthat/test-initialization-failure-policy.R`, covering invalid
  settings, warned fallback, strict mode, failed PCA, downstream errors, and
  continuation provenance. Evidence: `/tmp/mpcurver-r3-tests.R`, `.log`, `.rds`.
- An independent public fit with duplicated observations triggers an actual
  Rtsne initialization error, emits the PCA fallback warning, and returns a fit
  recording requested `tSNE`, actual `PCA`, and the original cause. Strict mode
  stops on that error. Evidence: `/tmp/mpcurver-r3-evidence.R` and `.log`.
- README executes. The initialization-control reference, home page, news,
  search index, and sitemap regenerate from canonical sources. Roxygen 7.3.2
  emits the existing version advisory for the recorded 7.3.3 RoxygenNote;
  generation completes using the pinned pkgdown tools.
- Actual Firefox checks pass for home, initialization-control reference, and
  news: the new option and warning policy appear, no math errors occur, and
  no page-width overflow occurs at a 1,053-pixel viewport. The rendered argument
  description screenshot was inspected. Evidence: `/tmp/mpcurver-r3-browser.py`,
  `.log`, and `/tmp/mpcurver-r3-on-failure.png`.
- Site audit passes **33 HTML pages and 1,131 local links/assets/anchors**,
  including existing figure/table numbering checks. `git diff --check` passes.

The preceding full `R CMD check` remains historical; it was not repeated for
this change. Existing user changes and historical InferOrder results are
preserved. No commit or push.
