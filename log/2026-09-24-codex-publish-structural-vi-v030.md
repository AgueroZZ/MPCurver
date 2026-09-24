# MPCurver 0.3.0 Structural VI Publication

- **Agent:** codex
- **Update title:** publish-structural-vi-v030
- **Update date:** 2026-09-24

## Key updates

- Prepared the fixed-`M` structural VI implementation, unified package API,
  generated help pages, public vignettes, and mathematical derivations for
  publication.
- Rebuilt the curated pkgdown site with `scripts/build_public_site.R`; the public
  surface contains only the home page, curated reference pages, and the
  `mpcurve_intro`, `partition`, and `fitness` articles.
- Kept all local experiments and results, local package libraries, temporary
  PDFs, generated output, and experiment-specific plans/logs out of the
  release.

## Validation

- `devtools::document()`: completed.
- `devtools::test()`: passed; four expected skips cover retired greedy cross-`M`
  and active/compaction semantics.
- `R CMD build --no-manual`: completed for `MPCurver_0.3.0.tar.gz`.
- `R CMD check --no-manual`: `Status: OK`.
- Curated-site audit: 26 HTML files checked; all local links and assets resolve,
  and no internal project materials are exposed.
- Real-browser deployment QA found raw Rd math delimiters on reference pages.
  Setting `template.math-rendering: katex` in `_pkgdown.yml` and rebuilding the
  curated site supplies the required renderer for inline and display equations.

## Publication

- Main release commit: `d6a104659d972dd221bd2d328d4f6998c9558603`
  (`Release MPCurver v0.3.0 structural VI`).
- The release commit was pushed to `origin/main`; its first GitHub Pages build
  completed successfully. A follow-up display-only commit adds KaTeX rendering
  after the browser QA finding above.
- KaTeX follow-up commit: `7f5f635445b94eb2b157233b5fc9d06c3fbd539a`
  (`Fix pkgdown math rendering`).
- GitHub Pages workflow run `36057414794` completed successfully. Final browser
  QA confirmed version `0.3.0`, rendered structural-VI display and inline math,
  the fixed-`M` partition article, and the structural `plot_type = "mu"` help.
