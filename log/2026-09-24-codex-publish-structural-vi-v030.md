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

## Publication

- Commit: pending.
- GitHub Pages deployment: pending.
