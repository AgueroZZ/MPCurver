# MPCurver 0.3.1 publication

- **Agent:** codex
- **Update title:** publish-v031
- **Update date:** 2026-09-25

## Key updates

- Published `effective_intrinsic_dim` reporting for adaptive partition priors and `select_mpcurve_dimension()` for fixed uniform-prior forward/backward ELBO selection.
- Updated `DESCRIPTION` to 0.3.1 and rebuilt the curated pkgdown website, including the new selector reference page and revised partition article.
- Kept experimental results and unrelated untracked plans/logs outside the release commit.

## Validation and deployment

- `devtools::test(reporter = "summary", stop_on_failure = TRUE)`: passed with four expected legacy skips.
- `devtools::check(document = FALSE, args = "--no-manual")`: 0 errors, 0 warnings, and one existing top-level-file note for `README.Rmd` and `elife-61271-fig2-data1-v2.csv`.
- Curated site audit: 27 HTML pages, 583 local links/assets checked, no missing targets or internal-content references.
- Release commit `df0fd82a1b136b09d259154c6bea2d51aa0fada5` was pushed to `origin/main` and verified as the remote head. GitHub Pages run `36171600649` completed successfully.
- Browser verification confirmed the live home page reports 0.3.1, the new `select_mpcurve_dimension` reference page loads, and the partition article displays the effective-M guidance and rendered mathematical expression.
