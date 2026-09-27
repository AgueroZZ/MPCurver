# Homepage convergence cleanup

- Agent: codex
- Update title: homepage-convergence-cleanup
- Update date: 2026-09-26

## Key updates

- Removed the detailed convergence section from README.Rmd, README.md, and the package homepage at the user's request.
- Kept the existing Changelog navigation link and detailed convergence documentation in NEWS.md and the function reference.
- Rebuilt the curated public site with scripts/build_public_site.R and refreshed its search index.
- Updated the site audit to check the simplified homepage and retained Changelog link.

## Validation

- Site build completed successfully.
- Site audit passed: 27 pages, 898 local targets, and 138 search entries; no errors.
- git diff --check passed.
- Documentation-only update; package version remains 0.3.2.
