# Publish Structural VI v0.3.0 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Publish the completed MPCurver 0.3.0 structural-VI implementation,
documentation, and curated pkgdown website to GitHub while keeping local
experiments, experimental environments, and temporary reference material out
of the public repository update.

**Architecture:** Regenerate package metadata and the curated `docs/` site from the current checkout, validate the package and public-site surface, then stage an explicit allowlist of release files. Commit and push the release directly to `main`, whose committed `docs/` directory is the GitHub Pages source, and verify both the remote repository and deployed website.

**Tech Stack:** R, roxygen2, testthat, pkgdown, Git, GitHub Pages, shell-based release validation.

---

### Task 1: Regenerate package documentation

**Files:**
- Modify: `NAMESPACE`
- Modify: `man/*.Rd`

- [x] Run `Rscript --vanilla -e 'devtools::document()'` from the repository root.
- [x] Confirm `DESCRIPTION` remains version `0.3.0` and `R/10_partition_structural_cavi.R` remains the fixed-`M` implementation selected by the public wrapper.
- [x] Run `git diff --check` and stop on whitespace errors.

### Task 2: Validate the package

**Files:**
- Test: `tests/testthat/*.R`

- [x] Run `Rscript --vanilla -e 'devtools::test(reporter = "summary", stop_on_failure = TRUE)'`.
- [x] Build a source tarball in a temporary directory with `R CMD build`.
- [x] Run `R CMD check --no-manual` on the source tarball and review every warning, note, and error before publication.

### Task 3: Rebuild and audit the public website

**Files:**
- Modify: `docs/**`
- Use: `scripts/build_public_site.R`

- [x] Run `Rscript --vanilla scripts/build_public_site.R`; do not call `pkgdown::build_site()`.
- [x] Verify the home page and public reference/article pages report version `0.3.0` and describe fixed-`M` structural VI, shared `sigma2`, intrinsic `ridge = 0`, canonical state, and retired freeze/drop controls.
- [x] Verify committed `docs/` does not contain or reference `AGENTS`, `CLAUDE`, `PACKAGE_OVERVIEW_FOR_AGENTS`, `Intro`, `internal/`, `important_derivations`, `log/`, or `plan/`.
- [x] Check every local link and asset reference in the curated HTML site.

### Task 4: Stage the release allowlist

**Files:**
- Stage package source and metadata: `.Rbuildignore`, `.gitignore`, `DESCRIPTION`, `NAMESPACE`, `R/`, `README.Rmd`, `README.md`, `_pkgdown.yml`, `man/`, `scripts/build_public_site.R`, `tests/testthat/`, `vignettes/`
- Stage the curated site: `docs/`
- Stage internal mathematical sources: `important_derivations/*.Rmd`
- Stage relevant project records: the July structural-VI implementation and
  documentation plans/logs, plus this publication plan/log
- Exclude: `experiments/`, all local `.library/` directories, temporary PDFs,
  `output/`, experiment-specific plans/logs, and unrelated August PFPCA
  plans/logs

- [x] Use explicit `git add` paths rather than `git add .`.
- [x] Inspect `git diff --cached --stat`, `git diff --cached --name-status`, and the complete staged path list.
- [x] Confirm no PDF, local package library, temporary artifact, or unrelated PFPCA experiment is staged.

### Task 5: Commit, push, and verify deployment

**Files:**
- Create: `log/2026-09-24-codex-publish-structural-vi-v030.md`

- [x] Record package checks, site checks, and staged scope in the update log; add the commit and deployment verification after publication.
- [ ] Commit the staged release as `Release MPCurver v0.3.0 structural VI`.
- [ ] Push `main` to `origin`.
- [ ] Verify GitHub `main` reports version `0.3.0` and contains `R/10_partition_structural_cavi.R`.
- [ ] Poll the GitHub Pages URL until it reports version `0.3.0`, then verify the online structural partition article and primary help pages.
