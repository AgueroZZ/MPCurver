# Reader-facing documentation and website

- Agent: codex
- Update title: reader-facing-documentation
- Update date: 2026-09-24

## Key updates

- Rewrote the home page and package overview around sample orderings, smooth
  feature trajectories, uncertainty, and practical analysis workflows.
- Reorganized the feature-partitioning tutorial around a worked analysis,
  assignment interpretation, result extraction, and convergence, followed by
  the statistical model and structural variational updates.
- Revised the introductory and fitness tutorials to explain modeling choices
  and results directly. Clarified intrinsic-prior notation and the fitness
  example's measurement-error calibration assumption.
- Rewrote public fitting, plotting, summary, and object help to explain what
  users obtain and how to interpret it. Confined compatibility details to the
  relevant parameter and object entries.
- Regenerated README, Rd help, and the curated public site using
  `Rscript scripts/build_public_site.R`.

## Validation

- All executable R expressions in `R/` are unchanged from the preceding commit.
- All fitting and plotting chunks retain their original settings. The partition
  tutorial adds an assignment-probability display and revises the displayed
  convergence diagnostic.
- The complete testthat suite passed, with four existing skips for superseded
  greedy/active-state behavior.
- All 37 Rd files passed `tools::checkRd()`.
- README and all three public tutorials rendered successfully; generated figures
  are unchanged.
- Checked all 26 generated HTML pages: no broken local links or assets, and no
  references to excluded internal documents or directories.
- `git diff --check` passed.

## Scope

This is an editorial update to version 0.3.0, not an algorithm change. The
implementation plan is `plan/2026-09-24-reader-facing-documentation.md`.
Experiments and their local notes are not part of this publication.
