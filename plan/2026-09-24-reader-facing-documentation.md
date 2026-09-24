# Reader-facing MPCurver documentation

## Goal

Explain the statistical model, analysis workflow, and interpretation of results
to researchers and software users. Replace development-status and migration
narratives in introductory material with substantive explanations of capabilities.

## Scope

- Revise `README.Rmd`, `DESCRIPTION`, and package/object documentation in
  `R/MPCurver-package.R`.
- Revise public roxygen documentation in `R/07_mpcurve.R`; retain accurate,
  concise descriptions of supported and deprecated arguments in the reference.
- Revise `vignettes/mpcurve_intro.Rmd`, `vignettes/partition.Rmd`, and
  `vignettes/fitness.rmd`. Lead the partition tutorial with a worked analysis,
  then explain the model and variational updates.
- Regenerate README, Rd help, and the curated site using
  `scripts/build_public_site.R`.
- Verify that executable package code is unchanged, render all public examples,
  check links and public-site exclusions, and inspect the deployed pages.
- Commit and publish the documentation update with a project update log.

## Editorial principles

- State what is estimated and how users interpret it.
- Explain structural VI through dependence between assignments and trajectories.
- Keep mathematical qualifications close to the quantities they describe.
- Put compatibility details in parameter/object reference entries; keep migration
  history and implementation bookkeeping out of tutorials and the home page.
- Preserve the scientific examples and their fitting settings.
