# Structural VI Documentation and Plot Implementation Plan

> **For agentic workers:** Execute each checkbox in order. Keep all code,
> comments, and collaborator-facing documentation in English. Do not commit or
> publish changes without an explicit user request.

**Goal:** Make every current package-facing explanation, mathematical
derivation, help page, vignette, and plot method consistent with the fixed-M
structural variational implementation, including a genuine structural
`plot_type = "mu"` visualization.

**Architecture:** Treat the canonical structural state in
`R/10_partition_structural_cavi.R` as the sole implementation authority.
Partition plotting reads `conditional_posterior$mean`, `gamma`, and
`partition$pi_weights` from the normalized `mpcurve` object; help pages and
derivations describe exactly those fields and the factorization
`q(C) prod_j q(Z_j) q(U_j | Z_j)`. Historical augmented mean-field material is
removed from current derivations rather than mixed into the implementation
description.

**Tech Stack:** R, base graphics, roxygen2, testthat, R Markdown, pkgdown via
`scripts/build_public_site.R`, pypdf for read-only PDF reference extraction.

---

### Task 1: Implement a genuine structural trajectory plot

**Files:**
- Modify: `R/07_mpcurve.R`
- Modify: `tests/testthat/test-partition-m-cavi.R`
- Modify: `tests/testthat/_snaps/public-api.md` only if the documented public
  signature changes (it should not)

- [x] Add a test that fits a small fixed-M structural model, opens a null PDF
  device, and verifies that both one-dimensional and two-dimensional
  `plot(fit, plot_type = "mu", dims = ...)` calls return invisibly without
  changing the fit.
- [x] Add an internal plotting helper which validates the canonical structural
  posterior and draws one panel per fixed ordering. For two dimensions, plot
  `conditional_posterior$mean[[m]][dims[1], ]` against
  `conditional_posterior$mean[[m]][dims[2], ]`; for one dimension, plot the
  selected feature against `(0:(K - 1))/(K - 1)`.
- [x] Include ordering labels and selected-feature assignment probabilities in
  panel titles. Never read an independently updated per-ordering fit as the
  source of the trajectory.
- [x] Keep `plot_type = "elbo"` as the fixed-M structural objective trace and
  route `plot_type = "mu"` separately before that branch.
- [x] Run:
  `Rscript -e 'devtools::test(filter="partition-m-cavi", stop_on_failure=TRUE)'`.
  Expected: zero failures and zero warnings.

### Task 2: Replace obsolete multi-ordering mathematical derivations

**Files:**
- Modify: `important_derivations/partition_coherent_elbo_derivation.Rmd`
- Regenerate: `important_derivations/partition_coherent_elbo_derivation.html`
- Modify: `important_derivations/known_noise_cavi_derivation.Rmd`
- Regenerate: `important_derivations/known_noise_cavi_derivation.html`
- Preserve as current single-ordering reference:
  `important_derivations/cavi_math_details.Rmd`
- Delete after reference:
  `important_derivations/MPCurve - 推导讲解-关于 replace potential trajectories.pdf`

- [x] Rewrite the partition derivation from the current joint model and the
  structured family
  `q(C) prod_j q(Z_j) q(U_j | Z_j)` without retaining an obsolete
  fully-factorized implementation section.
- [x] Derive the conditional Gaussian update exactly:
  `P_jm = lambda_jm Q_K + D_jm`, `P_jm mu_jm = b_jm`, with neither term
  multiplied by `w_jm`; cover estimated shared variance and vector/matrix known
  measurement SD.
- [x] Derive the ordering-specific `q(C_m)` update, explicitly showing where
  `w_jm` weights each feature's expected likelihood and uncertainty term.
- [x] Derive the shared variance update
  `sigma2_j = n^{-1} sum_m w_jm E_jm`, the ordering-specific lambda update,
  the complete conditional local block, the `q(Z_j)` softmax, and the total
  fixed-M objective.
- [x] Explain intrinsic-RW pseudo-determinants, null-space identifiability,
  why `w_jm = 0` no longer creates singularity, and why truly unidentified
  likelihood information can still leave an intrinsic posterior singular.
- [x] Explain annealing: objectives at different temperatures are not directly
  comparable; exact ELBO monotonicity is assessed within the T=1 segment.
- [x] Explain why the raw fixed-M structural ELBO is not a cross-M selection
  criterion and why greedy selection is temporarily unavailable.
- [x] Rewrite the known-noise derivation so its single-ordering and structural
  partition sections share notation and no longer refer to ordering-specific
  estimated variances, weighted conditional trajectory precision, or obsolete
  helper names/line numbers.
- [x] Render both Rmd files to HTML and search the Rmd/HTML outputs for obsolete
  claims including `lambda Q + w D`, `sigma_j^{2(m)}`, freeze/drop state, and
  statements that augmented mean-field is the current implementation.
- [x] Delete the internal Chinese PDF after its formulas have been independently
  checked against the implementation. Verify that it is absent from
  `git status` and cannot enter the source tarball.

### Task 3: Make S3 help pages implementation-complete

**Files:**
- Modify roxygen in: `R/07_mpcurve.R`
- Modify package object documentation in: `R/MPCurver-package.R`
- Regenerate: `man/fit_mpcurve.Rd`, `man/do_mpcurve.Rd`,
  `man/plot.mpcurve.Rd`, `man/print.mpcurve.Rd`,
  `man/summary.mpcurve.Rd`, `man/fitted_prior.Rd`, `man/mpcurve.Rd`

- [x] Document the normalized structural `mpcurve` schema: named position-prior
  and gamma lists, shared `params$sigma2`, `conditional_posterior`,
  `lambda_mat`, `partition$pi_weights`, fixed M, and derived `$fits`.
- [x] Correct `summary.mpcurve` documentation so partition summaries are
  described as canonical structural summaries rather than delegated
  per-ordering summaries. List the actual returned fields.
- [x] Correct `plot.mpcurve` documentation: structural scatter panels,
  structural objective for `"elbo"`, and conditional posterior mean paths for
  `"mu"`; document one- versus two-dimensional `dims` behavior.
- [x] Expand `print.mpcurve` documentation to state the fixed-M structural
  reporting semantics and absence of active/frozen/drop state.
- [x] Ensure `fitted_prior()` documents named position-prior lists and the
  fixed-M partition prior with optional ordering extraction.
- [x] Run `Rscript -e 'devtools::document()'`, then run `git diff --check` and
  inspect every regenerated Rd file for agreement with runtime behavior.

### Task 4: Synchronize README and public vignettes

**Files:**
- Modify: `README.Rmd`
- Regenerate: `README.md`
- Modify as needed: `vignettes/partition.Rmd`
- Modify as needed: `vignettes/fitness.rmd`
- Verify single-ordering scope: `vignettes/mpcurve_intro.Rmd`
- Leave archival-only: `vignettes/Intro.Rmd`

- [x] Add a concise README section stating that `intrinsic_dim = M >= 2` uses
  fixed-M structural VI, shared estimated feature variances, no freeze/drop
  state, and no current greedy cross-M selection.
- [x] Add and render a partition-vignette example of
  `plot(fit, plot_type = "mu")`, explaining that panels show
  `E[U_j | Z_j = m]` from the canonical conditional posterior.
- [x] Cross-check all partition and fitness examples against the current public
  argument defaults and object schema.
- [x] Regenerate README output with the repository's established rendering
  command and verify that generated figures remain intentional.

### Task 5: Rebuild and audit the public website

**Files:**
- Regenerate curated files under: `docs/`

- [x] Run `Rscript scripts/build_public_site.R`; do not call
  `pkgdown::build_site()` directly.
- [x] Confirm that the public site includes only the home page, curated
  reference pages, and `mpcurve_intro`, `partition`, and `fitness` articles.
- [x] Search committed `docs/` for obsolete user guidance: active/frozen/drop
  behavior, supported greedy forward/backward selection, weighted conditional
  precision, ordering-specific estimated sigma, or augmented mean-field as the
  current method.
- [x] Verify that `docs/` contains no references to `AGENTS`, `CLAUDE`,
  `important_derivations`, `internal/`, `log/`, or `plan/`.

### Task 6: End-to-end verification and update log

**Files:**
- Create: `log/2026-07-21-codex-structural-vi-docs-and-plot.md`
- Update checkboxes in this plan

- [x] Run focused plot and documentation tests.
- [x] Run `Rscript -e 'devtools::test(stop_on_failure=FALSE)'`; expected: zero
  failures and zero warnings, with only intentional retired-greedy skips.
- [x] Run `R CMD build .` and `R CMD check --no-manual MPCurver_0.3.0.tar.gz`;
  expected final status: `OK`.
- [x] Inspect the rebuilt partition article and structural-mu figures visually.
- [x] Run repository-wide searches over `R/`, `man/`, `README*`, `vignettes/`,
  `important_derivations/`, and `docs/` and classify any remaining historical
  terms as intentional legacy implementation references or remove them.
- [x] Record key changes, files, exact verification results, and remaining
  intentional legacy caveats in the update log.
- [x] Remove generated check artifacts and both explicitly authorized copies of
  the temporary Chinese reference PDF; preserve every other user file under
  `experiments/` unchanged.
