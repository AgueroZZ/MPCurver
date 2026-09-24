# Structural VI Partition Implementation Plan

> **For agentic workers:** Implement this plan task-by-task with test-first checkpoints. Keep all code, comments, documentation, and logs in English.

**Goal:** Replace the augmented fully-factorized multi-ordering CAVI path with structural VI for every fixed `M >= 2`, using a shared feature variance vector and no active/freeze/drop numerical machinery.

**Architecture:** New partition fits use `q(C) prod_j q(Z_j) q(U_j | Z_j)` with one canonical top-level state. Conditional trajectory precisions are `lambda[j, m] * Q + D[j, m]` and never contain `w[j, m]`; feature weights enter only ordering-responsibility updates, shared-parameter updates, and the structural ELBO. Ordering-specific `cavi` objects remain derived compatibility views, while continuation reads only canonical state.

**Tech Stack:** R, Matrix, matrixStats, testthat, roxygen2.

---

### Task 1: Add structural state and coordinate updates

**Files:**
- Create: `R/10_partition_structural_cavi.R`
- Test: `tests/testthat/test-partition-structural-cavi.R`

- [x] Add strict conditional-Gaussian updates for estimated and known noise using `A = lambda * Q + D`, with an identifiability error on Cholesky failure and no weight floor, jitter, freeze, or ridge fallback.
- [x] Add exact weighted `q(C)`, shared `sigma2`, ordering-specific `lambda`, complete local branch ELBO, `q(Z)`, assignment-prior, and total structural-objective updates.
- [x] Add a canonical state containing shared `params$sigma2`, `position_pi`, `gamma`, conditional `posterior`, `lambda_mat`, `pi_weights`, objective/temperature histories, and `variational_family = "structured"`.
- [x] Generate `$fits` only from canonical state as compatibility views; never read them during structural continuation.

### Task 2: Route fixed-M fitting and continuation

**Files:**
- Modify: `R/10_partition_cavi.R`
- Modify: `R/07_mpcurve.R`
- Test: `tests/testthat/test-partition-structural-cavi.R`

- [x] Route `soft_partition_cavi()` and `soft_two_trajectory_cavi()` to the structural backend for new fits.
- [x] Make `sigma2_init` a shared explicit input, reject it with known `S`, and make `lambda_init`/`fix_lambda` apply consistently to partition fits.
- [x] Route `do_mpcurve()` continuation through canonical structural state and reject continuation of legacy augmented partition objects while retaining read-only conversion.
- [x] Keep `freeze_*` and `drop_unused_ordering` formals for one transition release; warn when explicitly supplied and ignore them.
- [x] Reject `greedy != "none"` with a clear message until a cross-M selection criterion is designed.

### Task 3: Normalize the public object schema

**Files:**
- Modify: `R/07_mpcurve.R`
- Test: `tests/testthat/test-mpcurve-cavi.R`

- [x] Expose one shared `mpcurve$params$sigma2` vector, named ordering-specific `params$pi`/`params$mu`, `gamma`, `conditional_posterior`, and `lambda_mat`.
- [x] Report one fitted `intrinsic_dim = M`; keep `active_intrinsic_dim` and `displayed_intrinsic_dim` only as compatibility aliases equal to `M`.
- [x] Remove active/frozen/drop state from new raw and wrapped objects and ensure print, summary, plot, and fitted-prior access work with the canonical schema.

### Task 4: Replace tests and document the model transition

**Files:**
- Create: `tests/testthat/test-partition-structural-cavi.R`
- Modify: `tests/testthat/test-partition-m-cavi.R`
- Modify: `tests/testthat/test-mpcurve-cavi.R`
- Modify: `vignettes/partition.Rmd`
- Modify: `vignettes/fitness.rmd`
- Modify: `DESCRIPTION`
- Create: `log/2026-07-21-codex-structural-vi-partition.md`

- [x] Test analytical conditional posterior equality for `w = 1`, `1e-12`, and `0`; shared-sigma formula; local ELBO; assignment update; and full objective.
- [x] Test `M = 2, 3`, RW2, `ridge = 0`, estimated noise, vector/matrix known `S`, continuation, and `M = 1` conditional-update reduction.
- [x] Test transition warnings, canonical object fields, read-only `$fits`, legacy continuation rejection, and greedy rejection.
- [x] Rewrite the partition derivation as structural VI, remove freeze/drop guidance, explain intrinsic pseudo-ELBO semantics, and record that fixed `M` is currently required.
- [x] Bump the package to `0.3.0`, regenerate roxygen documentation, run focused tests, the full test suite, and `R CMD check`; record exact results and any pre-existing failure in the update log.
