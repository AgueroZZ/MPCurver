# Simplify the fit_mpcurve interface

The user requested a component-by-component review of MPCurver's public
interface, beginning with `fit_mpcurve()`. The first decision removes the
algorithm selector because CAVI is the sole supported algorithm, together
with the freeze/drop controls that structural VI already ignored.

## Historical first stage: six argument removals

Removed these six arguments from the fitting signature, roxygen documentation,
generated Rd, and generated pkgdown reference:

- `algorithm`
- `freeze_unused_ordering`
- `freeze_unused_ordering_threshold`
- `freeze_feature`
- `freeze_feature_weight_threshold`
- `drop_unused_ordering`

The signature now has 44 named arguments, including `X`, plus `...`, down from
50 named arguments. A check rejects the removed names in `...`, including
explicit `NULL` values, before they can reach either inference backend. The
error instructs callers to omit them. Result objects retain
`algorithm = "cavi"` as provenance. The development NEWS entry documents the
required call changes.

This change concerns `fit_mpcurve()`. The continuation interface and lower-level
legacy implementations retain their existing compatibility controls pending
their own review. No other option was selected for removal in this step.

## Verification of the historical first stage

- Compared the complete fit objects against the wrapper from commit
  `dbb10542594ff6e8a93726ec3041040211243675`, and compared both objects after
  two additional continuation iterations. Objects agree at zero numeric
  tolerance for `intrinsic_dim = 1`, `2`, and `"auto"`. The common data use
  `simulate_dual_trajectory(n = 60, d1 = 4, d2 = 4, d_noise = 0,
  sigma = 0.15, seed = 102)`, fitting seed 42, PCA, K = 6, five single-ordering
  sweeps, two annealing steps, one inner sweep, five final partition sweeps,
  and zero stopping tolerances.
- Source-loaded full test suite: 681 passing expectations, zero failures,
  zero warnings, and four existing retired-semantics skips across 86 test
  cases. Updated the existing interface assertions and added rejection checks
  for removed arguments in single-ordering, partition, and automatic-M modes.
- Regenerated `man/fit_mpcurve.Rd` from roxygen sources and rebuilt only the
  fitting reference with pkgdown 2.2.0, matching the existing site builder.
  Roxygen2 7.3.2 emitted an older-version advisory relative to the retained
  `RoxygenNote: 7.3.3`; the generated Rd change contains only the requested
  argument removals.
- Inspected the actual HTML: every documented argument matches the source,
  the removed controls are absent, the KaTeX resources are retained, and all
  96 relevant local links, assets, and incoming/outgoing anchors resolve.
  The HTML diff contains only the requested removals. Browser screenshot
  inspection remains unverified because headless Firefox timed out.
- Source package build completed; it emitted two environment warnings that
  large numeric user/group IDs were normalized in the archive. No full
  `R CMD check` was performed for this interface-only change.
- `git diff --check` passed. No commit, push, or publication was performed.

## Historical follow-up: audit all 44 remaining arguments

At the user's request, traced all remaining arguments through the fitting
wrapper, initialization, prior updates, and stopping loops. The complete
inventory is `plan/2026-10-01-fit-mpcurve-argument-audit.csv`; its 44 argument
names were checked against the current signature with no missing or duplicate
entries. The recommendations below are proposals, not additional API changes.

- Rename the requested grid size to `num_bins`, with internal mathematical
  notation K retained. A scalar `initial_method` would identify the normal
  initialization choice; per-ordering overrides and multiple candidate starts
  should have separate documented interfaces. Currently `method` vectors have
  different meanings for one versus several orderings, and `num_cores > 1`
  determines whether the single-ordering wrapper fits all candidates rather
  than just changing execution speed.
- `tol` and `tol_outer` govern the same default quantity, absolute ELBO change
  divided by N times D. They differ by inference path: single-ordering versus
  joint partition updates after annealing at T = 1. Propose a single `tol`.
- `iter` caps single-ordering updates and supplies the default final partition
  budget when `max_converge_iter` is NULL. Propose one `max_iter` for that final
  fitting phase. Partition annealing performs `n_outer * inner_iter` additional
  joint sweeps, followed by at most `max_converge_iter` final sweeps. The default
  joint budget is 25 + 100, excluding initialization. Rename annealing controls,
  and consider fixing `inner_iter = 1` and `T_end = 1` in the ordinary interface.
- `greedy` has no supported nontrivial value and can be removed. Its eventual
  removal must also update `select_mpcurve_dimension()`'s internal candidate
  call, which still explicitly passes `greedy = "none"`.
- `assignment_prior` and `ordering_alpha` are deprecated but are not no-ops in
  new fits: the legacy Dirichlet branch uses posterior concentration updates,
  expected log weights, and a KL term, unlike adaptive EB. Retiring that branch
  is a model-support decision. The deprecated uniform alias maps to the current
  fixed partition prior. Continuation's no-op behavior must not be generalized
  to the fitting entry point.
- Propose an ordinary interface with 12 named arguments including `X`: `X`,
  `S`, `num_bins`, `num_orderings`, `initial_method`, `position_prior`,
  `partition_prior`, `max_iter`, `tol`, `verbose`, `init_control`, and `control`.
  Preserve important scientific alternatives in explicit advanced settings,
  combine precision/variance bounds, move effective-count thresholds to
  reporting, and retire obsolete controls. Current fitting defaults and
  supported models have not been changed by this audit.

Two targeted, read-only checks corroborated the source inspection. First,
correlation similarity retains all input columns but zeros off-diagonal
similarities for columns whose raw sample standard deviation is less than
`similarity_min_feature_sd`; diagonal similarities remain one. This protects
initialization and does not remove features from the likelihood. The threshold
is expressed in the input's units.

Second, spline similarity has another absolute guard:
`total_sum_squares > sqrt(.Machine$double.eps)`. For the 20-point sequence
`t = seq(0, 1, length.out = 20)`, the similarity between `t` and `2*t` is one,
but becomes zero when both columns are multiplied by `1e-5`, even though neither
column triggers the `1e-8` feature-SD threshold. This demonstrates a scale
sensitivity in the current numerical guards. Recorded for a focused correction
and validation before replacing them with an internal rule. No similarity or
fitting code was changed during this follow-up.

## Approved redesign implemented

Following the user's decisions, replaced the ordinary fitting signature with
12 named arguments and no ellipsis. `R/12_interface.R` defines that interface
and the exported, validated `mpcurve_control()` and `mpcurve_init_control()`
constructors. Each accepts a named settings list through the fitting wrapper;
unknown or duplicate options fail before initialization. The audit CSV retains
the historical proposals and adds the implemented disposition of all 44
arguments.

- Renamed `K`, `intrinsic_dim`, and `method` to `num_bins`, `num_orderings`,
  and `initial_method`. Exactly one initialization method is accepted for
  single, specified multiple, and automatic ordering counts. Removed
  `num_cores` and the candidate-list branch; comparisons use separate calls.
- Model controls include `rw_order`, `lambda_init`, `fix_lambda`,
  `sigma2_init`, `lambda_bounds`, and `sigma2_bounds`. For RW3 with smoothness
  precision held at 5, use
  `control = mpcurve_control(rw_order = 3, lambda_init = 5, fix_lambda = TRUE)`.
  Omitting `fix_lambda = TRUE` retains empirical-Bayes precision updates.
- Unified `tol` and the final-phase cap `max_iter`. Advanced partition
  annealing uses `anneal_steps`, `anneal_start`, and `anneal_sweeps`, with an
  endpoint of 1. The default joint budget remains 25 annealing sweeps plus
  at most 100 final sweeps, excluding initial fits.
- Initialization controls retain the existing similarity metrics and defaults,
  automatic-cut eligibility rules, discretization, optional starting fits,
  and PCA component settings. `method_args` reaches the selected ordering
  helper through both fitting engines and automatic one-ordering fallback,
  independently of similarly named fitting settings. Invalid helper argument
  names are rejected before a subset initializer can fall back to PCA.
- Removed `greedy` and rewired `select_mpcurve_dimension()` to the new
  arguments/constructors. Its comparison still requires soft assignments
  and a fixed uniform partition prior.
- Removed the exploratory Dirichlet posterior/ELBO branch. The ordinary
  entry points reject `assignment_prior` and `ordering_alpha`. Historical
  internal callers retain explicit retirement checks; old result fields and
  historical method code are not treated as supported new-fit options.
- `do_mpcurve()` now uses `tol` for both model types and removes the seven
  ignored freeze/drop and assignment-prior arguments. Other continuation
  controls are retained for the separate review of that component.
- Scientific defaults, numerical low-variance rules, and returned mathematical
  field names remain unchanged. The scale-sensitivity finding above remains
  open. Version 0.3.4 is the base version; this is a local development API
  change rather than a published release.

Current verification of the implemented redesign:

- 696 passing test expectations, zero failures, errors, warnings, or skips.
  Dedicated checks cover RW3, fixed versus estimated lambda, precision/noise
  bounds, shared stopping tolerance, validation, helper argument forwarding,
  and continued fitting.
- Fit and continuation results agree at zero numerical tolerance with the
  source at `dbb10542594ff6e8a93726ec3041040211243675` for one ordering, two
  orderings, and automatic ordering initialization, using matched settings.
  Compared responsibilities, model parameters, ELBO/objective traces, and
  temperature traces.
- Roxygen/Rd, README, and all vignette call sites are migrated. Added reference
  topics for both control constructors and numbered figure/table captions in
  edited articles. Full source-build, installed-package checks, and generated
  HTML verification are being completed; final outcomes are recorded below.

No commit, push, or publication is authorized or performed.

Final validation outcomes:

- Source package build and full `R CMD check --no-manual` completed with
  `Status: OK` (zero errors, warnings, or notes). Installed-package tests
  passed 678 expectations, with one expected CRAN snapshot skip. All three
  public vignette code runs and output rebuilds passed. The checked archive
  and final checkout have identical parsed runtime R code in all 19 R files;
  subsequent changes concern documentation comments and site generation.
- Rebuilt README and the curated pkgdown site with pkgdown 2.2.0, then refreshed
  the final reference/article/home outputs from their canonical sources.
  HTML verification covers all 29 pages and 973 local links, images, scripts,
  and anchors, with no missing targets. The fitting reference documents exactly
  the 12 public arguments, and both control constructors link back to fitting.
- Audited the edited public pages: the homepage has three numbered figures;
  the getting-started article has four, partitioning has three figures and
  one table, and fitness has nine figures and four tables. Numbered captions
  and text references are present. Inspected representative actual plot images
  for readable labels, panel titles, and clipping. Browser screenshot layout
  verification was not repeated after the earlier Firefox timeout.
- README's temporary HTML preview initially warned about fetching MathJax in
  the network-restricted environment; the final Markdown rendering disables
  that unneeded preview and completed without that warning. Public pages retain
  their KaTeX configuration. Package repository-index connectivity messages did
  not produce check warnings or notes.
- `git diff --check` passes. Checked archives, temporary comparisons, and
  diagnostics remain under `/tmp`; no commit, push, or publication occurred.

## Follow-up: explain advanced controls and isolate spline scale sensitivity

Expanded the canonical roxygen documentation for all 16 `mpcurve_control()`
settings, with explicit defaults, dimensions, units, applicable fitting paths,
and effects on inference versus reporting. Grouped the explanation into
smoothness/noise, position/partition probabilities, and optimization/reporting.
Added examples for precision bounds and annealing budgets. Regenerated Rd and
pkgdown reference/example outputs and refreshed search. Constructor examples and
Rd validation passed; 29 generated pages and 979 local links/assets/anchors
resolve. The parsed runtime code is unchanged from the checked archive, so the
prior full package check remains applicable to runtime behavior.

Reconfirmed that `num_cores` is absent from the fitting signature and that a
method vector fails before fitting for one, multiple, and automatic ordering
counts. Resource-dependent candidate fitting is removed.

Reproduced the unresolved spline issue using 20 points and features `t`, `2*t`:

- At original scale, feature SDs are approximately 0.311373 and 0.622745,
  total sums of squares are 1.842105 and 7.368421, and similarity is 1.
- At scale `1e-5`, SDs are approximately `3.113726e-6` and `6.227452e-6`,
  total sums of squares are `1.842105e-10` and `7.368421e-10`, and similarity
  is 0. Neither column triggers the separate `1e-8` feature-SD guard.
- In `R/10_partition_cavi.R`, valid responses must also have total sum of
  squares greater than `sqrt(.Machine$double.eps)`, approximately `1.490116e-8`.
  That absolute threshold uses squared data units. Failed responses retain the
  initial directional score of zero, and the symmetric off-diagonal similarity
  becomes zero when both responses fail.
- Since SSE and SST both scale quadratically, their ratio and hence R-squared
  should retain their values under this rescaling. At 20 samples, the absolute
  SST threshold corresponds to an SD cutoff of approximately `2.800485e-5`,
  substantially above the separate SD guard. It can therefore change feature
  grouping and automatic ordering-count initialization when units change.

This follow-up diagnoses the cause; it does not change spline inference or
similarity guards. A focused fix should evaluate dimensionless similarity using
centered/scaled nonconstant responses and a scale-aware constant-column rule,
without adding an ordinary fitting option. No commit or push was performed.

## Spline follow-up resolved

Following the user's variance-one normalization proposal, implemented internal
centering and unit sample variance before spline scoring, moved its SD guard to
standardized units, removed the raw-SST threshold, and computed SSE directly
from residuals. The scale-sensitivity issue recorded above is now resolved for
representable nonconstant inputs. Full regression suite: 731 passing checks;
ordinary-scale score differences are at roundoff level. Details and current
verification: `log/2026-10-01-spline-r2-normalization.md`.
