# Remaining refactor review decisions and fixes

The user authorized fixes for R5 through R10, removal of final hard assignment
for R11, and closure of R4 without compatibility work for older fitted objects.
The review report retains the original reproductions and records these final
dispositions.

## Changes

- **R4:** old-version fit continuation is outside the supported scope. No
  migration or compatibility layer was added. Frozen InferOrder analyses and
  their versioned libraries are unchanged.
- **R5:** sinusoidal generators sample an index into the admissible frequency
  vector. A scalar frequency now specifies that exact frequency; vector inputs
  retain their previous sampling behavior.
- **R6:** GP draws are explicitly shaped as feature-by-sample matrices before
  transposition. A block containing one feature remains an `n x 1` matrix,
  including through row permutation and noise addition.
- **R7:** single- and multi-ordering scatterplots map posterior-mean pseudotime
  directly to the common `[0,1]` palette. The 2D embedding path receives these
  precomputed colors instead of rescaling the observed range. Trajectory
  segments too short for the graphics device's arrow renderer are skipped,
  avoiding undefined-angle warnings for nearly constant posterior curves.
- **R8:** plotting branches merge named user graphics arguments with defaults
  before calling the plotting function. Titles and axis labels work for 1D/2D
  scatterplots, trajectory plots, and single/partition objective traces. Single
  objective plots apply axis arguments to the plot rather than passing them
  into `lines()`; line styling is forwarded separately.
- **R9:** the shared default grid is
  `max(rw_order + 1, min(50, floor(nrow(X) / 5)))`. Too few samples for that
  automatic minimum, or an explicitly insufficient bin count, raises a clear
  error. Single, grouped, and automatic fitting use the same rule. Initialization
  cuts can still collapse on tied scores, as documented.
- **R10:** both spiral manual examples now load the complete development package
  with `pkgload::load_all(".")` and fit through `fit_mpcurve()`. They no longer
  source selected internal files. Their local simulation-scoring and plotting
  helpers remain in each example; no external analysis helper is required.
- **R11:** `hard_assign_final` and final weight-hardening branches are removed
  from the public control and implementation. Fits retain posterior soft
  assignments. `fit$partition$assign` remains a most-probable label for reporting
  and does not mutate those probabilities. No postprocessing point is appended
  to the CAVI iteration history. The active derivation now describes this behavior.

## Verification

- Full source suite: **1,182 passing expectations**, zero failures, errors,
  unexpected warnings, or skips. Regression coverage includes scalar-frequency
  recovery from noiseless generated trajectories, GP singleton dimensions,
  random-walk grid requirements, captured plot arguments/colors, and soft
  partition iteration accounting. Evidence: `/tmp/mpcurver-r5-r11-tests.R`,
  `.log`, and `.rds`.
- Independent minimum-sample checks pass for RW1 with two bins/samples, RW2
  with three, and RW3 with four. The narrow-pseudotime reproduction now uses
  midpoint palette colors in both 1D and 2D. The generated plot was inspected.
  Evidence: `/tmp/mpcurver-r5-r11-evidence.R`, `.log`, and
  `/tmp/mpcurver-r7-r8-plots.png`.
- Both manual scripts run to completion with their full default simulation
  sizes and fitting budgets, including the optional PNG output path. Their
  printed ELBO monotonicity diagnostics pass. Evidence:
  `/tmp/mpcurver-r10-spiral.log`, `/tmp/mpcurver-r10-largek.log`, and matching PNGs.
- Roxygen reference generation, README execution, and the complete curated
  pkgdown build succeed, including all four tutorials. The existing roxygen
  version advisory remains (installed 7.3.2 versus recorded 7.3.3).
- Browser checks pass for ten pages: home, controls, fitting, plotting,
  dimension selection, four tutorials, and news. All inspected images load;
  no math errors or page-width overflow occur at a 1,053-pixel viewport.
  The retired option is absent from the control reference, and the revised
  grid/color descriptions are present. Evidence:
  `/tmp/mpcurver-r5-r11-browser.py` and `.log`.
- Site audit passes **33 HTML pages and 1,131 local links/assets/anchors**,
  including figure/table numbering checks. Whitespace checks pass.

- Full source build and `R CMD check --no-manual`: **Status: OK**, with no
  check errors, warnings, or notes. Installed tests pass 1,164 expectations,
  with one standard CRAN snapshot skip. Examples, all four vignette code runs,
  and vignette output rebuilding pass. All 21 runtime R files in the checked
  archive match the working tree byte for byte. Evidence:
  `/tmp/mpcurver-r5-r11-check/build.log`, `check.log`, and
  `MPCurver.Rcheck/00check.log`. Archive creation emits the existing large
  numeric user/group ID normalization warnings; the package check is clean.

The first build attempt encountered an unrelated legacy `RevoUtils` startup
profile; the rerun disables site/user R startup profiles while retaining the
required library paths. No package code was changed for this environment issue.

No commit or push.
