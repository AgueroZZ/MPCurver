# Public component interface cleanup, 2026-10-01

Reviewed the remaining exported functions and registered inspection methods
following the `fit_mpcurve()` redesign. The ordinary entry points keep common
scientific settings; advanced settings have documented defaults and explicit
scopes. Removed arguments are rejected instead of silently ignored.

## Resulting interfaces

| Function | Previous arguments | Current arguments | Main decision |
| --- | ---: | ---: | --- |
| `do_mpcurve` | 12 | 5 | Additional `max_iter`, optional `tol`, `verbose`, continuation control; inherit saved data and errors. |
| `select_mpcurve_dimension` | 4, including `...` | 12 explicit | Rename cap to `max_num_orderings`; expose accepted fit settings; retain fixed uniform comparison model. |
| `fiedler_ordering` | 6 | 3 | `X`, `num_neighbors`, advanced graph control. |
| `isomap_ordering` | 12 | 5 | `X`, `num_neighbors`, `component`, `seed`, advanced control; internal stabilization constant. |
| `pcurve_ordering` | 8 | 5 | `X`, `smoother`, `max_iter`, `tol`, advanced control. |
| `tSNE_ordering` | 8 | 6 | Keep component, perplexity, iteration cap, and seed; move embedding conventions into control. |
| `simulate_cavi_toy` | 9 | 5 | `n`, `d`, `num_bins`, `seed`, advanced generating-model control. |
| `simulate_intrinsic_trajectories` | 14 | 7 | Keep sample/block counts, `noise_sd`, family, seed; trajectory shapes and positions in control. |
| `simulate_dual_trajectory` | 15 | 9 | Keep block counts, `noise_sd`, family, crossing, seed; trajectory shapes in control. |
| `simulate_two_order_gp_dataset` | 11 | 5 | `n`, `d`, `noise_sd`, `seed`, advanced GP/permutation control. |
| `PCA_ordering` | 5 | 5 | Retain coherent coordinate/preprocessing settings; reject fractional components. |
| `simulate_spiral2d` | 6 | 6 | Retain geometry; rename observation noise to `noise_sd`. |
| `simulate_swiss_roll_1d_2d` | 4 | 4 | Retain geometry; rename observation noise to `noise_sd`. |
| `fitted_prior` | 4 | 4 | Retain standard S3 signature; reject unused extra arguments. |
| `plot.mpcurve` / `print.mpcurve` / `summary.mpcurve` | unchanged | unchanged | Existing inspection interfaces are concise; graphical `...` remains useful and S3 signatures are retained. |

The selector's former `...` concealed its real fitting surface. Its new
signature explicitly lists accepted settings and omits prior modes controlled
by the comparison procedure. One ordering method is still selected by
`initial_method`; trying two partition initialization paths for each M > 1
is part of the selector's documented comparison procedure.

The complete argument-by-argument decisions are in
`plan/2026-10-01-public-component-argument-audit.csv`. The preceding 44-argument
fit audit remains in `plan/2026-10-01-fit-mpcurve-argument-audit.csv`.

## Advanced settings and migration

`mpcurve_continue_control()` has six settings, all defaulting to NULL/inherit:
`lambda_init`, `fix_lambda`, `lambda_bounds`, `sigma2_bounds`,
`lambda_sd_prior_rate`, and `convergence`. A zero rate removes the penalty.
Do not pass the new-fit constructor `mpcurve_control()` to continuation:
its defaults specify a new model rather than inherited settings.

```r
fit <- do_mpcurve(fit, max_iter = 100)
fixed <- do_mpcurve(fit, max_iter = 100,
  control = mpcurve_continue_control(lambda_init = 5, fix_lambda = TRUE))
estimated <- do_mpcurve(fixed, max_iter = 100,
  control = mpcurve_continue_control(fix_lambda = FALSE))

selected <- select_mpcurve_dimension(X, max_num_orderings = 4,
  num_bins = 20, initial_method = "PCA", max_iter = 100)

fit <- fit_mpcurve(X, initial_method = "isomap",
  init_control = mpcurve_init_control(method_args = list(
    num_neighbors = 10, control = list(num_landmarks = 100))))

sim <- simulate_intrinsic_trajectories(n = 100, d_signal = c(5, 5),
  noise_sd = 0.1, trajectory_family = "linear", seed = 1,
  control = list(linear_slope_range = c(1, 2)))
```

Ordering and simulation helpers accept uniquely named control lists. Their
reference pages enumerate every supported setting, default, and interpretation.
They retain separate method/generator scopes rather than sharing irrelevant
options through a universal control object. Unknown settings are rejected.
No obsolete argument aliases were added.

## Continuation corrections

The old single-ordering `lambda` override implicitly fixed precision and
ignored supplied bounds; the partition override only reset its value. New
continuation consistently separates resetting precision from fixing it.
Inherited fixed precisions are clipped when bounds tighten, and estimated
noise variances are constrained. A revised partition objective is evaluated
under changed precision/prior/variance settings before testing the next
objective increment. Historical trace values are retained unchanged.

Removed an unreachable legacy augmented-partition continuation branch after
its unconditional rejection. Historical fits remain readable. New CAVI
continuations retain ordering-helper provenance and record the requested
continuation budgets/overrides in `$continuation_history`. Selection records
retain their candidates/history and mark continuation explicitly. New selector
records use `max_num_orderings`; inspection still accepts the historical
`max_intrinsic_dim` record field.

## Verification

Validation results are recorded below after the final checks. No commit or
push is authorized. Experiment libraries and completed research results are
unchanged.

- Local suite: **798 passing expectations**, zero failures, errors, warnings,
  or skips (`/tmp/mpcurver-component-final-tests.log`, results RDS).
- Independent before/after examples: all five ordering helpers, four larger
  simulation generators, single/partition continuation parameters and traces,
  and selector candidate/history tables match with zero numerical tolerance.
  Baseline seed 778, 40-by-4 input; helper seeds 3 and generator seeds 2.
  Evidence: `/tmp/mpcurver-compare-components-final.log` and the comparison
  script/baseline RDS in `/tmp`. Defaults retain the generating distributions.
- Full curated pkgdown rebuild completed with the pinned pkgdown 2.2.0 builder,
  executing all three public tutorials. Reviewed new function references and
  actual figures. Increased the tutorial comparison figure's title margin
  and rerendered it successfully after the package archive was built; this
  final edit affects layout only.
- **30 HTML pages and 996 local links/assets/anchors** checked without issues;
  figure/table numbering and caption/text references pass; `git diff --check`
  passes. Refreshed `fitted_prior` reference/search after its final wording edit.
- The check archive's parsed code matches all **20 current R source files**.
  Archive/check paths: `/tmp/mpcurver-component-check-20261001/`. Installed
  tests pass 780 expectations, with one standard CRAN snapshot skip. Final
  package-check status is recorded below after vignette rebuilding completes.

- Full source build and `R CMD check --no-manual`: **Status: OK**, with zero
  errors, warnings, or notes. Examples, installed tests, all three vignette
  executions, and vignette-output rebuilding passed. `00check.log` and the
  source archive remain under `/tmp/mpcurver-component-check-20261001/`.
