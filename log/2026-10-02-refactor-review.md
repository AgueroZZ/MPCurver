# MPCurver 0.4.0 refactor review

This review checks the local 0.4.0 software after the interface and implementation
changes made on October 1. It identifies reproducible defects and behavior that
needs a decision before repair. No software fixes, test changes, commits, or
publication are part of this review.

The review baseline is commit `dbb10542594ff6e8a93726ec3041040211243675`
(0.3.4); the reviewed implementation includes the uncommitted changes and new
`R/12_interface.R`, `R/13_helpers.R`, and `R/14_dimension_finalization.R`.
Findings distinguish new regressions from existing problems. Intended breaking
argument removals and the default removal of annealing are not defects.

Follow-up status: R1 is addressed in
[Numerical-zero ordering prior weights](2026-10-02-prior-weight-pruning.md).
Its original evidence below is retained; oversized settings now warn, and
zero retained probability raises an error instead of returning NaN. R2 has been
reassessed: its counterexample concerns an oversized pruning threshold, and no
convergence change is proposed for numerical-zero cleanup. The user has agreed
to close R1 as fixed and R2 as requiring no runtime change. R3 is now fixed in
[Explicit initialization failure handling](2026-10-02-initialization-failure-policy.md).
R4 is closed by the decision not to support old-version fits. R5 through R11
are addressed in [Remaining review fixes](2026-10-02-remaining-review-fixes.md).
Verification below describes the original review snapshot; follow-up logs
record checks of the subsequent fixes.

## Verification

- R 4.3.3; BLAS and OpenMP each limited to one thread.
- Current source suite: **1,036 passing expectations**, zero failures, errors,
  warnings, or skips. Evidence: `/tmp/mpcurver-review-tests.log` and
  `/tmp/mpcurver-review-tests.rds`.
- Independent continuation comparisons: six uninterrupted sweeps agree with
  three sweeps followed by three continuation sweeps for all **24 combinations**
  of one/two orderings, RW1/RW2/RW3, estimated/known observation errors, and
  adaptive/fixed position priors. Parameters, responsibilities, objectives, and
  precision traces agree within `1e-10`. Data seed 849, initialization seed 850,
  40 samples, six features, six bins, and `tol = 0`. Evidence:
  `/tmp/mpcurver-review-continuation.R` and `.log`.
- Current site audit: **33 HTML pages and 1,130 local links, assets, and anchors**
  pass, including the existing figure/table caption checks. Canonical README,
  the four curated tutorials, and public reference signatures are consistent
  with the current API. This was a read-only audit, not a new site render.
- Public summary and plot dispatch also works with `pkgload::load_all()` using
  `export_all = FALSE`, `helpers = FALSE`, and `attach_testthat = FALSE`.
- Full source build and `R CMD check --no-manual`: **Status: OK**, zero check
  errors, warnings, or notes. Examples, installed tests, all four vignette code
  runs, and vignette rebuilding pass. Installed tests report 1,018 expectations
  and one standard CRAN snapshot skip. All 21 runtime R source files in the
  checked archive match the working tree byte for byte. Evidence:
  `/tmp/mpcurver-review-check/MPCurver.Rcheck/00check.log` and
  `/tmp/mpcurver-review-check.log`. Archive creation emitted two environment
  warnings about normalizing large numeric user/group IDs; the package check
  itself is clean.
- `git diff --check` passes for the reviewed working tree.

The passing suite does not cover the counterexamples below. This review did not
rerun the historical scientific simulation studies or prove correctness for
every supported numerical input.

## Issue list

P2 denotes a functional or scientific-output issue under a supported input or
setting. P3 denotes a narrower usability or maintenance issue. None of the
confirmed findings establishes a failure of ordinary default fitting on the
tested data. IDs are stable so that repair decisions can be discussed separately.

| ID | Priority | Finding | Origin |
| --- | --- | --- | --- |
| R1 | Closed (fixed) | Automatic pruning can return a non-finite objective; addressed in follow-up | New finalization code |
| R2 | Closed (no change) | No default-threshold defect established; oversized-threshold example retained for context | New finalization code |
| R3 | Closed (fixed) | Invalid helper settings now error; numerical fallback warns and can be disabled | New forwarding exposed existing fallback; PCA validation gap repaired |
| R4 | Closed (out of scope) | Continuation compatibility for old-version fits is not required | New retirement checks interact with legacy metadata |
| R5 | Closed (fixed) | Scalar sinusoid frequencies retain their supplied value | Existing generator bug |
| R6 | Closed (fixed) | Singleton GP blocks retain sample-by-feature dimensions | Existing generator bug |
| R7 | Closed (fixed) | Pseudotime colors use the absolute [0,1] scale | Existing plotting issue |
| R8 | Closed (fixed) | Named plotting arguments override defaults consistently | Existing plotting issue |
| R9 | Closed (fixed) | Default bin counts respect the random-walk order | Existing default-setting issue |
| R10 | Closed (fixed) | Manual examples load the package and use the public fitting API | Existing maintenance debt compounded by file splitting |
| R11 | Closed (removed) | Final hard assignment removed; posterior weights remain soft | Existing trace-accounting issue |

## R1 Automatic pruning can return a non-finite objective

**Location:** `R/14_dimension_finalization.R:63-64`.

After selecting retained orderings, the code divides each feature's remaining
assignment probabilities by their row sum. If a removed ordering contained all
of one feature's representable probability, that row sum is zero and division
produces `NaN`. Hard assignments cause exact zeros; sufficiently concentrated
soft assignments can also underflow to zero.

The threshold below is accepted because `4.5 < ncol(X) / 2 = 5`. This is a
nondefault pruning threshold; it deliberately removes a supported feature group.

```r
sim <- simulate_dual_trajectory(
  n = 40, d1 = 6, d2 = 4, d_noise = 0, noise_sd = 0.01,
  trajectory_family = c("quadratic", "monotone"), seed = 1
)
fit <- fit_mpcurve(
  sim$X, intrinsic_dim = "auto", num_bins = 10, max_iter = 20, tol = 0,
  init_control = mpcurve_init_control(max_intrinsic_dim = 2),
  control = mpcurve_control(effective_count_tol = 4.5)
)
fit$intrinsic_dim
tail(fit$elbo_trace, 1)
fit$dimension_estimation$expected_feature_counts
# 1; NaN; NaN
```

A second reproduction uses noise SD `0.05`, eight bins, ten sweeps, and
`hard_assign_final = TRUE`, with the other settings unchanged. It also returns
a one-ordering fit with a `NaN` objective and expected count rather than an error.
R1 and R2 were independently rerun from
`/tmp/mpcurve_dimension_review/summary.R`; output is retained in
`/tmp/mpcurver-review-pruning.log`.

**Repair decision:** require positive retained mass for every feature before
normalization. Decide whether unsupported pruning should fail clearly or use
an explicitly defined reassignment policy. Merely replacing `NaN` with arbitrary
probabilities would change the statistical operation. There is no evidence here
that the default `1e-8` cutoff produces this failure on the tested data.

## R2 Oversized pruning thresholds can invalidate inherited convergence

**Location:** `R/14_dimension_finalization.R:107` and `:116` after the R1 changes.

**Reassessment (2026-10-02 UTC):** the original counterexample does not establish
a defect in numerical-zero cleanup. It sets `effective_count_tol = 4.5` for ten
features, allowing an estimated prior weight as large as 0.45 to be removed.
That operation can discard a supported ordering. The original recommendation
to reset convergence after every pruning event was too broad.

For the adaptive prior, `omega_m = sum_j w_jm / P`. Thus an exactly zero
estimated prior implies zero feature-assignment weight for that ordering. Its
weighted local blocks and assignment terms vanish. The remaining position
contribution is

\[
\mathcal C_m = -\sum_i \operatorname{KL}(q(C_i^{(m)})\|\pi^{(m)}).
\]

At a coordinate-stationary solution, zero assignment weights make the position
update equal its prior, so this term vanishes as well. Removing that ordering
therefore leaves the mathematical objective and retained responsibilities
unchanged. For small positive prior weights, this is an approximation whose
error also depends on the local scores and position distributions; the prior
threshold alone is not a universal ELBO error bound for arbitrary states.

A follow-up check uses the public automatic fit and the default `1e-8 / P`
prior-weight cutoff, without editing any fitted responsibilities:

```r
set.seed(2)
X <- matrix(rnorm(240), 30, 8)
fit <- fit_mpcurve(
  X, intrinsic_dim = "auto", num_bins = 5, max_iter = 200, tol = 1e-6,
  init_control = mpcurve_init_control(max_intrinsic_dim = 3)
)
event <- fit$dimension_estimation$pruning[[1]]
event[c("prior_weights", "objective_before", "objective_after")]
```

The fit converges and reduces three orderings to one. Removed prior weights
are approximately `3.8911e-27` and `1.9316e-27`. Both recorded objectives are
`-337.4580809164822`; their difference is exactly zero in double precision.
Seeds 1 through 4 were checked with these settings; only seed 2 prunes, while
the other three initialize with one ordering. Evidence:
`/tmp/mpcurver-r2-numerical-zero.R`, `.log`, and `.rds`. This is a targeted
check, not a universal guarantee or a new full-suite/package check.

**Original oversized-threshold counterexample:** the finalizer copies
`fit$converged` into the reduced model. When pruning materially changes the
model, that flag need not describe the reduced state:

```r
sim <- simulate_dual_trajectory(
  n = 40, d1 = 6, d2 = 4, d_noise = 0, noise_sd = 0.05,
  trajectory_family = c("quadratic", "monotone"), seed = 1
)
fit <- fit_mpcurve(
  sim$X, intrinsic_dim = "auto", num_bins = 8, max_iter = 200, tol = 1e-4,
  init_control = mpcurve_init_control(max_intrinsic_dim = 2),
  control = mpcurve_control(effective_count_tol = 4.5)
)
more <- do_mpcurve(fit, max_iter = 1, tol = 1e-4)
c(returned_converged = fit$converged,
  before = tail(fit$elbo_trace, 1), after = tail(more$elbo_trace, 1))
abs(tail(more$elbo_trace, 1) - tail(fit$elbo_trace, 1)) / length(sim$X)
```

The two-ordering optimization converges in 25 sweeps, then pruning returns one
ordering with `converged = TRUE`. Its objective is approximately `-1128.258`;
one continuation sweep raises it to `136.255`. The normalized change is
`3.161282`, far above `tol = 1e-4`.

**Current recommendation:** retain the existing convergence behavior for
numerical-zero cleanup. Do not add fitting sweeps or automatically reset the
flag based on this example. Keep R2 as a conditional concern for aggressive
nondefault pruning, covered in part by R1's oversized-threshold warning, rather
than a confirmed defect at the default cutoff. No runtime code was changed in
this reassessment.

**Disposition:** closed with the user's agreement; no algorithm change is
needed for the intended numerical-zero cleanup behavior.

## R3 Invalid helper settings silently invoke another initializer

**Disposition:** fixed following the user's decision. Invalid settings error
before fallback, ordering computation failures warn before trying PCA, and
`mpcurve_init_control(on_failure = "error")` stops instead. PCA failure stops;
requested and actual methods and the failure cause remain in initialization
metadata. The original evidence below is retained. See the linked follow-up
log for the implementation and current checks.

**Original review locations:** `R/12_interface.R:420-427`,
`R/10_partition_cavi.R:1107-1134`, and `:1018`.

The wrapper validates helper argument names and nested `control` values, but
some top-level helper values are validated only inside grouped fitting. A broad
error handler treats their validation failures as numerical initialization
failures and falls back to PCA. The resulting fit has no warning under ordinary
verbosity.

```r
set.seed(614)
X <- matrix(rnorm(200), 40, 5)
fit <- fit_mpcurve(
  X, num_bins = 5, intrinsic_dim = 2, max_iter = 0,
  initial_method = "fiedler",
  init_control = list(method_args = list(num_neighbors = 1))
)
fit$fit$init_info$A[c("method_requested", "method_used", "fallback_reason")]
# requested: fiedler; used: PCA
# reason: num_neighbors must be a single integer >= 2.
```

The equivalent single-ordering call correctly errors. Consequently, comparing
fits labelled by the requested initializer can compare different methods from
those specified. The fallback predates the refactor, while the new `method_args`
path exposes the inconsistent validation behavior.

There is a related, older validation gap: `init_control =
list(pca_components = 1.9)` is rejected for one ordering but silently converted
to component 1 for two orderings. The grouped path coerces with `as.integer()`
before checking; `mpcurve_init_control()` itself does not validate this field.

**Repair decision:** validate deterministic option errors before the fallback
handler and use the same integer checks across dimensions. Decide separately
whether genuine numerical fallback should be reported or remain available.

## R4 Supported fixed-uniform saved fits fail during continuation

**Disposition:** closed by user decision. Compatibility with old-version
fits is outside the current software scope; no migration layer was added.

**Locations:** `R/07_mpcurve.R:210-223` and
`R/10_partition_cavi.R:1457-1458`.

A structural fit made with 0.3.4's deprecated `assignment_prior = "uniform"`
stores that name in its metadata. The current continuation wrapper accepts the
structural object and runs updates, but result conversion forwards the retained
metadata into a validator that now rejects the retired argument. No retired
argument was supplied by the continuation caller.

Reproduction uses an actual fit created by loading the exact HEAD archive:

```r
# Under the baseline 0.3.4 checkout:
set.seed(842)
X <- matrix(rnorm(160), 40, 4)
old_fit <- fit_mpcurve(
  X, K = 5, intrinsic_dim = 2, iter = 2, n_outer = 1,
  inner_iter = 1, max_converge_iter = 2, tol = 0, tol_outer = 0,
  assignment_prior = "uniform"
)
saveRDS(old_fit, "old-uniform-fit.rds")

# Under the reviewed 0.4.0 checkout:
do_mpcurve(readRDS("old-uniform-fit.rds"), max_iter = 1, tol = 0)
# Error: assignment_prior and ordering_alpha have been removed ...
```

An otherwise identical old adaptive fit continues successfully. Evidence:
`/tmp/mpcurver-review-compat.R`, `.log`, and
`/tmp/mpcurver-review-legacy-fits.rds`.

**Repair decision:** normalize the old uniform alias to the supported fixed
uniform model before computation. Old Dirichlet fits also error, but retirement
of that statistical model is intentional and is not the defect claimed here;
their continuation policy can be documented separately.

## R5 A scalar sinusoid frequency generates other frequencies

**Disposition:** addressed in the remaining-review follow-up linked above.
The reproduction below describes the original reviewed version.

**Location:** `R/08_SoftSmoothEM.R:1372`.

The shared generator uses `sample(sinusoid_freq, 1L)`. R interprets a single
positive integer such as `4` as the sequence `1:4`, rather than sampling the
sole admissible value. This affects both intrinsic and dual trajectory helpers.

```r
t <- seq(0, 1, length.out = 80)
sim <- simulate_intrinsic_trajectories(
  n = 80, d_signal = 8, d_noise = 0, noise_sd = 0, seed = 1,
  control = list(sinusoid_freq = 4, signal_range = c(1, 1),
                 latent_positions = matrix(t, ncol = 1))
)
sse <- sapply(1:4, function(f) {
  colSums(qr.resid(qr(cbind(sin(f * pi * t), cos(f * pi * t))), sim$X)^2)
})
apply(sse, 1, which.min)
# 3 1 3 1 2 2 3 1, in returned feature-column order
```

The reconstructed frequencies have residual sums of squares below `1.5e-29`;
forcing frequency 4 produces errors around 21 to 40. The baseline documentation
defined an admissible-frequency vector, so scalar input should specify one
frequency. This can change a simulation design silently.

**Repair direction:** sample an index into the supplied frequency vector.

## R6 Two-feature GP simulation has the wrong matrix dimensions

**Disposition:** addressed in the remaining-review follow-up linked above.
The reproduction below describes the original reviewed version.

**Location:** `R/benchmarking_curves.R:166-172`.

```r
simulate_two_order_gp_dataset(n = 10, d = 2, seed = 1)
# Error: subscript out of bounds
simulate_two_order_gp_dataset(
  n = 10, d = 2, seed = 1,
  control = list(permute_rows_block2 = FALSE)
)
# Error: non-conformable arrays
```

Two features is accepted by the public validator. With one feature per block,
`MASS::mvrnorm(n = 1)` returns a vector, so `t()` produces a `1 x n` matrix
instead of `n x 1`. Row permutation or subsequent noise addition then fails.

**Repair direction:** preserve the sample-by-feature dimensions explicitly for
singleton draws. This bug exists in the baseline implementation.

## R7 Pseudotime colors disagree with the displayed scale

**Disposition:** addressed in the remaining-review follow-up linked above.
The reproduction below describes the original reviewed version.

**Locations:** `R/07_mpcurve.R:2351-2353`, `:2392-2393`,
`:2484-2486`, and `:2522`.

Colors are rescaled from each panel's observed minimum and maximum to the full
palette, but the legend is always labelled zero to one. The effect can be large
when the inferred positions are almost indistinguishable.

```r
set.seed(44)
X <- matrix(rnorm(120), 40, 3)
fit <- fit_mpcurve(
  X, S = rep(100, 3), num_bins = 5, max_iter = 20,
  position_prior = "fixed"
)
range(fit$locations$mean$pseudotime)
# 0.499999986332260 0.500000018843998
plot(fit, dims = 1)
```

These positions nevertheless span palette indices 1 to 255. The displayed
legend therefore misstates their values and may visually exaggerate weak
ordering. Separate ordering panels also use separate implicit scales.

**Repair decision:** map the absolute `[0,1]` pseudotime directly to color, or
clearly label the actual within-panel range and rescaling convention. This
finding concerns displayed posterior means, not a claim about rank recovery.

## R8 Ordinary plot labels through the ellipsis error or are ignored

**Disposition:** addressed in the remaining-review follow-up linked above.
The reproduction below describes the original reviewed version.

**Locations:** `R/07_mpcurve.R:2514-2516` and `:2437-2439`, with equivalent
argument conflicts in the partition scatter and trajectory branches.

Using `fit` from R7:

```r
plot(fit, dims = 1, xlab = "My x axis")
# formal argument "xlab" matched by multiple actual arguments
plot(fit, plot_type = "mu", dims = 1, main = "Custom title")
# formal argument "main" matched by multiple actual arguments
```

Plot help advertises additional base graphics arguments. The implementation
supplies the same arguments explicitly before forwarding `...`, which prevents
ordinary user overrides. For partition ELBO plots the ellipsis is omitted
instead, so the arguments are silently ignored.

**Repair direction:** merge user graphics arguments with defaults once, and
consistently forward them across plotting branches. Default plots still work.

## R9 Automatic bin count conflicts with default RW2 on small inputs

**Disposition:** addressed in the remaining-review follow-up linked above.
The reproduction below describes the original reviewed version.

**Location:** `R/09_cavi.R:475-485`.

```r
set.seed(1)
fit_mpcurve(matrix(rnorm(40), 10, 4), max_iter = 0)
# rw_q must be an integer in {1, ..., K-1}.
```

The automatic rule selects two bins for fewer than 15 samples, whereas default
RW2 requires at least three. RW3 has an analogous conflict below 20 samples.
The interface permits these sample counts, so a normal default call fails with
an error expressed in old internal argument names. This is an existing issue.

**Repair decision:** choose the automatic grid jointly with the random-walk
order, or state and validate the required sample/grid limits clearly before
initialization. A sensible grid also depends on the number of available samples.

## R10 Two standalone manual scripts omit required source files

**Disposition:** addressed in the remaining-review follow-up linked above.
The reproduction below describes the original reviewed version.

**Locations:** `tests/manual/cavi_spiral_example.R:1` and
`tests/manual/cavi_spiral_largek_example.R:1`.

These scripts source a selected subset of R files. The new public wrappers
depend on the split interface/helper files, so direct execution now fails on
missing `.mpcurve_integer_setting`; additional ordering wrappers are also
missing. Other updated manual scripts load the complete checkout.

```sh
Rscript tests/manual/cavi_spiral_example.R
```

This is maintenance debt, not evidence that a previously working formal
analysis was broken: the baseline script already failed later because it did
not source `.mpcurve_elbo_change`. The scripts are excluded from package builds,
so `R CMD check` does not exercise them.

**Repair direction:** load the complete package checkout consistently in these
manual examples. Keep this separate from immutable historical results.

## R11 Final hard assignment is counted as an extra fitting iteration

**Disposition:** addressed in the remaining-review follow-up linked above.
The reproduction below describes the original reviewed version.

**Locations:** `R/10_partition_structural_cavi.R:1040-1054` and `:790`.

Final hard assignment appends an objective record without executing a CAVI
sweep. The iteration count is nevertheless derived from trace length, so the
reported fitting budget is one too large.

```r
set.seed(1)
X <- matrix(rnorm(160), 40, 4)
fit <- fit_mpcurve(
  X, num_bins = 5, intrinsic_dim = 2, max_iter = 0,
  control = mpcurve_control(hard_assign_final = TRUE)
)
fit$iter
# 1, although no annealing or final-phase CAVI sweeps were requested
```

The same call with `hard_assign_final = FALSE` reports zero. This does not mean
the optimizer actually exceeds its sweep budget; the problem is the displayed
iteration count and associated provenance. The trace-based counting and final
hardening record also exist in the baseline implementation.

**Repair direction:** record postprocessing separately from the number of
fitting sweeps. Evidence: `/tmp/mpcurver-review-hard-count.R` and `.log`.

## Scope exclusions and next decisions

InferOrder's formal studies explicitly load frozen 0.3.2 or 0.3.4 libraries;
their retained old calls are part of reproducibility and should not be migrated
as a consequence of this audit. Old exploratory pages are not current main
studies. No active production-doc argument mismatch was confirmed.

All eleven review items now have recorded dispositions. R2 requires no algorithm
change under numerical-zero cleanup, R4 is outside the supported compatibility
scope, and the other items are fixed or removed as agreed. Historical reproduction
snippets are retained above; current verification is recorded in the follow-up logs.
