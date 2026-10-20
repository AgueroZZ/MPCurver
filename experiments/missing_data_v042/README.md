# Missing-data verification for MPCurver 0.4.2

This study checks that observed-entry CAVI agrees with an independently
evaluated variational objective and that complete-data results retain the
0.4.1 numerical behavior. Input matrices have samples in rows and features in
columns. The development baseline is commit
`4f8b583d7549d4068d0e810e93f90653cbcb2ed9` (MPCurver 0.4.1).

For a minimal worked example, see [the rendered toy report](toy_example.html)
and its [R Markdown source](toy_example.Rmd). It fits 60 samples and four
features with 35 missing entries and displays the fitted trajectory and traces.

## Design

`verify.R` generates 40 samples and six features from linear, quadratic, and
sinusoidal functions plus independent Gaussian noise with SD 0.25. All fits
use six position bins and disable early stopping with `tol = 0`.

The complete-data comparison uses seed 4218 and eight sweeps, crossing one or
two orderings, estimated noise or known feature/entry SDs, fixed or adaptive
position/partition priors, and intrinsic or ridge-0.1 trajectory priors. It
compares position probabilities, trajectory means and full covariances,
assignment probabilities, variance and smoothness estimates, and objective
histories against a separately installed copy of the baseline source.

The missing-pattern study uses seeds 4218 and 4219 and 25 sweeps. Its 72 fits
cross three masks, one or two orderings, three noise modes, and two trajectory
priors. The masks are independent 35-percent missingness, feature-specific
missingness probabilities from 5 to 75 percent, and complementary rectangular
blocks. The first row is retained and empty rows receive their original first
feature value, so every public fit satisfies the input contract. Even seeds
use fixed position/partition priors and odd seeds use adaptive priors. Known
feature SDs equal 0.3; known entry SDs range from 0.2 to 0.5, with NA at missing
observations.

The reference objective loops over individual observed entries and calls
`dnorm`, adding the posterior-variance correction, independently computed
prior quadratic forms and determinants, and position/assignment entropies.
It does not call the implementation's likelihood or objective helpers. The
checks require absolute objective agreement below `1e-6`, objective increments
above `-1e-6`, finite objectives, and unchanged original data. Existing
quantile-to-equal-width initialization fallback warnings are recorded; any
other warning fails the study.

A separate boundary case retains only one value of the last feature. The
intrinsic fit must stop with a singular or ill-conditioned precision error.
The ridge-0.1 fit must remain finite, monotone, and agree with the
independent objective.

## Results on 2026-10-18

These results are for the source identified in `source-sha256.txt`. They were
regenerated after the per-feature loops became masked matrix products, and
again after three changes: precision initialization from data-carrying bins,
singular-precision errors in all fits, and pairwise-complete similarities.
On 2026-10-19 they were regenerated four more times: after the unused
SmoothEM and csmoothEM code was removed and random starts were routed through
the common ordering initialization, after unreachable CAVI-era partition code
was removed, after compatibility code for objects saved by earlier versions
was removed, and after constant compatibility fields were removed from fitted
objects. The CSV results were byte-identical each time.

- All 24 complete-data configurations were exactly identical to 0.4.1 for
  the compared outputs; maximum absolute difference was zero.
- All 72 missing-pattern configurations passed. Maximum absolute scalar
  objective discrepancy was `1.023e-12`; minimum objective increment was
  `1.284e-7`.
- Twelve documented initialization fallback warnings occurred; there were
  no numerical or other warnings in those 72 fits.
- The intrinsic single-observation boundary case stopped with the singular
  precision error: its one observed value occupies one position bin. The
  proper-prior case passed with objective discrepancy `5.684e-14`.
- The full package source test suite passed 1,885 assertions, with zero
  failures, errors, warnings, or skips. This includes independent posterior
  linear-system checks, variance derivatives, RW3 and scale-prior cases,
  continuation equivalence, input validation, initialization, dimension
  selection, extraction, and plotting. Scalar per-entry references also
  cover a common offset of `1e4` and non-unit feature weights.
- The source package built successfully; `R CMD check --no-manual` finished
  with zero errors, warnings, or notes, including all vignette executions and
  output rebuilds. Generated documentation passed the local-resource and
  caption audits.

The CSV files retain all configuration-level results. `session-info.txt`
records the regression environment. `source-sha256.txt` identifies the tested
source files before a release commit exists. Package-check and documentation
results are recorded in `VALIDATION.md`.

## Reproduction

Use R with the package Imports and Suggests, plus `pkgload`. Install the
baseline source into `tmp/baseline-library`; keep development dependencies in
a separate library. From the development package root, run:

```sh
Rscript experiments/missing_data_v042/verify.R baseline
Rscript experiments/missing_data_v042/verify.R current
Rscript experiments/missing_data_v042/verify.R compare
Rscript -e 'testthat::test_local(".")'
R CMD build .
R CMD check --no-manual MPCurver_0.4.2.tar.gz
Rscript scripts/build_public_site.R
```

Fit snapshots and build logs belong under ignored `tmp/`; only the compact
verification evidence is tracked. These checks establish agreement with the
specified observed-data CAVI objective. They do not establish propriety of the
exact intrinsic latent-variable posterior or validity under nonignorable
missingness; see the accompanying missing-data derivation.
