# MPCurver 0.3.2 simulation regression

This release check reuses earlier simulation data and initialization seeds to
verify the new default `abs(delta_ELBO) / (N * D) < 1e-6` stopping rule.

## Reproduce

From the MPCurver repository root, with MPCurver 0.3.2, `digest`, and `clue` installed:

```r
system('Rscript --vanilla experiments/v032_convergence_regression/run.R')
system('Rscript --vanilla experiments/v032_convergence_regression/summarize.R')
```

`inputs.rds` freezes the required inputs and reference summaries. The
`prepare.R` script documents their extraction from the sibling InferOrder
checkout; preparation is unnecessary when using the supplied bundle. Full
new fits are saved in the ignored `fits/` directory, and reruns reuse them.
Delete only that generated directory to request fresh fits.

## Comparisons

- **Website examples:** M = 1 Swiss roll (N = 1000, D = 2, K = 60) and M = 2
  feature partition (N = 500, D = 20, K = 50), at their original iteration
  budgets and Fiedler initialization. The original M = 2 fit reached its
  iteration cap; this check compares its fixed-budget estimates and does not
  require that example to converge within 100 post-annealing sweeps.
- **ND stopping:** 27 M = 1 datasets spanning N = 100, 200, 400 and D = 10,
  20, 40 with three seeds. Require the release default to reproduce the
  earlier experimental normalized-rule traces, stopping iterations and
  fitted parameters within 1e-7.
- **Intrinsic M:** true M = 3, 4, 5 at SNR = 4, one saved pilot dataset each,
  N = 300, D = 60, K = 50. Refit adaptive EB at M = 8 and the original
  similarity-only uniform forward sequence, using Isomap, Spearman
  similarity and single linkage. Require convergence, unchanged selected or
  effective M, ARI loss at most 0.02, and ordering-recovery loss at most 0.005
  relative to the stricter 0.3.1 baseline. Compare each new objective trace
  with the old trace at common iterations (maximum error below 1e-6).

All fits check finite values, normalized probabilities and nondecreasing
fixed-temperature objectives. New pilot fits use the package default
thresholds and at most 3000 post-annealing sweeps. The old strict pilot
baselines used relative 1e-8 stopping. The previously observed true-M = 3
adaptive overestimate is retained as a known baseline outcome.

`results.csv` contains comparisons, timing and pass/fail flags; `run.log`
records execution; `validation.txt` is written only when every check passes.
`results.md` gives the concise summary, and `expected_iterations.csv` verifies
all pilot stopping indices against the frozen old trajectories.
`provenance.txt` records the source commit and input-bundle hash. This small
regression check does not replace the comprehensive intrinsic-M study.
