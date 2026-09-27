# MPCurver 0.3.2 simulation regression and release

## Authorization

The user requested basic regression checks on the earlier simulations and
authorized pushing the package if those checks show no update-related problems.

## Checks

1. Refit the saved website M = 1 Swiss roll and M = 2 feature-partition examples
   using the same input, initialization seed and iteration budgets. Compare
   recovery and assignments; verify finite, normalized posterior probabilities
   and monotonic fixed-temperature objectives.
2. Refit all 27 simple ND-stopping datasets with the installed 0.3.2 default.
   Compare each trace, stopping iteration and fitted estimates with the saved
   experimental normalized-rule fit from 0.3.1 (absolute tolerance 1e-7).
3. Refit the saved true M = 3, 4, 5, SNR = 4 pilot datasets using adaptive EB
   at M = 8 and the original similarity-only uniform forward protocol. Keep
   data, seeds, K = 50 and annealing fixed. Use normalized 1e-6 stopping and
   compare with the strict relative 1e-8 baselines. Require convergence,
   unchanged effective/selected M, ARI loss <= 0.02, and mean ordering-recovery
   loss <= 0.005. Check new traces against saved old traces at matched sweeps.
   The known M = 3 adaptive overestimate is a baseline result, not a new failure.
4. Save scripts, results and provenance in a dedicated experiment directory;
   preserve original fit artifacts. Run one fitting process with one CPU thread.
5. If checks pass, rebuild the curated package website, audit links and content,
   and commit only this release's source, documentation, validation and records.
   Push main without force and verify the remote head and Pages deployment.

## Scope

This is a release regression check, not the larger intrinsic-M simulation
study. The prior complete package test/check results remain applicable unless
implementation changes are needed. Do not publish unrelated InferOrder work.

## Validation completed

All 47 fresh fits passed; all 18 pilot stopping indices match the frozen
reference predictions. The curated site audit passed (27 pages, 900 local
targets, 139 search entries). See the release log and experiment report.

## Completed

Release `fcc29de29ae8431ced4a3ce38808dd2e47fb8552` was pushed to origin/main.
Pages deployment `36288638786` succeeded, and the live homepage and fitting
reference were verified in Chrome. All requested release steps are complete.
