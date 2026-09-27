# MPCurver 0.3.2: normalized ELBO convergence

## Requested outcome

Make ELBO change per observed matrix entry the default stopping criterion,
record the supporting experiment, and update package version and documentation.

## Design

- Add `convergence = c("normalized", "relative")` to `fit_mpcurve()`.
  Normalized change is `abs(ELBO_new - ELBO_old) / (N * D)`, where rows are
  samples and columns are features. Default `tol` and `tol_outer` are `1e-6`.
- Keep existing objective-decrease guards and restrict partition stopping to
  the final `T = 1` phase. Do not change coordinate updates or ELBO traces.
- Preserve the previous relative denominators (`abs(previous) + 1` for one
  ordering; `abs(previous) + 1e-12` for partition fits).
- Store the rule in fit controls. Continuation inherits it; objects without
  the field use the legacy relative rule unless explicitly overridden.
- Forward the rule through multiple initializations, dimension selection,
  and the two-ordering compatibility wrapper.

## Implementation and validation

1. Update active fitting and continuation paths, preserving positional APIs.
2. Test stopping against independent trace calculations for estimated and
   known noise, single and multiple orderings, continuation and selection.
3. Bump DESCRIPTION to 0.3.2; update help, README, vignette and NEWS.
4. Run focused and full tests, package checks, and install a local build.
5. Record results and validation in `log/`; retain frozen 0.3.1 experiments.

## Supporting evidence

InferOrder `experiments/elbo_nd_stopping_v031/` contains the reproducible
comparison on 27 M = 1 datasets (three seeds over N = 100, 200, 400 and
D = 10, 20, 40). At tolerance 1e-6, relative versus normalized stopping took
50.753 versus 46.632 seconds, with median 633 versus 591 iterations. All
normalized fits used fewer iterations, and 23/27 took less time. Ordering
correlation with the stricter reference was at least 0.9999625. This supports
the default for these simple examples; it does not establish speed or
model-selection accuracy across all multi-ordering problems.

## Completion

All five steps completed on 2026-09-26. Full tests passed (578 assertions;
four existing legacy skips). Package check reported 0 errors, 0 warnings and
one existing NOTE. Version 0.3.2 was installed and verified in a fresh R
process. See `log/2026-09-26-codex-normalized-elbo-convergence.md`.
