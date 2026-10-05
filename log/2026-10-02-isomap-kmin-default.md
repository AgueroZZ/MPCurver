# Isomap defaults to the smallest connected neighborhood

At the user's request, changed Isomap's automatic neighborhood from a fixed
15 to `k_min`, the smallest integer from 1 through `nrow(X)-1` connecting all
samples in the undirected union-kNN graph. The public `num_neighbors` default
is now `NULL`; an explicit integer remains fixed. This applies to Isomap
initialization within each feature group independently. The package development
version is 0.4.0.9000, based on commit `f1511a0`.

Connectivity search brackets the first connected graph by doubling the count,
then finds the exact minimum by binary search. Only the selected graph is
embedded, and only one MPCurve fit is run. Exact Euclidean neighbor queries
exclude the query sample explicitly. Equal-distance neighbors are ordered by
sample row index, expanding the query to include every boundary tie before
truncation; the resulting graphs are nested even with duplicate observations.
Ordering results expose `k_used`, and fitting carries realized graph counts
into single/grouped initialization metadata and through continuation.

The change is motivated by the fixed 300-by-12 B-group example from
`InferOrder/experiments/isomap_elbo_screen_v040/`. Connectivity selects k=2.
With the original 50-bin RW2 controls and seed 261925351, the new automatic
default exactly matches the saved explicit k=2 fit: one-sweep ELBO
-3916.29106229426, final ELBO -3444.4225542047, absolute Spearman recovery
0.996795964399604, and convergence after 176 sweeps. ELBO history, positions,
noise variances, precisions, and position probabilities have zero maximum
absolute difference. The earlier k=15 recovery is 0.75861420682452. The selected
B-group matrix SHA256 is
`e1ce6e64c93b0a51e31cb45888dfdd17e5a7daae5ba21d95c7ad86990e892afa`.

Connectivity does not establish optimal neighborhood geometry or the best
converged ELBO. Here k=2 recovers positions as well as the narrow-range
one-sweep selection k=4, with a final ELBO 1.814885 below it. Generalization
across trajectories and sampling densities remains an empirical question.
The graph upper-bound heuristic and early-ELBO screening remain internal
experiments; this change uses the lower-bound connectivity rule alone.

## Verification and reproduction

- `tests/testthat/test-isomap-connectivity.R` checks the minimum against an
  independent exhaustive pairwise-distance reference on connected, separated,
  random, and three-sample inputs. It also checks all neighborhoods on tied
  grids and duplicate observations, explicit disconnected-graph handling,
  group-specific counts, and continuation metadata. Invalid-count tests now
  accept `NULL` and reject zero, fractional, and out-of-range counts.
- The full source test suite passes 1,525 expectations with no failures,
  errors, warnings, or skips. Evidence: `/tmp/mpcurver-kmin-tests.log/.rds`.
- Run `Rscript tests/manual/isomap_default_kmin.R` from the package root for
  the unchanged-input historical comparison above. This requires the adjacent
  InferOrder input and saved reference, `pkgload`, and `digest`. Evidence:
  `/tmp/mpcurver-kmin-case.log`.
- README, roxygen reference pages, and the curated six-article pkgdown site
  execute against the current checkout. The site audit checks all 38 HTML
  pages and 1,403 local links/assets/anchors. Five affected pages pass browser
  checks for math errors, image loading, and page-width overflow. Shortened
  the fitness example's four-panel plot titles after inspecting their rendered
  clipping; this is a display change. Evidence: `/tmp/mpcurver-kmin-docs.log`,
  `/tmp/mpcurver-kmin-site-audit.json`, and browser logs/screenshots.
- `R CMD check --no-manual` finishes with **Status: OK**, without errors,
  warnings, or notes; examples and all six vignette code runs/rebuilds pass.
  All 22 runtime R files, 44 Rd files, and NAMESPACE match the checked archive.
  The archive predates the final fitness plot-title and default-description
  display edits; the revised article was separately executed and rebuilt,
  its nine figures visually inspected in the browser, and the final site
  audited again. Evidence: `/tmp/mpcurver-kmin-check/`,
  `/tmp/mpcurver-kmin-fitness-docs.log`, and the final figure browser records.

Temporary checks, installed check packages, archives, and browser diagnostics
are under `/tmp/mpcurver-kmin-*`. Historical InferOrder inputs, results, fitting
libraries, and public research pages are preserved. No commit or push.
