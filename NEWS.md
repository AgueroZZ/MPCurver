# MPCurver 0.3.3

- `fit_mpcurve(intrinsic_dim = "auto")` can now use the feature-similarity
  tree to choose an initialization dimension by mean silhouette. Candidate
  cuts must satisfy `similarity_min_cluster_size` (default 2); if no
  multi-ordering cut is eligible, fitting falls back to one ordering.
- The new fast `similarity_metric = "spline_r2"` uses a fixed-degree-of-freedom
  natural-spline variance-explained score (`spline_r2_df = 5` by default).
  Adaptive fits initialized from an automatic cut start their global ordering
  probabilities at the observed cluster-size proportions. The minimum-size
  rule supplies multi-feature evidence for an initialized ordering; it does
  not classify excluded singleton features as noise or garbage.

# MPCurver 0.3.2

- New fits default to `convergence = "normalized"`: the absolute ELBO change
  per sample-feature entry, `abs(delta_ELBO) / (N * D)`, must be below the
  tolerance, subject to the existing numerical nondecrease checks.
- The default tolerance is `1e-6` for both single-ordering (`tol`) and
  multi-ordering (`tol_outer`) fits. Partition stopping is checked at `T = 1`
  after annealing. A zero tolerance disables early stopping.
- `convergence = "relative"` retains the previous scaling. Reproducing the
  old defaults also requires `tol = 1e-6` and `tol_outer = 1e-5`.
- `do_mpcurve()` inherits the saved convergence rule, with relative stopping
  for objects saved before 0.3.2. An explicit `convergence` argument overrides
  this choice. Multiple initializations and dimension-selection candidates
  use the selected rule. ELBO values and coordinate updates are unchanged.
