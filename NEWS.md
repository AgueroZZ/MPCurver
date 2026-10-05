# MPCurver 0.4.1

- Isomap now defaults to `num_neighbors = NULL`, selecting the smallest
  neighbor count that connects the full undirected graph (`k_min`) separately
  within each feature group. Explicit counts remain fixed; counts of one are
  supported. Ordering results and fit initialization metadata record `k_used`.
  Neighbor queries exclude each sample itself and break distance ties by row
  index, keeping graphs nested during the connectivity search.

# MPCurver 0.4.0

- Added `fitted_positions()`, `fitted_trajectories()`, and
  `fitted_assignments()` with consistent dimensions across ordering counts,
  posterior uncertainty extraction, and preserved sample/feature names.
- Plotting accepts feature names, uses named axis labels, and adapts its default
  to one-feature fits. Invalid numeric/nonfinite data and measurement SDs are
  rejected before initialization, with locations for nonfinite observations.
- Tutorials now present models and software examples; single- and multi-ordering
  CAVI algorithms have separate pages. The README starts with a small runnable
  workflow. Dimension guidance recommends adaptive fitting first when the
  dimension is unknown and explains the fixed-uniform comparison alternative.

- Removed `hard_assign_final`: fitting and continuation retain posterior soft
  feature assignments. Most-probable labels remain available for reporting.
- Default position grids now contain at least `rw_order + 1` bins. Automatic
  grids report a clear error when there are too few samples for that order.
- Pseudotime scatterplots use a common absolute `[0,1]` color scale. Plot labels
  and other named graphics arguments consistently override defaults.
- Fixed scalar sinusoid-frequency sampling and singleton GP feature blocks.
- Spiral manual examples load the package checkout and use `fit_mpcurve()`
  instead of sourcing selected internal implementation files.

- Invalid ordering-helper arguments and fractional PCA component indices now
  error consistently across single, grouped, and automatic initialization.
  Ordering computation failures warn before falling back to PCA component 1;
  `mpcurve_init_control(on_failure = "error")` disables this fallback. A failed
  PCA attempt stops fitting. Single-feature groups run the requested helper,
  and initialization metadata preserves the actual method and fallback cause.
- Redesigned `fit_mpcurve()` around 12 arguments: `X`, `S`, `num_bins`,
  `intrinsic_dim`, `initial_method`, `position_prior`, `partition_prior`,
  `max_iter`, `tol`, `verbose`, `init_control`, and `control`.
- Added validated `mpcurve_control()` and `mpcurve_init_control()` constructors.
  Model settings include `rw_order`, `lambda_init`, `fix_lambda`, combined
  `lambda_bounds` and `sigma2_bounds`, and `anneal_steps`, `anneal_start`,
  `anneal_sweeps`.
- Disabled annealing by default (`anneal_steps = 0`): multi-ordering fits run
  directly at temperature 1 within the `max_iter` budget. A positive
  `anneal_steps` enables the optional schedule, which finishes at temperature 1.
- Fixed scale sensitivity in `spline_r2` similarity. Nonconstant features are
  centered and standardized to unit sample variance before spline scoring;
  numerical guards operate in standardized units. Constant columns have zero
  off-diagonal similarity. Residual sums of squares are computed directly from
  standardized residuals. The data used for model fitting are unchanged.
- Initialization now accepts exactly one method for any ordering count. Removed
  resource-dependent candidate selection and `num_cores`; compare initializations
  using separate fitting calls. Ordering-helper settings use `method_args`.
- Dimension selection fits each candidate dimension once, honoring the chosen
  `initial_method`. Removed the automatic
  comparison of similarity-based and complete-matrix initializations.
- Removed `partition_init` and the complete-matrix `"ordering_methods"`
  strategy. Multiple orderings are initialized from feature-similarity groups,
  applying the single selected `initial_method` within every group. The default
  PCA uses each group's own PC1.
- Unified `tol` across single-ordering and partition fits, including
  `do_mpcurve()`. `max_iter` replaces the ordinary fitting caps; partition
  optional annealing has a separate advanced budget.
- Removed `algorithm`, `greedy`, the ignored freeze/drop arguments, and the
  exploratory Dirichlet assignment model (`assignment_prior`, `ordering_alpha`).
  Dimension selection uses the new interface. Continuation no longer accepts
  the retired freeze/drop or assignment-prior arguments.
- Automatic adaptive fits remove orderings whose estimated prior weights are
  at or below `effective_count_tol / ncol(X)`. The tolerance defaults to `1e-8`,
  replacing `effective_weight_tol`, and identifies numerical-zero prior weights.
  Automatic adaptive pruning warns for tolerances above `1e-6` and stops with
  an error if pruning leaves a feature with no positive remaining assignment
  probability. Remaining probabilities are normalized and
  the objective is recomputed without further CAVI sweeps. A result reduced
  to one ordering uses the ordinary single-ordering schema. Estimation and
  pre-removal fitting records are preserved in `dimension_estimation`.
- Unified dimension reporting: `intrinsic_dim` is the actual number of
  orderings in the returned model; an explicitly specified dimension is kept.
  `model_intrinsic_dim` remains a compatibility alias for current fits.
- This is a breaking API change: retired named fitting arguments are rejected.
  For example, replace `K = 20, rw_q = 3, lambda = 5, fix_lambda = TRUE` with
  `num_bins = 20, control = mpcurve_control(rw_order = 3, lambda_init = 5,
  fix_lambda = TRUE)`. The returned grid-size field `K` is retained.

- Simplified `do_mpcurve()` to `object`, `max_iter`, `tol`, `verbose`, and
  `control`. Added `mpcurve_continue_control()`, whose omitted settings inherit
  fitted values. `lambda_init` resets precision and `fix_lambda` selects its
  update mode consistently for one or multiple orderings. Continuation honors
  both precision and variance bounds, including fixed precisions, and retains
  initialization and continuation provenance. Measurement errors are inherited.
- Dimension selection now exposes fitting settings explicitly and uses
  `max_intrinsic_dim`; comparisons retain fixed uniform partition priors.
- Ordering helpers use `num_neighbors`, `max_iter`, and `tol` consistently.
  Method-specific advanced options are documented under a named `control`
  list. Isomap's stabilization constant is internal. PCA's concise interface
  is retained, with strict integer component validation.
- Simulation helpers use `n`, `d`, `num_bins`, and `noise_sd`; random-walk,
  trajectory-shape, and Gaussian-process options use named `control` lists.
  Default generating distributions and truth fields are preserved. Unknown
  options and invalid counts are rejected. Plotting and summary interfaces
  remain concise; `fitted_prior()` rejects unused extra arguments.
- Renamed `simulate_cavi_toy()` to `simulate_mpcurve()`: data are generated
  from a single-ordering MPCurve model; CAVI names the fitting algorithm.
  The generating distribution and returned fields are unchanged.

# MPCurver 0.3.4

- Similarity-based partition initialization now applies a length-one ordering
  method independently within every selected feature block. The default PCA
  initialization uses each block's own PC1 instead of assigning successive
  components (PC1, PC2, ...) across unrelated blocks. Explicit per-block
  methods and `pca_components` remain supported, and
  `partition_init = "ordering_methods"` retains its successive-PC behavior on
  the common full feature matrix.

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
