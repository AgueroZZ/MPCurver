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
