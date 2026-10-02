# Default fitting without annealing — 2026-10-01

The default `mpcurve_control()` now uses `anneal_steps = 0L`. Multi-ordering
fitting starts at temperature 1 and performs only the joint coordinate sweeps
allowed by `max_iter`, subject to the existing tolerance rule. Initial objective
and assignment information are evaluated at temperature 1; the annealing
schedule is empty. This removes the previous 25 extra annealing sweeps from
ordinary fits. Initialization itself is unchanged.

Annealing remains an explicit advanced choice, for example
`mpcurve_control(anneal_steps = 25, anneal_start = 5)`. The start temperature
and sweeps per temperature apply only when the step count is positive. Enabled
annealing precedes the `max_iter` phase and preserves its earlier numerical
behavior. The structural backend and dual-ordering wrapper share the new
default. The unused augmented legacy engine retains its historical defaults.

The existing `n_anneal` backend field includes the initial objective record;
its value is 1 when annealing is disabled, representing the initial trace index
and zero annealing sweeps. Continuation uses temperature 1 as before. No new
public argument or output alias is introduced.

Roxygen sources, the partition and introductory tutorials, release notes, and
the final argument-audit mappings explain the default and the opt-in controls.
Tutorial outputs are regenerated from current source. Historical research fits
and earlier verification logs retain their original settings and provenance.

## Verification

- Full local suite: 856 passing expectations; zero failures, errors, warnings,
  or skips. New checks cover zero-step fitting, exact iteration budgets,
  known and estimated noise, inactive annealing controls, continuation matching
  uninterrupted fitting, and explicit annealing budgets and temperature traces.
- An explicitly enabled three-temperature, two-sweeps-per-temperature fit is
  completely identical to the saved before-change result, including all
  numerical state and metadata. Default fits intentionally follow the new
  schedule and can differ from earlier default fits.
- README and all three public tutorials execute against current source;
  reference pages, Changelog, search index, and sitemap regenerate through the
  existing builder. Browser inspection confirms the new tutorial paragraph,
  convergence example, objective and trajectory figures, and control reference.
- Generated-site audit passes for 31 HTML pages and 1056 local links, assets,
  and anchors; all version badges are 0.4.0. Figure and table numbering and
  captions pass the existing audit. `git diff --check` passes.
- The earlier full package check remains historical verification; it was not
  repeated for this change. The current suite and executed documentation are
  the verification reported above.

Version remains 0.4.0. No commit, push, tag, or release.
