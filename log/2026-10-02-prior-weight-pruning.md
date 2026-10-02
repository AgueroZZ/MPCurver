# Numerical-zero ordering prior weights

The user clarified the intended R1 behavior: describe automatic pruning through
the estimated ordering prior weight, with a very small tolerance divided by the
number of features, and warn when users make the tolerance too large.

For P feature columns, an automatic adaptive fit removes ordering m when its
estimated prior weight omega_m is at or below effective_count_tol / P. The
default numerator remains 1e-8. Adaptive EB gives omega_m = sum_j q(Z_j=m) / P,
so the existing count comparison is retained internally to preserve its exact
floating-point boundary. The argument name and stored fields remain unchanged.

## Changes

- Automatic adaptive partition finalization warns once per fitting or
  continuation call when effective_count_tol exceeds 1e-6. This is an advisory
  boundary on the numerator, before division by P, chosen at 100 times the
  default; it is not a mathematical threshold for model validity. The warning
  explains that large tolerances can remove supported orderings and identifies
  the 1e-8 default. Explicit dimensions, fixed partition priors, and
  single-ordering states do not invoke this warning because they are not
  subject to automatic adaptive partition pruning.
- Before normalizing retained feature-assignment probabilities, finalization
  checks that each feature retains positive finite mass. Unsupported pruning
  now stops with an error identifying effective_count_tol instead of returning
  a model with a NaN objective. No reassignment rule or extra fitting sweeps
  were introduced.
- Fitting/control/object help, README, release notes, and the dimension tutorial
  lead with estimated prior weights. The tutorial presents prior weights first
  in its support table and states the posterior-average identity explicitly.
  Canonical roxygen/R Markdown sources generate the updated Rd and site pages.

The review's R1 example used a prior weight of 0.4 and an equivalent cutoff of
0.45. It therefore tested an oversized setting rather than removal of a
numerical-zero prior. The default criterion remains the same. R2 and the other
review findings await separate discussion and are not repaired by this change.

## Verification

- Complete source suite: 1,050 passing expectations, zero failures, errors,
  warnings, or skips. Focused finalization checks pass 166 expectations.
  Tests cover the exact 1e-6 warning boundary, constructor/named-list controls,
  warning scope and multiplicity, continuation, and both soft-underflow and
  hard-assignment R1 reproductions. Expected warning/error conditions are
  captured explicitly. Logs: `/tmp/mpcurver-r1-full-tests.log` and
  `/tmp/mpcurver-r1-focused-tests.log`.
- Before/after default single, specified two-ordering, and automatic fits, their
  continuations, and a numerical-residue pruning fixture are identical as
  complete R objects. Evidence: `/tmp/mpcurver-r1-baseline.R`,
  `/tmp/mpcurver-r1-baseline.rds`, and `/tmp/mpcurver-r1-comparison.log`.
- README and the edited dimension tutorial execute during regeneration. The
  existing pkgdown builder's normalization and pinned documentation tools are
  used for the affected home, reference, article, news, search, and sitemap
  outputs. Roxygen2 7.3.2 emits the existing advisory relative to RoxygenNote
  7.3.3; only the three affected Rd topics are regenerated.
- Browser checks on six affected pages confirm the wording, prior-weight-first
  table, rendered equations, zero KaTeX errors, and no page-width overflow at
  a 1053-pixel viewport. Screenshots of the edited formula/table and control
  description are inspected. Browser evidence is under `/tmp/mpcurver-r1-*`.
- Site audit passes all 33 HTML pages and 1,130 local links/assets/anchors,
  including the existing sequential figure/table captions. Whitespace checks
  pass. The earlier full package check in the review remains historical;
  it was not repeated for this focused patch.

No commit, push, or publication was performed. Historical scientific results
remain unchanged.
