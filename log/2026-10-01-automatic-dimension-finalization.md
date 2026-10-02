# Return the estimated model dimension after automatic fitting

## Behavior

Current fits use intrinsic_dim as the actual number of orderings in the
returned model. An explicitly specified integer M retains that dimension.
Automatic fitting initializes M from the feature-similarity tree. Under an
adaptive partition prior, finalization removes orderings whose expected
feature count N_m = sum_j q(Z_j = m) is at or below effective_count_tol
(default 1e-8). The equivalent empirical-Bayes prior-weight cutoff is
1e-8 / P. A fixed partition prior retains the initialized dimension.

R/14_dimension_finalization.R performs this finalization without any additional
CAVI sweeps. It subsets the complete structural state, preserves the retained
conditional trajectory posteriors, positions, smoothness parameters, noise
variances, and measurement errors, normalizes remaining assignment probabilities,
and recomputes prior weights and the objective. A reduction to M = 1 returns
the ordinary cavi/mpcurve single-ordering schema and can be plotted or continued.
The model_intrinsic_dim compatibility field equals intrinsic_dim for current
fits. Saved objects with the earlier effective-count schema remain readable.

## Provenance and continuation

Automatic fits include dimension_estimation with initial and final dimensions,
the threshold, expected feature counts, and pruning events. Each event records
retained/removed labels, pre-removal counts and prior weights, objective before
and after removal, refitted = FALSE, and the preceding fitting history. Those
traces remain associated with their original dimension. When pruning occurs,
the returned model's trace starts at its recomputed objective. Iteration offsets
preserve the total number of fitting sweeps through subsequent continuation.
Continuation retains the estimation record and applies the same finalization
to any newly empty orderings.

This updates the previous reporting-only definition: negligible orderings were
formerly retained while intrinsic_dim reported a smaller effective count. The
new behavior removes them for automatic adaptive fits, and explicit dimensions
are now reported as specified. Historical checks and outputs remain historical.

## Documentation

Updated fitting/control/object/summary/continuation roxygen help, README source,
NEWS, and the two dimension-related tutorial sources. The estimation tutorial
reports one dimension in Table 1. Initialization history remains accessible in
process records. Ordinary print and summary output present one dimension for
current fits. Anneal_steps and anneal_sweeps remain unchanged as requested.

Regenerated README, the home page, affected reference pages, all four public
articles, news, search, and sitemap through the existing pkgdown builders.
The current tutorial still returns M = 3: its three expected feature counts
are 6, 3, and 3, so no numerical-residue ordering is removed in that example.
The independent uniform-prior search still selects M = 2.

## Current verification

- The complete package suite passes 1028 expectations with zero failures,
  errors, warnings, or skips. New tests exercise removal to two and one
  orderings, exact-zero and positive numerical residue, the strict cutoff,
  fixed-prior and explicit-M retention, single-ordering initialization,
  multiple pruning events, and continued fitting. Noise cases include estimated
  variances, feature-level known standard deviations, and observation-level
  matrices. Tests verify unchanged retained posteriors/parameters, probability
  normalization, an independently reconstructed objective, single-ordering
  score agreement, ordinary plotting, metadata, and iteration counts.
- Public fitting is tested with an intentionally larger cutoff to trigger
  pruning at initialization and confirm that no further fitting is invoked.
  Default-threshold residue cleanup is tested on complete structural fixtures.
- Loaded browser checks on the estimation article and fitting/control help
  confirm the single-row dimension table, section/TOC order, updated definitions,
  zero KaTeX errors, and no page-width overflow at 1053 pixels. Actual table
  layout and equations are visually inspected.
- The site audit passes 32 pages and 1122 local links/assets/anchors, with
  sequential numbered figure/table captions and references on all four public
  tutorials. Whitespace checks pass. The preceding package-check Status OK is
  historical and was not rerun for this change.

No commit, push, tag, or release.
