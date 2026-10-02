# Similarity-only initialization — 2026-10-01

Removed `partition_init` from `mpcurve_init_control()` and the complete-matrix
`ordering_methods` strategy from the fitting initialization chain. Multiple
orderings now use feature-similarity clustering and apply the single selected
`initial_method` independently within each feature group. Default PCA uses
each group's own PC1. A method and a feature-grouping strategy are different
concepts; the redundant strategy selector is now retired.

Updated the public wrapper, private fitting dispatcher, structural initializer,
shared fixed-dimension initialization helper, and dimension selector. Removed
the complete-matrix branch, sequential-PC defaults, and its unused validation
helpers. The archival augmented implementation shares the remaining grouping
helper. A fixed `partition_init = "similarity"` provenance field remains in
fitted control metadata; it is not a user option. Supplied warm-start fits and
within-group PCA-component controls retain their existing advanced behavior.

Updated active tests to use the supported initialization, with regression
coverage rejecting the retired option in the constructor, named control lists,
single/multiple/automatic fitting, and dimension selection. Updated help,
the estimation tutorial, and the 0.4.0 changelog. Historical changelog entries
and experiment outputs retain their original provenance.

Current verification:

- Single-ordering, fixed multiple-ordering, and automatic default fits are
  identical to saved pre-change objects, including their numerical states
  and histories.
- Complete suite: 1036 passing expectations; no failures, errors, warnings,
  or skips.
- All 22 exported functions and public S3 display methods have documented
  argument names and defaults matching the current function definitions.
  Their signatures contain no retired public arguments, and their help has
  no obsolete partition initialization strategy.
- Regenerated all 23 curated reference topics, the reference index,
  executed intrinsic-dimension article, news, search, and sitemap in place.
- Browser checks confirm the current initialization usage and grouping
  descriptions on the initialization, fitting, selection, and estimation
  pages, with no math errors or overflow at a 1053-pixel viewport.
- Site audit passes 33 pages and 1124 local links, assets, and anchors;
  existing figure and table numbering and whitespace checks pass.

No commit or push.
