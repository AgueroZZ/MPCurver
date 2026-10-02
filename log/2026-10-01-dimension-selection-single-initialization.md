# Dimension selection uses one initialization per candidate

## Correction

The redesigned public interface accepted one initial_method, but the
selector still fit every M > 1 twice: similarity-group initialization and
complete-matrix ordering-method initialization. It overrode the supplied
init_control$partition_init and selected the larger finite objective. This
legacy behavior conflicted with the agreed single-initialization interface.

## Change

R/11_dimension_selection.R now fits each candidate dimension exactly once,
passing the validated initialization settings unchanged to fit_mpcurve().
The selected initial_method and init_control$partition_init apply throughout
the search. Defaults use similarity grouping and PCA within each group.
A failed candidate reports its selected strategy and error; it does not
retry another initialization. The common-grid check, fixed uniform prior,
soft temperature-one objective, acceptance rule, and continuation metadata
are retained. Candidate records retain selected_for_M for compatibility;
the introductory table omits that redundant column.

Updated the selector's roxygen help, NEWS, and the canonical intrinsic-
dimension tutorial. Regenerated the Rd/reference page, executed article,
changelog, and search index. Historical results in earlier logs remain
historical; this is a correction to the current selector behavior.

## Current verification

- The dimension-selection test file passes 86 expectations with zero failures,
  errors, warnings, or skips. Tests check one record per candidate in both
  forward and backward searches, scores against direct fits, the explicitly
  selected complete-matrix strategy, the fiedler method with similarity
  grouping, known measurement errors, and fixed precision. Existing checks
  cover continuation metadata, uniform priors, bounds, and invalid controls.
- The executed tutorial still selects M = 2. Its M = 1, 2, 3 candidate scores
  are 1179.850, 2353.147, and 2136.357, respectively. The first two candidates
  report convergence; M = 3 reaches its iteration limit. There is one fit per
  dimension. The former M = 2 score of 2389.838 came from the removed automatic
  comparison and is superseded by the default-strategy score above.
- Loaded browser checks on the tutorial and selector reference confirm the
  single-fit description, unique candidate dimensions, zero KaTeX errors,
  and no page-width overflow at 1053 pixels. Actual equations and tables are
  visually inspected.
- The generated-site audit passes 32 HTML pages and 1121 local links/assets/
  anchors, with figure/table numbering and reference checks. git diff --check
  passes. The preceding 856-expectation full suite is historical; it was not
  repeated for this focused selector change.

No commit, push, tag, or release.
