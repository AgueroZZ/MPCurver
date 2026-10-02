# Intrinsic dimension terminology — 2026-10-01

The public interface uses intrinsic-dimension terminology throughout:

- fit_mpcurve(intrinsic_dim = ...) specifies the model dimension or "auto".
- mpcurve_init_control(max_intrinsic_dim = ...) bounds automatic initialization.
- select_mpcurve_dimension(max_intrinsic_dim = ...) bounds dimension search.
- fit$intrinsic_dim reports the effective dimension using the selected
  feature-count-aware numerical occupancy rule.
- fit$model_intrinsic_dim records the retained model dimension, with one latent
  ordering and its trajectory state per dimension.

These replace the reviewed draft's num_orderings and max_num_orderings names.
The fitting entry still has 12 arguments. Reporting remains separate from model
array dimensions. Saved draft objects with num_orderings and dimension-selection
records with max_num_orderings remain readable. No duplicate output aliases are
introduced; continuation returns the current schema.

README, vignettes, reference sources, release notes, examples, tests, and final
argument-audit mappings use the new names. Prior logs retain the naming history;
the current names above supersede those earlier draft names. Version remains
0.4.0. There is no commit, push, tag, or release.

## Verification

- Before/after fits for one dimension, three dimensions, and automatic
  initialization have exactly identical numerical state, parameters,
  responsibilities, locations, and effective counts. Continuation from current
  and saved draft objects is also numerically identical.
- Dimension-search candidate scores, decisions, and selected numerical fit are
  exactly identical to the saved before-change baseline.
- Full local suite: 826 passing expectations; zero failures, errors, warnings,
  or skips. All three public tutorials and README execute with the current
  names, references and Changelog regenerate, and actual Firefox rendering of
  the intrinsic-dimension section is inspected.
- Generated-site audit passes for 31 HTML pages and 1058 local links, assets,
  and anchors. All page version badges remain 0.4.0. git diff --check passes.
