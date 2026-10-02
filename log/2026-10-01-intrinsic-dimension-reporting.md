# Intrinsic dimension reporting — 2026-10-01

The public intrinsic dimension now denotes the effective number of orderings.
Previously, `intrinsic_dim` denoted the fitted ordering slots and
`effective_intrinsic_dim` denoted the effective count. New fits and summaries
report that count once, as `intrinsic_dim`; `num_orderings` records retained
ordering slots. The duplicate `effective_intrinsic_dim` output is removed.

The estimator is unchanged. With an adaptive partition prior, the count is the
number of fitted global ordering weights strictly above `effective_weight_tol`
(default `1e-12`). With a fixed prior, the count equals the retained ordering
count. This is a summary of fitted weights, rather than a posterior over the
number of orderings. No trajectory slots are removed by reporting the count.

Conversion, printing, summaries, trajectory plotting, and continuation distinguish
these quantities. Saved objects with the old schema remain readable; summaries
use their effective count, while plotting and continuation use their saved model
slots. Continuation returns the current schema. Legacy requested/active/displayed
metadata remains available for compatibility but is not used as the estimate.

The partition tutorial result table now contains one intrinsic-dimension row.
README, object/reference documentation, and release notes describe the meaning
and the field migration. Website changes are generated from canonical sources.

## Verification

- Full local suite: 813 passing expectations; zero failures, errors, warnings,
  or skips. Tests cover two retained slots with reported dimension one, plots,
  summaries, and identical continuation state from current and old schemas.
- Before/after comparison for one and three ordering slots: underlying fit,
  parameters, responsibilities, sample locations, and effective count are exactly
  identical. Only the public reporting schema changes.
- All three public articles execute successfully through the site builder.
  Generated-site audit passes for 30 HTML pages and 996 local links, assets,
  and anchors. Numbered figure/table captions are checked on all three articles.
- Corrected a duplicated figure number found during visual inspection. The
  partition article now has four uniquely numbered figure captions and one
  numbered table, with matching text references. This final caption-only edit
  was rendered and inspected after the package-check archive was built.
- Actual Firefox rendering of the revised tutorial/table is inspected. Local
  review images are stored under `/tmp/`; they are not package artifacts.
- `git diff --check` passes. Full `R CMD build` and `R CMD check --no-manual`
  finish with `Status: OK`: zero errors, warnings, or notes. Installed tests
  have 795 passing expectations and the standard CRAN snapshot skip. Check
  artifacts are under `/tmp/mpcurver-dimension-check-20261001/`. All current
  R source files match the checked archive exactly; the subsequent change
  affects only tutorial figure captions and their generated webpage.

All changes remain local and uncommitted. The version is still 0.3.4 pending the
user's release-version decision; no commit, push, or release is performed.

## Subsequent definition correction

The user subsequently selected a feature-count-aware threshold. The current
rule counts total soft feature assignments above 1e-8, with control named
`effective_count_tol`. The prior-weight threshold described earlier in this log
is superseded. See `2026-10-01-p-aware-intrinsic-dimension.md` for the definition,
compatibility behavior, and current verification.
