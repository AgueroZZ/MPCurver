# MPCurver 0.3.4 common per-partition ordering initialization

- Agent: Codex
- Date: 2026-09-28 UTC
- Request: update the default partition initialization so one selected ordering
  method is applied independently to every feature partition.

## Software behavior

For `partition_init = "similarity"`, a length-one `method` is now repeated
across all selected feature blocks and fitted independently within each block.
The default remains PCA, and every block uses its own first principal component.
Thus the default no longer assigns successive components PC1, PC2, ..., PCM to
different blocks.

Explicit method vectors and explicit `pca_components` are still honored. The
independent-block rule applies equally to scalar Isomap, Fiedler, principal
curve, random, and other supported ordering methods. The distinct
`partition_init = "ordering_methods"` path retains successive PCA components by
default because all candidate orderings there use the same full feature matrix.

## Documentation and generated site

The package version is 0.3.4. NEWS, README source and rendered Markdown,
`fit_mpcurve()` roxygen/Rd documentation, and the partition vignette describe
the corrected semantics. The curated pkgdown site was rebuilt with pkgdown
2.2.0 and fontawesome 0.5.3, matching the previous site's asset versions. The
relevant API/article text and four changed partition figures were inspected;
the transparent-background objective figure was additionally checked to be
nonempty.

The site audit passed for 27 HTML pages, 903 local links/assets/anchors, and 138
search entries, with no errors.

## Validation

- Source-loaded test suite: zero failures or warnings and four existing skips.
- Packaged `R CMD check --no-manual`: 644 passing expectations, zero failures or
  warnings, and five expected skips (the four existing retired-semantics skips
  plus one CRAN-only skip).
- Complete package check: zero errors, zero warnings, and one environment NOTE
  because optional Suggests `devtools` and `fields` were unavailable. Examples,
  packaged tests, vignette execution, and vignette rebuild all passed.
- `git diff --check` passed.
- Source archive: `/tmp/MPCurver_0.3.4.tar.gz`, SHA-256
  `58d0c99c6170994c82eedba190fb1c28a63eb47530c575705323c3cf3992b7c6`.

## Publication verification

The release commit is `15f2b0bbe5dfa61cd46da5160b2bc251e75a0475` on
`origin/main`. GitHub Pages deployment
[36499279949](https://github.com/AgueroZZ/MPCurver/actions/runs/36499279949)
completed successfully for that exact commit. The live `fit_mpcurve()`
reference and partition article both show version 0.3.4 and document that a
length-one method is applied independently to every selected feature block;
the live pages also state that default PCA uses each block's own PC1.

## UTF-8 site correction

The fitness article was rebuilt under a valid UTF-8 locale after the initial
deployment was found to render the micro symbol as malformed HTML tags in
seven environment labels. The corrected article preserves `μg/ml` and `μM`.
`scripts/build_public_site.R` now refuses a non-UTF-8 locale and checks the
rendered fitness article for the malformed tag pattern before completing.
This correction changes generated documentation only; the validated 0.3.4
package source archive and statistical results are unchanged.
