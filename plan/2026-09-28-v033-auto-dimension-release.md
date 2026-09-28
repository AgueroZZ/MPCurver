# MPCurver 0.3.3 automatic-dimension release

## Goal

Release the package-level automatic ordering-count initialization after
validating it on fixed inputs from the one-monotone-anchor ordering-count
study. Preserve explicit numeric-dimension behavior.

## Method gate

- Use fixed-df natural-cubic-spline variance explained with `spline_r2_df = 5`.
- Use single linkage and maximum mean silhouette among cuts whose clusters all
  contain at least two features.
- Initialize adaptive global ordering probabilities from cluster-size
  proportions without imposing hard feature assignments.
- Confirm automatic and explicit numeric paths, fallback behavior,
  continuation, documentation, and positional API compatibility.

## Release gate

- Pass the package test suite and `R CMD check --no-manual`.
- Rebuild and audit the curated pkgdown site.
- Build and install the 0.3.3 source archive in an isolated simulation library.
- Run the complete 90-dataset paired study extension before publishing the
  research webpage.
- Commit and push package and research repositories only after their respective
  validation gates pass.
