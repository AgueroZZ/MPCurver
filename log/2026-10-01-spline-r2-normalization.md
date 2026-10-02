# Make spline R-squared similarity independent of feature units

## Problem and scientific intent

Spline R-squared initialization previously screened responses using an
absolute total-sum-of-squares threshold of approximately 1.49e-8, in squared
input units, in addition to the raw feature-SD threshold. For 20 samples and
proportional features `t` and `2*t`, scaling both by 1e-5 changed similarity
from 1 to 0 despite preserving the relationship and passing the SD guard.

The user proposed unit-variance normalization before R-squared calculation.
R-squared is dimensionless: centering and positive rescaling preserve its
population/least-squares definition while making its numerical calculation
independent of feature units.

## Implementation

- Added `.cavi_standardize_similarity_features()` in `R/10_partition_cavi.R`.
  Nonconstant columns are centered and divided by their sample SD. Centered
  values are first divided by their maximum absolute magnitude, allowing
  variance normalization without squaring extremely small or large input
  values. This supports tested unit changes from 1e-200 to 1e200.
- The spline score uses this standardized copy. For each usable column,
  sample variance is 1 and total sum of squares is n - 1. Removed the absolute
  raw-SST eligibility threshold and applied the internal SD guard in
  standardized units. Constant/unusable columns receive zero off-diagonal
  similarity; diagonal entries retain the existing value of one.
- Preserve predictor ordering from the original feature values. Compute SSE
  directly from the standardized spline residuals, avoiding cancellation from
  subtracting two large sums of squares. Retain the existing fixed-df basis,
  directional-score clipping, and maximum-of-two-directions symmetrization.
- Model fitting continues to use the original X and supplied S. Existing
  `feature_info$sd` records raw-unit SD; `standardized_sd` records the normalized
  SD. Similarity results and partition initialization record
  `centered_unit_variance` as their normalization convention.
- This correction applies to `spline_r2`. Correlation and smooth-fit similarity
  rules are unchanged. Saved fits and frozen experiment libraries are preserved.
- Updated initialization reference documentation, the partitioning tutorial,
  and development NEWS from canonical sources. No ordinary option is added.

## Current verification

- 731 passing test expectations; zero failures, errors, warnings, or skips.
  Added regressions for the original proportional-feature example, scales
  from 1e-200 to 1e200, independent positive feature rescaling, large offsets,
  constant columns, normalized sample variance/SST, and original-data retention
  in an automatic partition fit.
- The proportional example now has similarity 1 at both original and 1e-5
  scales. The score remains 1 below the former raw-SD threshold as well.
- On a seed-903 ordinary-scale 40-by-10 Gaussian matrix, the largest directional
  score difference from the preceding implementation is 3.330669e-16. All 100
  directional scores agree within 1e-12 with independent `lm.fit()` residual
  R-squared calculations.
- Rebuilt initialization Rd/reference examples and rendered the partition
  article. Inspected the resulting standardization paragraph and retained
  figure/table captions; all 979 local links/assets/anchors across 29 HTML pages
  resolve. `git diff --check` passes.
- The earlier full package check applies to the interface-refactor snapshot.
  Current checks for this scientific correction are the full regression suite,
  independent score comparisons, and relevant documentation builds above.

The previously recorded spline scale-sensitivity issue is resolved for finite,
representable nonconstant features. Scaling/offsets that erase variation through
input rounding cannot recover information absent from the supplied numbers.
No commit, push, or publication was performed.
