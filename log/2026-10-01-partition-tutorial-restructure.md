# Multiple-ordering and dimension-estimation tutorials — 2026-10-01

The multiple-ordering tutorial keeps its existing `articles/partition.html`
URL and now has three main sections in the requested order:

1. Intrinsic dimension: why one sample ordering can be insufficient and how
   multiple feature groups motivate latent orderings.
2. Statistical model: the standard fixed-M Gaussian model, intrinsic RW2
   trajectory prior, structured variational family, conditional trajectory,
   sample-position, and feature-assignment updates, empirical-Bayes parameter
   updates, and the default temperature-1 stopping rule.
3. Fit a model with M = 2: the executable simulation and fit, feature
   assignments, trajectories, result extraction, convergence, continuation,
   and optional fitting settings.

The fitting example uses the ordinary defaults with `intrinsic_dim = 2` and
50 bins. The local assignment-score expression and known-measurement-error
formulas remain available in expandable details. Sample-position bin counts
use n_mk, distinct from the expected feature counts N_m on the estimation page.
The one-feature plot uses a Viridis palette and its pseudotime axis instead of
a redundant gradient legend; the caption specifies colors and fitted means.
The model-results table records `model_intrinsic_dim`, with a link to effective
dimension reporting on the separate page.

The new `vignettes/intrinsic_dimension.Rmd` / `articles/intrinsic_dimension.html`
tutorial explains automatic silhouette-based initialization, adaptive expected
feature counts and the existing P-aware threshold, and uniform-prior adjacent
dimension selection. It executes both approaches on a reproducible 200-sample,
12-feature simulation with two generating orderings, seed 1, twenty bins, and
200 fitting iterations. Adaptive fitting retains and reports three orderings;
uniform-prior forward selection chooses two. Candidate objectives and
convergence flags are shown, including candidates that reach the iteration
limit. The page explains the different criteria and the need to check fit
quality. These are tutorial outputs, not a new method-comparison study.

Four numbered tables report the results. Reporting-only table construction
is hidden, while fitting and result-access calls remain visible. Wide tables
use responsive containers and descriptive column labels. README, the getting
started tutorial, article navigation, and the curated public-site builder link
to the separate estimation page. All four public tutorials are included in
the builder. Package fitting code, dimension definitions, and defaults are
unchanged; version remains 0.4.0.

## Verification

- README and all four public tutorials execute through the existing builder.
  The estimation page is subsequently rendered for the final table layout.
- Browser checks confirm the requested main-section and live table-of-contents
  order on both pages, no KaTeX errors, and no document-width overflow at a
  1053-pixel viewport. The early blank TOC screenshot preceded completion of
  navigation scripts; loaded navigation works and the builder already performs
  site initialization. No extra initialization fix was needed.
- Final browser inspection confirms readable model equations, working main
  and nested navigation, the M=2 example, the trajectory figures and their
  captions, and the compact dimension-comparison tables without clipped
  headers. Reporting-only construction code is absent from the estimation
  page.
- Generated-site audit passes for 32 HTML pages and 1120 local links, assets,
  and anchors. All four tutorial pages pass figure/table caption and reference
  checks; the partition page retains four figures and one table, and the
  estimation page has four tables. All page version badges are 0.4.0.
  `git diff --check` passes. The new page is present in article navigation,
  metadata, search, and sitemap.
- Runtime code is unchanged; the preceding 856-expectation suite remains the
  latest runtime verification and was not repeated for this documentation edit.

No commit, push, tag, or release.

## Follow-up: define sample-count notation before use

The posterior-propriety paragraph previously used the expected sample-count
vector before its definition in the trajectory update. Moved the condition
after the Gaussian precision is introduced and expanded the first definition:
n_mk is the sum of sample-position probabilities for bin k of ordering m;
the vector n_m has K entries, and their sum is the total sample count n.
These are soft sample counts, separate from the expected feature counts used
in dimension reporting. The canonical article and search index regenerate;
the site audit passes for 32 pages and 1120 links/assets/anchors. Diff whitespace
checks pass after rendering completes. No fitting code changes or commit/push.


## Follow-up: concise initialization guidance

Reduced the introductory initialization section to two short paragraphs:
feature similarity and default single linkage form M initial groups, with
one independently initialized sample ordering per group; the four metric
names, a constructor setting, and a link to exported function help follow.
Specified-M fitting defaults to Spearman similarity and each group's own
PC1. Removed the repeated initialization description from the fitting
example and moved alternative initialization and metric definitions to the
roxygen help for mpcurve_init_control(). The helper now documents absolute
correlations, directional spline R-squared, smooth-fit evidence scores,
and the distinct automatic-initialization default.

Regenerated the constructor Rd/reference page, partition article, and search
index. Browser inspection confirms the concise section and four metric
entries, zero KaTeX errors, and no width overflow at 1053 pixels. The site
audit passes 32 HTML pages and 1121 local links/assets/anchors, including the
new metric-help anchor; numbered figure/table checks and diff whitespace
checks pass. No fitting behavior changes; runtime tests were not repeated
for this documentation edit. No commit or push.


## Follow-up: reuse defined parameter notation

The Noise and smoothness settings now identify the estimated per-feature
noise variance as sigma_j squared and its initialization with the same
notation. Smoothness initialization and the fixed-at-five example likewise
use the already defined lambda_jm. Regenerated the partition article and
search index. Actual browser inspection confirms correct math rendering,
zero KaTeX errors, and no width overflow at 1053 pixels. The site audit passes
32 pages and 1121 local links/assets/anchors; diff whitespace checks pass.
No fitting code changes or commit/push.


## Follow-up: simplify known measurement-error notation

Removed the auxiliary variance symbol v_ij from the known-error formulas.
The text explicitly defines S_ij as the supplied measurement standard
deviation for sample i and feature j, and uses S_ij squared directly in
the likelihood and trajectory updates. A feature-specific vector supplies
S_ij = S_j for all samples. Regenerated the article and search index;
expanded browser inspection confirms correct equations, zero KaTeX errors,
and no width overflow. The audit passes 32 pages and 1121 local links/assets/
anchors, with figure/table and diff whitespace checks. No fitting behavior
changes or commit/push.


## Follow-up: remove advanced tuning from the introduction

Removed the paragraph covering RW3, initial/fixed smoothness precision,
parameter bounds, and the optional smoothness prior from the model tutorial.
These options and examples remain in exported mpcurve_control() help; the
introduction retains a single link to advanced fitting settings.
Regenerated the article and search index. Actual browser inspection confirms
the paragraph is absent, the helper link remains, and there are no math
errors or width overflow. The audit passes 32 pages and 1121 local links/
assets/anchors, with figure/table and diff whitespace checks. No runtime
behavior changes or commit/push.


## Follow-up: one initialization per dimension-selection candidate

Corrected the selector's legacy automatic comparison of two initialization
strategies. It now honors the supplied initial_method and partition_init,
fitting each candidate dimension once. The estimation article and its tables
regenerate from the corrected code; the example still selects M = 2, with
its default-strategy objective now 2353.147. Removed the redundant retained-
fit column from Table 4. See 2026-10-01-dimension-selection-single-initialization.md
for the correction, current 86-expectation checks, and updated results.


## Follow-up: return the estimated model itself

Automatic adaptive fits now return only the supported orderings, with no
additional CAVI sweeps. The estimation article's Table 1 reports the single
intrinsic_dim of that returned object; initialization and removal history
are available in dimension_estimation. Explicit integer dimensions are kept.
README, references, all four articles, news, search, and sitemap regenerate.
See 2026-10-01-automatic-dimension-finalization.md for runtime details and the
current 1028-expectation full suite. Browser checks pass; the site audit
passes 32 pages and 1122 links/assets/anchors. No commit or push.
