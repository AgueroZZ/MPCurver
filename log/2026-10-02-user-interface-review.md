# User-facing interface review

This follow-up reviews the current local 0.4.0 interface after R1 through R11
were closed. It separates observed behavior from proposed design changes.
The initial review made no runtime changes. The user subsequently approved
U1 through U5; the implementation status and verification are recorded below.

The main fitting signature is reasonably compact: twelve arguments, with model,
initialization, and continuation settings separated into constructors. The
highest-value remaining work concerns result extraction, identity preservation,
and the first-use workflow rather than another broad argument redesign.

## U1 — Stable result extraction across ordering counts (high priority)

**Observed:** single-ordering positions use
`fit$locations$mean$pseudotime`; multiple orderings use
`fit$locations[[m]]$mean$pseudotime`. Trajectory means are respectively a matrix
and a list of matrices. Full trajectory posteriors also use different access
paths. Automatic dimension estimation can return either schema. Summary lists
likewise place diagnostics under different fields. The public object still
exposes dimension aliases and internal views such as `model_intrinsic_dim`,
`requested_intrinsic_dim`, `active_intrinsic_dim`, `displayed_intrinsic_dim`,
`fit`, and, for multiple orderings, `fits`.

**Implication:** downstream code must branch on the fitted dimension to perform
ordinary extraction; users must distinguish public results from internal and
compatibility fields. This behavior is documented, but remains a usability cost.

**Proposal:** establish a small stable result-extraction contract, potentially
with accessors that return sample-by-ordering positions and feature-by-ordering
probabilities even for one ordering. Keep detailed state available for advanced
use without making it the ordinary extraction interface. A full internal object
rewrite is not necessary to achieve this.

## U2 — Preserve sample/feature identities and support named plotting (high priority)

**Observed:** with explicitly named input rows and columns, posterior-mean
pseudotime vectors have no sample names, trajectory matrices have no feature
row names, and multi-ordering assignment probabilities have no feature row
names. The original names remain in `fit$data`, so correspondence can be
reconstructed by position; no data reordering defect is established here.
`plot(fit, dims = "gene_1")` emits a coercion warning followed by
`missing value where TRUE/FALSE needed`. Default axis labels use dimension
numbers rather than available feature names. For a valid one-feature fit,
`plot(fit)` also requires an explicit `dims = 1` because the default requests
two features.

**Proposal:** carry input identities into ordinary outputs, accept feature names
in `dims`, use feature names as default labels, and make default plotted
dimensions adapt to the number of available features. This reduces manual
index matching when joining results to annotations.

## U3 — Validate the data contract at the public boundary (medium priority)

**Observed:** a matrix with one missing value fails with a message framed as
`PCA initialization failed: infinite or missing values in 'x'`; a matrix
containing text fails as `PCA initialization failed: 'x' must be numeric`.
These are unsupported inputs, not evidence of a numerical failure that a
different initializer should repair. The messages omit the offending named
sample/feature.

**Proposal:** validate numeric and finite input before initialization and state
the input requirement directly, ideally identifying the first offending
sample/feature. Add a short data-preparation contract covering orientation,
missing values, and the units/alignment of supplied measurement SDs. This does
not imply silently imputing, filtering, or standardizing data.

## U4 — Make the two dimension workflows easy to choose (medium priority)

**Observed:** `intrinsic_dim = "auto"` uses similarity-based initialization
followed by adaptive fitting and numerical-zero removal. In contrast,
`select_mpcurve_dimension()` compares candidate fits under a fixed uniform
partition prior. Current detailed documentation explains these distinctions;
they are not interchangeable names for one procedure.

**Proposal:** put a short choice guide near the first multi-ordering example:
specified dimension, automatic adaptive dimension, and fixed-uniform candidate
comparison. Retain the existing statistical distinction and avoid making users
infer it from function names or multiple reference pages.

## U5 — Put a short executable workflow before the derivation (medium priority)

**Observed:** the getting-started vignette begins with the model and extensive
CAVI derivations; its first actual fitting call appears around line 293. The
README's opening worked example uses 1,500 samples and follows the initial fit
with three additional initialization comparisons.

**Proposal:** lead with a compact simulate/load, fit, inspect convergence, plot,
and extract workflow; place the derivation and initialization comparison after
that first successful run. This is an information-ordering change, not a
request to remove the scientific explanation.

## Evidence and scope

Read the exported namespace, README source, public reference documentation,
fitting/control interfaces, result conversion and printing, plotting, and the
dimension selector. Ran small public fits with one, two, and automatic
orderings on a named 40-by-5 matrix (seed 614, five bins, one sweep), named-feature
plotting, missing/non-numeric input probes, and a one-feature default plot.
Evidence: `/tmp/mpcurver-user-interface-review.R`, `.log`, and
`/tmp/mpcurver-user-interface-extra.R`, `.log`.

These are interface observations and design recommendations, not a new claim
that valid numerical fitting is incorrect. The preceding package check remains
the latest full check; this review did not rerun it. Existing convergence status
printing and explicit initialization-fallback warnings are useful and do not
need another warning policy change on the evidence reviewed here.

Suggested order: U1 and U2 first, U3 next, then U4/U5 documentation organization.
At the initial review stage there were no implementation edits, commits, or pushes.


## Implementation after approval

U1 through U5 are implemented in the current working tree:

- **U1:** `fitted_positions()`, `fitted_trajectories()`, and
  `fitted_assignments()` provide consistent arrays/matrices across single,
  multiple, and automatic ordering counts. They expose means, uncertainty,
  full position probabilities, and trajectory covariance without backend
  branching. Trajectories for multiple orderings remain conditional on feature
  assignment. Internal state fields remain available for detailed inspection.
- **U2:** Public positions, trajectory matrices, assignment probabilities, and
  variance/precision outputs preserve input identities. Plotting accepts feature
  names, uses them on axes, and defaults to one dimension for one-feature data.
  Unknown/ambiguous names and invalid indices give explicit errors. Automatic
  pruning retains the surviving ordering label across result accessors and
  named `fitted_prior()` calls, including after continuation.
- **U3:** The public fit and dimension selector validate numeric, real-valued,
  finite sample-by-feature data before initialization. All-numeric data frames
  are accepted. Nonfinite data errors identify the first offending sample and
  feature, with their names when available. Measurement SDs are checked for
  supported dimensions, finite nonnegative values, and numeric type. Documentation
  states positional alignment, units, zero-SD flooring, and that data are not
  imputed or standardized by fitting.
- **U4:** The README, multi-ordering tutorial, and dimension tutorial recommend
  adaptive automatic fitting first when the dimension is unknown because it
  fits one adaptive model rather than a sequence of candidate dimensions.
  They distinguish specified dimensions and fixed-uniform candidate comparisons,
  whose prior and selection rule differ from adaptive estimation.
- **U5:** A small executable README workflow precedes the worked spiral.
  `mpcurve_intro.Rmd` and `partition.Rmd` present models and software examples;
  `cavi_single.Rmd` and `cavi_partition.Rmd` retain the algorithm details on
  separate, linked pages, including known-measurement-error updates. The site
  has a separate Algorithm details navigation group. Existing tutorial URLs
  and historical analyses are retained.

Current implementation evidence: `R/15_results.R`, the public fitting and plot
entry points, `tests/testthat/test-result-interface.R`, the vignette sources,
README source, `_pkgdown.yml`, and `scripts/build_public_site.R`.

Current source tests pass **1,368 expectations**, with zero failures, errors,
unexpected warnings, or skips. Checks cover actual single/multiple/automatic
fits, numerical-zero pruning to one ordering, continuation, covariance and
position-probability extraction, input validation, and named/one-feature plots.
The dedicated named-plot image has been inspected. Full site generation and
README execution pass. All 38 HTML pages and 1,401 local links/assets/anchors
pass the site audit, including figure/table numbering on all published articles.
Fifteen browser-checked pages have no math-rendering errors, unloaded images,
or page-width overflow. Final heading/wording edits were rendered again and
rechecked on the two affected pages; the result table and algorithm page
screenshots were inspected.

`R CMD check --no-manual` completes with **Status: OK**, with no check errors,
warnings, or notes. Installed tests pass 1,350 expectations with one standard
CRAN snapshot skip. Examples and all six vignette code runs and rebuilds pass.
All 22 runtime R sources and 44 Rd files match the checked source archive byte
for byte. The final plain-text algorithm headings and dimension-description
wording were edited after archive creation and verified through targeted site
rebuilds and browser checks; no runtime source changed after that archive.

Temporary evidence is under `/tmp/mpcurver-ui-*`; the final runtime archive
and package checks use `/tmp/mpcurver-ui-final-check/`. The earlier
`/tmp/mpcurver-ui-check/` run predates the final input and ordering-label checks
and is historical. No commit or push.
