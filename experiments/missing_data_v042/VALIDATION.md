# MPCurver 0.4.2 validation record

Validation dates: 2026-10-18 and 2026-10-19. Development branch:
`codex/v0.4.2-missing-data`. Baseline main commit:
`4f8b583d7549d4068d0e810e93f90653cbcb2ed9`.

## Numerical and interface checks

- Full source test suite: 1,885 assertions passed, zero failures, errors,
  warnings, or skips.
- Twenty-four complete-data configurations: exact equality to the frozen
  0.4.1 installation across the compared fitted outputs and histories.
- Seventy-two missing-pattern fits: all independent scalar-objective and
  monotonicity checks passed. Maximum objective discrepancy was `1.023e-12`.
- Observed-entry kernels use masked matrix products instead of per-feature
  loops. Against the loop implementation they replaced, kernel outputs agreed
  to at most `5e-14` relative. Eighteen public fitting configurations (one to
  three orderings, estimated noise and known SDs, annealing, continuation, and
  dimension selection) had identical iteration counts. Position probabilities
  agreed to `1.5e-11`, trajectory means to `8e-10`, and objectives to `6e-11`
  relative. With n = 400, K = 20, 5% missing entries, and 15 sweeps,
  single-ordering fits took 0.8/1.2/3.4 s at d = 100/400/1000, versus
  1.7/4.2/10.5 s before.
- Singular intrinsic systems are rejected in every fit, with or without
  missing values, when the bins carrying data do not span the prior null
  space. Numerically unreliable missing-data systems are rejected when the
  reciprocal 1-norm condition number is at most `1e-12`. The
  single-observation intrinsic stress case stops with the singular error; its
  explicitly proper-prior counterpart passes the objective check.
- Precision initialization uses only data-carrying bins, and feature
  similarities use pairwise-complete observations. After these changes the 24
  complete-data configurations above remain identical to 0.4.1.
- On 2026-10-19 the unused SmoothEM and csmoothEM code was removed, the two
  ordering-initialization helpers were renamed, and random starts were routed
  through the common ordering initialization. Before/after snapshots of 43
  non-random configurations (one and two orderings, automatic dimension,
  continuation, dimension selection, known noise, missing data, and the
  internal initialization helpers) were identical, including warnings. The 24
  complete-data configurations remain identical to 0.4.1, and the
  missing-pattern and boundary CSV files are byte-identical. The same holds
  after the later removal of unreachable CAVI-era code (the augmented
  two-trajectory partition, greedy ordering-count search, and the raw-fit
  print method). None of the removed functions was called by any Rd example,
  vignette, or public fitting, print, summary, plot, or extraction path.
- Removing compatibility code for objects saved by earlier versions and for
  removed internal arguments left 40 public configurations identical in their
  fitted objects, printed output, summaries, plots, priors, and extracted
  results, apart from the deleted empty fields (`assignment_prior`,
  `ordering_alpha`, `legacy_assignment_prior`, `intrinsic_dim_semantics`,
  summary `greedy_selection` and `ordering_events`). The 43 configurations
  above differ only by the same fields, and the 24 complete-data
  configurations remain identical to 0.4.1.
- Removing the constant compatibility fields (`model_intrinsic_dim`,
  `requested_`/`active_`/`displayed_intrinsic_dim`, `active_idx`,
  `assignment_mode`, `assignment_M`, `posterior_alpha`, `kl_q_p`, and the
  write-only `assignment_posterior`) left the same 40 and 43 configurations
  identical apart from those fields, including all printed output.

See `README.md`, the three CSV files, `verify.R`, and `source-sha256.txt` for
the designs, thresholds, individual results, and exact source provenance.

## Package checks

`R CMD build` completed successfully, including all vignettes. It was run
from a clean staging directory containing the package sources, because the
ignored task-local dependency library and runtime symlinks are large.

The source tarball contains the final implementation, new tests, and generated
Rd documentation, checked byte for byte against the development source.
It contains no local validation library or Git metadata. Its SHA-256 is:

```text
73953899fb4096cab546d08356aad782c2e0e82fbc5b910d155b82d7b07398e0
```

The complete `R CMD check --no-manual` finished with **Status: OK**, zero
errors, zero warnings, and zero notes. Installation, examples, source analysis,
documentation consistency, the installed-package test suite, all six vignette
code executions, and all vignette-output rebuilds passed. The build and check
were rerun after each source change. Superseded tarballs (SHA-256
`b2735cba7b12eb21e5d58a3ab97d447a4a7297f03b772949dc6e1bad35d2ce6c` before the
kernel change, `7278ee2c573771c3c44be9a4eb790932769948a29beb930aa2cf6fb7453649ea`
before the three changes,
`99624555b27d9c285714cf5aed66e98e3f5b6087012fe7e3fab90a0fb8534ae3` before the
SmoothEM/csmoothEM removal,
`360e287a958aa159e1b56de3fb034e8ab1d9b86f4160abd895ae5ef466c26379` before the
CAVI-era removal,
`5096dfeb208f529d010c0dad11fe0a12a2092c2c7e36b6fc2433504e8ffd716a` before the
compatibility removal,
`78262e813b85df2b8570de386d8f3bced6ebda7f8f6f150d5b30f9bf8efc370d` before the
constant-field removal,
`e0e50d7ed590caf314c6ca7f60c544bf7fdfe382f33f2cd5078eb4ab390628e0` before the
probability-floor removal and ordering-helper warnings,
`6c932216a6985cec70a568e90dc02544791f1079e38fb70b7f3703ffc15ba782` before the one-feature cluster
message) are kept under ignored
`tmp/independent-review-2026-10-18/superseded/`.

## Documentation checks

- Generated Rd files with roxygen2 7.3.3, matching the baseline.
- Rendered README and the missing-data mathematical derivation from R Markdown.
  The derivation's computation section was added with the kernel change and
  re-rendered; generated Rd and the public site do not depend on it.
- Rebuilt all six curated public articles and references with the existing
  `scripts/build_public_site.R`, pkgdown 2.2.0, and Font Awesome 6.5.2 assets.
- Audited all 38 generated HTML pages: no missing local resources or links,
  and no figure/table caption numbering errors.
- Compared all 49 existing tracked PNG figure files to the baseline: identical.
- Rebuilt the site after each code removal on 2026-10-19. Relative to the
  previous site, only the `fit_mpcurve`, `select_mpcurve_dimension`,
  `do_mpcurve`, `summary.mpcurve`, and `mpcurve` reference pages, the news,
  home, and introduction pages, the search index, and build metadata changed;
  all 33 PNG files under `docs/` were byte-identical. The 38 HTML pages again passed the local-resource and
  caption audits.
- Inspected the generated tutorial and its missing-data example in a local
  headless browser, and inspected its scientific figure files. Browser preview
  copies omit external CDN dependencies; the generated site retains its
  existing KaTeX configuration.
- `git diff --check` passed.

## Environment and isolation

Checks use R 4.3.3 on CentOS Linux 8 with one OpenBLAS thread. The regression
package versions are in `session-info.txt`. Dependencies, logs, browser
previews, baseline fits, and runtime workarounds remain under ignored `tmp/`.
The machine's existing Rscript wrapper did not handle `--vanilla`; package
build/check use a task-local copy of the repaired R installation with the
native Rscript binary and its `RHOME` override. Existing runtime and library
installations were not changed. Build-time tar ownership warnings reflect the
host UID/GID encoding and do not indicate package-source failures.

The original main checkout, existing installed software, and published site
remain unchanged. No commit, push, merge, or publication has been performed.
This record verifies observed-data CAVI under ignorable missingness; it does
not establish exact posterior propriety for the intrinsic latent model.
