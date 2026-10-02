# Simulation model naming — 2026-10-01

Renamed the public `simulate_cavi_toy()` helper to `simulate_mpcurve()` and
its private numerical implementation to `.simulate_mpcurve`. The reference
title now identifies a single-ordering MPCurve model. CAVI refers to the
fitting algorithm, including the corrected internal fitting title.

Updated the fitting example, namespace, reference configuration, curated site
builder, and active tests and manual examples. The manual examples now load
the complete checkout through `pkgload::load_all()` so their CAVI calls have
all current helpers available. Historical experiment outputs retain their
original names and provenance. The 0.4.0 changelog records the public rename;
the retired function name is not exported.

Generating settings and returned fields are unchanged. Corrected help text
to identify the existing precision and noise-SD draws as log-uniform. With
matching seeds, default and explicitly controlled draws are identical to
results saved before the rename.

Regenerated the roxygen namespace and manuals, the simulation and fitting
reference pages, reference index, news, search, and sitemap. The site builder
preserves the old reference URL with a redirect to the new page. This update
was rendered in place to keep the running local preview available.

Current verification:

- Complete suite: 1028 passing expectations; no failures, errors, warnings,
  or skips.
- Both affected manual examples run successfully, including the five-seed
  larger-grid example.
- Browser inspection confirms the old-URL redirect, corrected heading and
  usage, no math errors, and no horizontal overflow at a 1053-pixel viewport.
- Site audit: 33 HTML pages and 1125 local links, assets, and anchors checked
  without issues; existing article figure and table numbering passes.
- Whitespace checks pass. No commit or push.
