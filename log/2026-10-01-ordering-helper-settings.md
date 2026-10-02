# Ordering-helper settings — 2026-10-01

Clarified that arguments for the selected initialization method are passed
through `mpcurve_init_control(method_args = list(...))`. Both Fiedler and
Isomap use `num_neighbors` for their nearest-neighbor graphs. Added complete
fitting examples for both methods and links to their exported helper
documentation in the fitting and initialization-control help.

Runtime interfaces and algorithms are unchanged. Current targeted checks
execute both documented single-ordering examples and both methods with two
feature groups. All four fits accept `num_neighbors = 20`; the grouped fits
use the requested methods within their feature groups. The preceding full
suite result belongs to the similarity-only cleanup entry.

Regenerated the two reference pages and index, search, and sitemap. Browser
checks confirm the example text and helper links, no math errors, and no
horizontal page overflow. The site audit passes 33 pages and 1130 local links,
assets, and anchors; existing caption numbering and whitespace checks pass.
No commit or push.
