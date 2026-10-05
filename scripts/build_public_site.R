#!/usr/bin/env Rscript

override <- list(home = list(sidebar = FALSE))

if (!isTRUE(l10n_info()[["UTF-8"]])) {
  stop("The public site must be built in a UTF-8 locale.")
}

public_reference_topics <- c(
  "MPCurver",
  "mpcurve",
  "fit_mpcurve",
  "mpcurve_control",
  "mpcurve_init_control",
  "mpcurve_continue_control",
  "select_mpcurve_dimension",
  "do_mpcurve",
  "fitted_prior",
  "fitted_positions",
  "fitted_trajectories",
  "fitted_assignments",
  "print.mpcurve",
  "summary.mpcurve",
  "plot.mpcurve",
  "PCA_ordering",
  "fiedler_ordering",
  "isomap_ordering",
  "pcurve_ordering",
  "tSNE_ordering",
  "simulate_mpcurve",
  "simulate_dual_trajectory",
  "simulate_intrinsic_trajectories",
  "simulate_spiral2d",
  "simulate_swiss_roll_1d_2d",
  "simulate_two_order_gp_dataset"
)

public_articles <- c(
  "mpcurve_intro",
  "partition",
  "intrinsic_dimension",
  "fitness",
  "cavi_single",
  "cavi_partition"
)

# Execute vignette code against the current checkout rather than a potentially
# stale installed copy of MPCurver. This mirrors the source-loading step used by
# the ordinary pkgdown build while preserving the curated public-site surface.
pkgload::load_all(
  ".",
  quiet = TRUE,
  export_all = FALSE,
  helpers = FALSE,
  attach_testthat = FALSE,
  warn_conflicts = FALSE
)

rewrite_pkgdown_metadata <- function(path, public_articles) {
  if (!file.exists(path)) {
    return(invisible(NULL))
  }

  lines <- readLines(path, warn = FALSE)
  output <- character()
  in_articles <- FALSE

  for (line in lines) {
    if (identical(line, "articles:")) {
      in_articles <- TRUE
      output <- c(output, line)
      next
    }

    if (in_articles && grepl("^[^[:space:]]", line)) {
      in_articles <- FALSE
    }

    if (in_articles) {
      match <- regmatches(
        line,
        regexec("^  ([^:]+): [^[:space:]]+\\.html$", line)
      )[[1]]

      if (length(match) > 1L) {
        if (match[2] %in% public_articles) {
          output <- c(output, line)
        }
        next
      }
    }

    output <- c(output, line)
  }

  writeLines(output, path)
}

if (dir.exists("docs")) {
  unlink("docs", recursive = TRUE, force = TRUE)
}

dir.create("docs", recursive = TRUE, showWarnings = FALSE)

pkg_home <- pkgdown:::section_init(".", override = override)
pkgdown:::build_home_index(pkg_home, quiet = TRUE)
if (!pkg_home$development$in_dev) {
  pkgdown:::build_404(pkg_home)
}

pkg_ref <- pkgdown:::section_init(".", "reference", override = override)
pkg_ref$topics <- pkg_ref$topics[
  pkg_ref$topics$name %in% public_reference_topics,
  ,
  drop = FALSE
]
pkgdown:::build_reference_index(pkg_ref)
pkgdown::build_reference(
  pkg = pkg_ref,
  lazy = FALSE,
  examples = FALSE,
  preview = FALSE,
  devel = FALSE
)

pkg_art <- pkgdown:::section_init(".", "articles", override = override)
pkg_art$vignettes <- pkg_art$vignettes[
  pkg_art$vignettes$name %in% public_articles,
  ,
  drop = FALSE
]
pkgdown:::build_articles_index(pkg_art)
for (name in public_articles) {
  pkgdown:::build_article(
    name,
    pkg = pkg_art,
    lazy = FALSE,
    seed = 1014L,
    new_process = FALSE,
    quiet = TRUE
  )
}

fitness_html <- readLines("docs/articles/fitness.html", warn = FALSE)
if (any(grepl("<ce>|<bc>", fitness_html, fixed = FALSE))) {
  stop("The fitness article contains malformed UTF-8 markup.")
}

pkgdown::build_news(pkg_home, preview = FALSE)
pkgdown:::build_search(pkg_home)
pkgdown:::build_sitemap(pkg_home)
rewrite_pkgdown_metadata("docs/pkgdown.yml", public_articles)

# Preserve bookmarked reference URLs after the public simulation helper rename.
write_public_reference_redirects <- function() {
  writeLines(c(
    '<!doctype html>',
    '<html lang="en">',
    '<head>',
    '  <meta charset="utf-8">',
    '  <title>simulate_mpcurve — MPCurver</title>',
    '  <meta http-equiv="refresh" content="0; url=simulate_mpcurve.html">',
    '  <link rel="canonical" href="simulate_mpcurve.html">',
    '</head>',
    '<body>',
    '  <p>This reference has moved to <a href="simulate_mpcurve.html">simulate_mpcurve()</a>.</p>',
    '</body>',
    '</html>'
  ), "docs/reference/simulate_cavi_toy.html")
}
write_public_reference_redirects()

# Keep generated HTML free of trailing blank lines across pkgdown builds.
normalize_public_html <- function() {
  paths <- list.files("docs", pattern = "\\.html$", recursive = TRUE,
                      full.names = TRUE)
  for (path in paths) {
    lines <- readLines(path, warn = FALSE)
    lines <- sub("[[:blank:]]+$", "", lines)
    if (length(lines) && !nzchar(tail(lines, 1L))) {
      while (length(lines) && !nzchar(tail(lines, 1L))) {
        lines <- head(lines, -1L)
      }
    }
    writeLines(lines, path)
  }
}
normalize_public_html()
