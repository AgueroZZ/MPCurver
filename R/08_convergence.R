# Scale an ELBO increment without changing the stored objective. Under the
# relative rule each backend supplies its own offset for the denominator.
.mpcurve_elbo_change <- function(delta, previous, n, d, convergence,
                                 relative_offset, n_observed = NULL) {
  denominator <- if (identical(convergence, "normalized")) {
    n_observed %||% (as.double(n) * as.double(d))
  } else {
    abs(previous) + relative_offset
  }
  delta / denominator
}
