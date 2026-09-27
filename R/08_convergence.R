# Scale an ELBO increment without changing the stored objective. The relative
# offset is supplied by each backend to preserve its pre-0.3.2 behavior.
.mpcurve_elbo_change <- function(delta, previous, n, d, convergence,
                                 relative_offset) {
  denominator <- if (identical(convergence, "normalized")) {
    as.double(n) * as.double(d)
  } else {
    abs(previous) + relative_offset
  }
  delta / denominator
}
