
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MPCurver

<!-- badges: start -->

<!-- badges: end -->

`MPCurver` fits smooth latent orderings with Gaussian-mixture / GMRF
models.

The public workflow is **`mpcurve`-first and CAVI-first**: use
`fit_mpcurve()` to fit a model, `do_mpcurve()` to continue a fitted
model, and the usual `print()`, `summary()`, and `plot()` methods to
inspect results.

## Documentation

The package documentation is organized around the public `mpcurve`
workflow.

- `vignettes/mpcurve_intro.Rmd`: current single-ordering introduction.
- `vignettes/partition.Rmd`: current multi-ordering / feature-partition
  tools.
- `vignettes/fitness.rmd`: applied case study using the public
  interface.

## Installation

You can install the development version of `MPCurver` like so:

``` r
pak::pak("AgueroZZ/MPCurver")
```

## Example

This is a basic example using the recommended `mpcurve` interface:

``` r
library(MPCurver)
```

We will simulate a moderately easy 2D spiral and fit the default
MPCurver workflow.

``` r
sim <- simulate_spiral2d(n = 1500, turns = 1, noise = 0.08, seed = 123)
plot(sim$obs, pch = 16, cex = 0.35, col = "grey80",
     xlab = "x1", ylab = "x2", main = "Noisy 2D spiral")
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

Now we fit the model with `fit_mpcurve()`. The public wrapper uses the
CAVI backend and returns an `mpcurve` object.

``` r
set.seed(123)
fit <- fit_mpcurve(
  X = as.matrix(sim$obs),
  method = "isomap",
  K = 40,
  rw_q = 2,
  iter = 150,
  discretization = "equal",
  tol = 1e-8
)

plot(fit)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

``` r
plot(fit, plot_type = "elbo")
```

<img src="man/figures/README-unnamed-chunk-3-2.png" width="100%" />

We can compare several initialization strategies under the same `cavi`
model. The explicit serial loop below is reproducible on platforms where
forked parallel numerical libraries are unavailable:

``` r
init_methods <- c("PCA", "fiedler", "isomap")
fits <- setNames(
  lapply(init_methods, function(init_method) {
    fit_mpcurve(
      X = as.matrix(sim$obs),
      method = init_method,
      K = 40,
      rw_q = 2,
      iter = 150,
      discretization = "equal",
      tol = 1e-8
    )
  }),
  init_methods
)

data.frame(
  method = names(fits),
  elbo_last = unname(vapply(fits, function(x) tail(x$elbo_trace, 1), numeric(1)))
)
#>    method elbo_last
#> 1     PCA -424.1734
#> 2 fiedler -375.5500
#> 3  isomap -376.9430
```

Each fit keeps the unified `params`, `gamma`, `elbo_trace`,
`loglik_trace`, and the underlying backend fit in `$fit`.

For multi-ordering analyses, use `fit_mpcurve(..., intrinsic_dim = 2)`
or more generally `fit_mpcurve(..., intrinsic_dim = M)` with a fixed,
user-chosen `M >= 2`.

## Fixed-M structural partition models

Multi-ordering fits use the structural variational family

$$
q(C)\prod_j q(Z_j)q(U_j \mid Z_j),
$$

so each trajectory posterior is conditional on the feature assignment
rather than treated as an independent potential trajectory. Estimated
feature variances are represented by one shared length-$d$
`params$sigma2` vector across all orderings. The fitted `mpcurve` object
stores one canonical structural state, including the named
ordering-specific `gamma` matrices, `conditional_posterior`,
`lambda_mat`, and `partition$pi_weights`; `$fits` is a derived
compatibility view, not an independent fitting state.

There are no active/inactive, freeze, or drop states in this model.
Automatic `greedy = "forward"` and `greedy = "backward"` selection is
temporarily unavailable because the raw fixed-$M$ structural ELBO is not
a complexity-penalized criterion for comparing different values of $M$.
Choose `M` from scientific context and assess alternatives through
stability or external validation.

The same `print()`, `summary()`, and `plot()` interface works for both
single-ordering and structural partition fits. For a partition fit,
`plot(fit, plot_type = "mu", dims = j)` displays the conditional
posterior mean $E[U_j \mid Z_j=m]$ separately for every ordering $m$.
