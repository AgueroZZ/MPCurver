
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MPCurver

<!-- badges: start -->

<!-- badges: end -->

`MPCurver` estimates latent sample orderings and smooth feature
trajectories from multivariate data. It represents an ordering by a grid
of positions and models each feature as a smooth function along that
grid. Gaussian measurement models and random-walk priors allow it to
estimate trajectories together with their uncertainty.

For data with several sources of variation, MPCurver can learn multiple
orderings and the probability that each feature follows each ordering.
This allows different groups of features to describe different patterns
among the same samples.

Use `fit_mpcurve()` to fit a model, `summary()` and `plot()` to explore
the results, and `do_mpcurve()` to run additional fitting iterations.

## Documentation

- [Getting
  started](https://aguerozz.github.io/MPCurver/articles/mpcurve_intro.html):
  the model, variational inference, and a simulated trajectory example.
- [Feature
  partitioning](https://aguerozz.github.io/MPCurver/articles/partition.html):
  fitting multiple orderings and interpreting feature assignments.
- [Fitness data
  analysis](https://aguerozz.github.io/MPCurver/articles/fitness.html):
  a worked example with smoothing and measurement uncertainty.

## Installation

Install `MPCurver` from GitHub:

``` r
pak::pak("AgueroZZ/MPCurver")
```

## Example

The example below recovers a latent ordering from a noisy
two-dimensional spiral.

``` r
library(MPCurver)
```

First, generate the observations:

``` r
sim <- simulate_spiral2d(n = 1500, turns = 1, noise = 0.08, seed = 123)
plot(sim$obs, pch = 16, cex = 0.35, col = "grey80",
     xlab = "x1", ylab = "x2", main = "Noisy 2D spiral")
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

Fit a model with 40 grid positions and a second-order random-walk prior.
The first plot shows the inferred ordering and trajectory; the second
shows the variational objective over fitting iterations.

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

Initialization can affect the fitted ordering, especially for curved or
weakly separated patterns. Here we compare three starting orderings
under the same model and fitting settings:

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

For these fits with the same model and data, a larger final ELBO
indicates a better variational objective. Inspect the fitted
trajectories as well as this numerical comparison. Sample positions are
available in `fit$locations$mean$pseudotime`, and `fit$params$mu`
contains the estimated feature trajectories.

## Multiple orderings and feature groups

A single ordering may describe some features well while missing the
patterns in others. Set `intrinsic_dim = M` to fit $M$ orderings of the
same samples. MPCurver jointly estimates those orderings, smooth
trajectories, and a probability distribution over orderings for each
feature. For example, two groups of genes may vary along different
cell-state axes.

Inference uses structural variational inference: a feature’s trajectory
distribution is conditional on its ordering assignment. This preserves
the dependence between which ordering a feature follows and the shape of
its trajectory. Each feature has one noise variance shared across the
candidate orderings, while smoothness can vary with the feature and
ordering.

For a two-ordering analysis of your sample-by-feature data matrix `X`:

``` r
fit_multi <- fit_mpcurve(X, intrinsic_dim = 2, K = 40)
summary(fit_multi)
fit_multi$partition$pi_weights
plot(fit_multi, plot_type = "mu", dims = 1)
```

The assignment probabilities describe how strongly each feature supports
each ordering. The trajectory plot shows the feature’s fitted mean under
each ordering, with the corresponding assignment probability in the
panel title. Choose the number of orderings for the scientific question
and examine the stability and interpretability of the resulting feature
groups. The [feature-partitioning
tutorial](https://aguerozz.github.io/MPCurver/articles/partition.html)
develops this workflow and the statistical model in detail.

## Measurement uncertainty

By default, MPCurver estimates a noise variance for each feature. When
measurement standard deviations are available, supply them through `S`
as either a vector of feature-specific values or a matrix matching `X`.
The likelihood then accounts for the supplied precision of each
measurement. See the [fitness case
study](https://aguerozz.github.io/MPCurver/articles/fitness.html) for an
example.
