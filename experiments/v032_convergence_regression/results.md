# MPCurver 0.3.2 simulation regression results

All 47 fresh fits passed the predeclared regression checks. None of the fitting calls emitted warnings.

## Earlier examples

- M = 1 Swiss roll: ordering Spearman correlation 0.999758, unchanged.
- M = 2 partition: ARI 1; mean ordering correlation 0.998205, unchanged at the original iteration budget.
- All 27 ND-stopping fits exactly reproduced the experimental 0.3.1 ELBO traces and stopping indices; fitted parameter comparisons passed at tolerance 1e-7.

## Intrinsic-M pilots

| True M | Method | Old estimate | New estimate | New ARI | New ordering recovery | Recovery change |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 3 | Adaptive EB | 4 | 4 | 0.975200 | 0.993395 | +0.000030 |
| 3 | Uniform + forward | 3 | 3 | 1.000000 | 0.993600 | -0.000001 |
| 4 | Adaptive EB | 4 | 4 | 1.000000 | 0.989587 | +0.000004 |
| 4 | Uniform + forward | 4 | 4 | 1.000000 | 0.989587 | -0.000012 |
| 5 | Adaptive EB | 5 | 5 | 1.000000 | 0.989029 | +0.000090 |
| 5 | Uniform + forward | 5 | 5 | 1.000000 | 0.989014 | +0.000078 |

All 18 intrinsic-M candidate fits converged at the first eligible normalized threshold crossing predicted from the saved old traces. Their common trace segments had zero numerical difference. The known adaptive split for true M = 3 persists; it was already present before this update.

Fresh fitting time: 1027.7 seconds (one process, one CPU thread).

The original M = 2 website example reaches its fixed iteration budget in both versions. Its regression check concerns matched-budget recovery; the intrinsic-M candidates and all ND-stopping fits were required to converge.
