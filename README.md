# Khan

`Khan` is the R companion package for
**[The Wreaths of KHAN: Uniform Graph Feature Selection with False Discovery
Rate Control](https://arxiv.org/abs/2403.12284)** (Liang, Zhang, and
Neykov, 2024). It provides FDR-controlled inference on subgraph features
and persistent homology generators for Gaussian graphical models and
ferromagnetic Ising models.

## Features

- **BHq subgraph selection (Algorithm 1)** — `bhq_test()`,
  `bonferroni_test()` for any user-supplied feature list (triangles,
  cliques, trees, paths, ...).
- **Discrete Gram-Schmidt (Algorithm 2)** — `dgs()` selects linearly
  independent homology generators at a fixed filtration level. The
  `k = 1` path uses Union-Find for near-linear-time rank tracking.
- **KHAN (Algorithm 3)** — `khan()` adaptively walks the filtration line
  and applies DGS at each change point to produce a uFDR-controlled
  persistent homology barcode.
- **Debiased graphical lasso** — `debias_glasso()`,
  `estimate_variance_ggm()`, `compute_ggm_pvalues()` implement the
  estimator from equations (2.7)-(2.8).
- **Ising model support** — `generate_ising_trees()`,
  `estimate_ising_weights()`, `compute_ising_pvalues()`,
  `verify_ising_correlation()`.
- **Bootstrap baselines** — `bootstrap_feature_test()` and
  `bootstrap_simultaneous_test()` for comparison against the Gaussian
  multiplier max-statistic approach.
- **Subgraph enumeration** — `find_all_features()`, `find_cycles()`,
  `find_cliques_by_size()`, `find_trees_5node()` plus edge-extraction
  helpers.
- **Real-data wrappers** — `khan_realdata()` and `bhq_realdata()` chain
  estimation, debiasing, p-values, and selection in a single call.

## Installation

```r
# install.packages("remotes")
remotes::install_github("h7nian/Khan")
```

The package depends on `igraph`, `MASS`, `CVglasso`, and `IsingSampler`.

## Quick start

### Subgraph selection on a Gaussian graphical model

```r
library(Khan)

set.seed(1)
sim   <- generate_ggm_subgraph(m1 = 5, m2 = 3, m3 = 2)
data  <- generate_ggm_data(n = 400, omega = sim$omega)

emap     <- idx_map(sim$d)
fit      <- estimate_ggm(data$X)
theta_d  <- debias_glasso(fit$sigma_hat, fit$theta_hat)
var_mat  <- estimate_variance_ggm(theta_d)
p_edges  <- compute_ggm_pvalues(theta_d, var_mat, n = 400,
                                scenario = "a", emap = emap)

features   <- find_all_features(sim$theta, emap,
                                feature_types = c("triangle", "four_cycle",
                                                  "five_cycle"))
candidates <- unlist(features, recursive = FALSE)

result <- bhq_test(p_edges, candidates, q = 0.1)
str(result)

# Or the same pipeline in a single call:
bhq_realdata(X = data$X, feature_edges = candidates, q = 0.1)
```

### KHAN persistent homology

```r
sim  <- generate_ggm_homology(m1 = 4, m2 = 6, m3 = 4)
data <- generate_ggm_data(n = 500, omega = sim$omega)

khan_result <- khan_realdata(X = data$X, q = 0.1,
                             mu_range = c(0, 2), scenario = "a")
khan_result$barcode
khan_result$change_points
```

## Package layout

```
R/
  edge_index.R      Edge-index map and reverse lookup
  graph_features.R  Cycle / clique / tree / path enumeration
  homology.R        H_k rank via Union-Find (k=1) or QR (k=2)
  ggm.R             GGM data generation, GLasso, debiasing, p-values
  ising.R           Ising data generation, estimation, p-values
  bhq.R             Algorithm 1 (BHq) and Bonferroni
  dgs.R             Algorithm 2 (Discrete Gram-Schmidt)
  khan.R            Algorithm 3 (KHAN persistent homology)
  bootstrap.R       Gaussian multiplier bootstrap baselines
  metrics.R         FDP, Power, uFDP evaluation
  khan_realdata.R   End-to-end wrappers for real datasets
tests/testthat/     Regression tests for every algorithm
```

## Testing

```r
devtools::test()      # or: R CMD check Khan_*.tar.gz
```

The `tests/testthat/` directory contains regression tests covering BH
thresholding edge cases, the debiased estimator against equation (2.7),
variance scaling, DGS rank tracking versus a brute-force reference,
homology rank for small graphs, and edge-index consistency.

## Citation

```bibtex
@article{liang2024wreaths,
  title   = {The Wreaths of KHAN: Uniform Graph Feature Selection
             with False Discovery Rate Control},
  author  = {Liang, Jiaqi and Zhang, Sinian and Neykov, Matey},
  journal = {arXiv preprint arXiv:2403.12284},
  year    = {2024}
}
```

## License

MIT (see [LICENSE.md](LICENSE.md)).
