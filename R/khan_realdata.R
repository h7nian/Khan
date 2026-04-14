# =============================================================================
# khan_realdata.R -- Real-data entry points for subgraph / homology selection.
#
# Thin wrappers that chain debias_glasso -> compute_ggm_pvalues -> bhq_test or
# khan(), so users can go from (sample covariance, glasso estimate) to
# FDR-controlled selections in a single call.
# =============================================================================

#' Run the full KHAN pipeline on real data.
#'
#' Given a dataset or a GLasso estimate, this function:
#' 1. estimates the precision matrix (if only `X` is provided),
#' 2. computes the debiased estimator via [debias_glasso()],
#' 3. estimates edge variances via [estimate_variance_ggm()],
#' 4. runs [khan()] for persistent-homology selection.
#'
#' @param X          Optional `n`-by-`d` data matrix.  If supplied, the
#'   sample covariance is used and [estimate_ggm()] fits the GLasso.
#' @param sigma_hat  Covariance estimate (required if `X` is `NULL`).
#' @param theta_hat  Precision estimate (required if `X` is `NULL`).
#' @param n          Sample size; inferred from `X` when possible.
#' @param q          FDR level.
#' @param mu_range   Filtration interval.
#' @param scenario   `"a"` (two-sided / GGM) or `"b"` (one-sided).
#' @param k          Homology dimension.
#' @param ...        Additional arguments forwarded to [estimate_ggm()].
#' @return A list combining the [khan()] result with the intermediate
#'   `theta_d`, `var_mat`, and `emap` so users can inspect or reuse them.
#' @export
khan_realdata <- function(X = NULL, sigma_hat = NULL, theta_hat = NULL,
                          n = NULL, q = 0.1,
                          mu_range = c(0, Inf),
                          scenario = "a", k = 1L, ...) {
  if (is.null(theta_hat) || is.null(sigma_hat)) {
    if (is.null(X)) stop("Provide either `X` or both `sigma_hat` and `theta_hat`.")
    fit       <- estimate_ggm(X, ...)
    theta_hat <- fit$theta_hat
    sigma_hat <- fit$sigma_hat
  }
  if (is.null(n)) {
    if (is.null(X)) stop("`n` must be provided when `X` is NULL.")
    n <- nrow(X)
  }

  d       <- nrow(theta_hat)
  emap    <- idx_map(d)
  theta_d <- debias_glasso(sigma_hat, theta_hat)
  var_mat <- estimate_variance_ggm(theta_d)

  # Edge-weight vector and sigma vector indexed by edge index.
  n_edges <- d * (d - 1L) / 2L
  w_hat   <- numeric(n_edges)
  sigma_e <- numeric(n_edges)
  for (u in seq_len(d - 1L)) {
    for (v in (u + 1L):d) {
      idx <- emap[u, v]
      w_hat[idx]   <- theta_d[u, v]
      sigma_e[idx] <- sqrt(max(var_mat[u, v], 1e-15))
    }
  }

  khan_result <- khan(w_hat, sigma_e, n, emap, d, q,
                      mu_range = mu_range, scenario = scenario, k = k)

  c(khan_result, list(theta_d = theta_d, var_mat = var_mat, emap = emap))
}

#' Run the BHq subgraph selection pipeline on real data.
#'
#' Computes edge p-values and applies [bhq_test()] for a user-supplied
#' candidate feature list.
#'
#' @param X              Optional data matrix.
#' @param sigma_hat      Covariance estimate (used if `X` is `NULL`).
#' @param theta_hat      Precision estimate (used if `X` is `NULL`).
#' @param n              Sample size; inferred from `X` when possible.
#' @param feature_edges  List of edge-index vectors for the candidate features.
#' @param q              FDR level.
#' @param scenario       `"a"` or `"b"`.
#' @param ...            Forwarded to [estimate_ggm()].
#' @return A list combining the [bhq_test()] result with `theta_d`, `var_mat`,
#'   `p_edges`, and `emap`.
#' @export
bhq_realdata <- function(X = NULL, sigma_hat = NULL, theta_hat = NULL,
                         n = NULL, feature_edges, q = 0.1,
                         scenario = "a", ...) {
  if (is.null(theta_hat) || is.null(sigma_hat)) {
    if (is.null(X)) stop("Provide either `X` or both `sigma_hat` and `theta_hat`.")
    fit       <- estimate_ggm(X, ...)
    theta_hat <- fit$theta_hat
    sigma_hat <- fit$sigma_hat
  }
  if (is.null(n)) {
    if (is.null(X)) stop("`n` must be provided when `X` is NULL.")
    n <- nrow(X)
  }

  d       <- nrow(theta_hat)
  emap    <- idx_map(d)
  theta_d <- debias_glasso(sigma_hat, theta_hat)
  var_mat <- estimate_variance_ggm(theta_d)
  p_edges <- compute_ggm_pvalues(theta_d, var_mat, n,
                                 scenario = scenario, emap = emap)

  bh <- bhq_test(p_edges, feature_edges, q)
  c(bh, list(theta_d = theta_d, var_mat = var_mat,
             p_edges = p_edges, emap = emap))
}
