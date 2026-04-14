# =============================================================================
# ggm.R -- Gaussian graphical model: data generation, estimation, p-values.
#
# Notation follows Liang, Zhang, and Neykov (2024):
#   Theta*   = (Sigma*)^{-1} (true precision),
#   theta_d  = debiased GLasso estimator (equation 2.7),
#   sigma_e  = standard deviation of sqrt(n) * theta_d entries.
# =============================================================================

# ---------------------------------------------------------------------------
# Graph / precision matrix generation
# ---------------------------------------------------------------------------

#' Generate a GGM precision matrix for subgraph feature testing.
#'
#' Builds a block-diagonal precision matrix comprising `m1` triangles,
#' `m2` four-cycles, and `m3` five-cycles, optionally adding extra diagonal
#' edges per block.  Edge weights are sampled uniformly from `weight_range`.
#'
#' @param m1            Number of triangles.
#' @param m2            Number of four-cycles.
#' @param m3            Number of five-cycles.
#' @param weight_range  Numeric length-2 vector: range of off-diagonal entries.
#' @param extra_prob_4  Probability of an extra diagonal edge per four-cycle.
#' @param extra_prob_5  Probability of each extra cross-edge per five-cycle.
#' @param v_diag        Constant added to the diagonal for positive definiteness.
#' @return List with `omega` (precision), `theta` (0/1 adjacency), `d`, and
#'   `feature_info`.
#' @export
generate_ggm_subgraph <- function(m1, m2, m3,
                                  weight_range = c(0.85, 1),
                                  extra_prob_4 = 1,
                                  extra_prob_5 = 1,
                                  v_diag = 0.1) {
  d     <- 3L * m1 + 4L * m2 + 5L * m3
  omega <- matrix(0, nrow = d, ncol = d)
  theta <- matrix(0L, nrow = d, ncol = d)
  offset <- 0L

  add_edge <- function(u, v) {
    w <- stats::runif(1L, weight_range[1L], weight_range[2L])
    omega[u, v] <<- w; omega[v, u] <<- w
    theta[u, v] <<- 1L; theta[v, u] <<- 1L
  }

  for (i in seq_len(m1)) {
    nodes <- offset + 1:3
    for (e in list(c(1, 2), c(2, 3), c(1, 3))) add_edge(nodes[e[1]], nodes[e[2]])
    offset <- offset + 3L
  }

  for (i in seq_len(m2)) {
    nodes <- offset + 1:4
    edges <- list(c(1, 2), c(2, 3), c(3, 4), c(1, 4))
    if (stats::runif(1L) < extra_prob_4) edges <- c(edges, list(c(1, 3)))
    for (e in edges) add_edge(nodes[e[1]], nodes[e[2]])
    offset <- offset + 4L
  }

  for (i in seq_len(m3)) {
    nodes <- offset + 1:5
    edges <- list(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(1, 5))
    for (ee in list(c(1, 3), c(2, 4), c(3, 5))) {
      if (stats::runif(1L) < extra_prob_5) edges <- c(edges, list(ee))
    }
    for (e in edges) add_edge(nodes[e[1]], nodes[e[2]])
    offset <- offset + 5L
  }

  min_eig <- min(eigen(omega, symmetric = TRUE, only.values = TRUE)$values)
  diag(omega) <- diag(omega) + abs(min_eig) + v_diag

  list(omega = omega, theta = theta, d = d,
       feature_info = list(m1 = m1, m2 = m2, m3 = m3))
}

#' Generate a GGM precision matrix for persistent homology simulations.
#'
#' Creates clique-based structures (triangles, four- and five-cliques) with
#' random edge reduction in the five-cliques.
#'
#' @inheritParams generate_ggm_subgraph
#' @param reduce_prob Probability of removing each edge in the five-cliques.
#' @return List with `omega`, `theta`, `d`.
#' @export
generate_ggm_homology <- function(m1, m2, m3,
                                  weight_range = c(0, 10),
                                  reduce_prob = 0.1,
                                  v_diag = 0.25) {
  d     <- 3L * m1 + 4L * m2 + 5L * m3
  omega <- matrix(0, nrow = d, ncol = d)
  theta <- matrix(0L, nrow = d, ncol = d)
  offset <- 0L

  add_edge <- function(u, v) {
    w <- stats::runif(1L, weight_range[1L], weight_range[2L])
    omega[u, v] <<- w; omega[v, u] <<- w
    theta[u, v] <<- 1L; theta[v, u] <<- 1L
  }

  for (i in seq_len(m1)) {
    nodes <- offset + 1:3
    for (a in 1:2) for (b in (a + 1L):3L) add_edge(nodes[a], nodes[b])
    offset <- offset + 3L
  }

  for (i in seq_len(m2)) {
    nodes <- offset + 1:4
    for (a in 1:3) for (b in (a + 1L):4L) add_edge(nodes[a], nodes[b])
    offset <- offset + 4L
  }

  for (i in seq_len(m3)) {
    nodes <- offset + 1:5
    for (a in 1:4) for (b in (a + 1L):5L) {
      if (stats::runif(1L) > reduce_prob) add_edge(nodes[a], nodes[b])
    }
    offset <- offset + 5L
  }

  min_eig <- min(eigen(omega, symmetric = TRUE, only.values = TRUE)$values)
  diag(omega) <- diag(omega) + abs(min_eig) + v_diag

  list(omega = omega, theta = theta, d = d)
}

# ---------------------------------------------------------------------------
# Data sampling
# ---------------------------------------------------------------------------

#' Sample n observations from `N(0, Omega^{-1})`.
#'
#' @param n     Sample size.
#' @param omega Precision matrix (positive definite).
#' @return List with `X` (n-by-d matrix), `sigma` (true covariance),
#'   `sigma_hat` (sample covariance `X'X / n`).
#' @export
generate_ggm_data <- function(n, omega) {
  sigma     <- solve(omega)
  X         <- MASS::mvrnorm(n, mu = rep(0, nrow(omega)), Sigma = sigma)
  sigma_hat <- crossprod(X) / n
  list(X = X, sigma = sigma, sigma_hat = sigma_hat)
}

# ---------------------------------------------------------------------------
# Estimation: GLasso + debiasing (Section 2.3, equations 2.7-2.8)
# ---------------------------------------------------------------------------

#' Estimate the precision matrix via cross-validated graphical lasso.
#'
#' Wraps [CVglasso::CVglasso()] on the sample covariance.
#'
#' @param X   n-by-d data matrix.
#' @param ... Additional arguments forwarded to [CVglasso::CVglasso()].
#' @return List with `theta_hat` and `sigma_hat`.
#' @export
estimate_ggm <- function(X, ...) {
  S   <- crossprod(X) / nrow(X)
  fit <- CVglasso::CVglasso(S = S, path = TRUE, trace = "none", ...)
  list(theta_hat = fit$Omega, sigma_hat = fit$Sigma)
}

#' Debiased GLasso estimator (equation 2.7).
#'
#' Computes
#' \deqn{\hat{\Theta}^d_{uv} = \hat{\Theta}_{uv} -
#'   \hat{\Theta}_{u \cdot}^\top
#'   (\hat{\Sigma} \hat{\Theta}_{\cdot v} - e_v) /
#'   (\hat{\Theta}_{u \cdot}^\top \hat{\Sigma}_{\cdot u}).}
#' Uses matrix products (`theta_hat %*% sigma_hat`, `sigma_hat %*% theta_hat`)
#' so the per-entry work is a dot product rather than the nested scalar loop
#' in the original implementation.
#'
#' @param sigma_hat Estimated (or sample) covariance matrix.
#' @param theta_hat Estimated precision matrix.
#' @return Debiased precision matrix (d-by-d).
#' @export
debias_glasso <- function(sigma_hat, theta_hat) {
  d <- nrow(theta_hat)
  theta_d <- matrix(0, nrow = d, ncol = d)

  ts <- theta_hat %*% sigma_hat
  st <- sigma_hat %*% theta_hat

  for (u in seq_len(d)) {
    denom <- ts[u, u]
    for (v in seq_len(d)) {
      # Numerator: theta_hat[u, ] %*% (sigma_hat %*% theta_hat[, v] - e_v)
      #         = theta_hat[u, ] %*% st[, v] - theta_hat[u, v]
      numer <- sum(theta_hat[u, ] * st[, v]) - theta_hat[u, v]
      theta_d[u, v] <- theta_hat[u, v] - numer / denom
    }
  }
  theta_d
}

#' Variance estimate for the debiased estimator (equation 2.7).
#'
#' `var_mat[u, v] = theta_d[u, u] * theta_d[v, v] + theta_d[u, v]^2`.
#'
#' @param theta_d Debiased precision matrix.
#' @return Symmetric d-by-d variance matrix (diagonal not used).
#' @export
estimate_variance_ggm <- function(theta_d) {
  d <- nrow(theta_d)
  var_mat <- matrix(0, nrow = d, ncol = d)
  for (u in seq_len(d - 1L)) {
    for (v in (u + 1L):d) {
      var_mat[u, v] <- theta_d[u, u] * theta_d[v, v] + theta_d[u, v]^2
      var_mat[v, u] <- var_mat[u, v]
    }
  }
  var_mat
}

# ---------------------------------------------------------------------------
# p-value computation (equation 2.5)
# ---------------------------------------------------------------------------

#' Edge-level p-values for the debiased GGM estimator.
#'
#' Scenario `"a"`: two-sided, `p_e = 2 - 2 Phi(|sqrt(n) theta_d_e / sigma_e|)`.
#' Scenario `"b"`: one-sided, `p_e = 1 - Phi(sqrt(n) theta_d_e / sigma_e)`.
#'
#' @param theta_d  Debiased precision matrix.
#' @param var_mat  Variance matrix from [estimate_variance_ggm()].
#' @param n        Sample size.
#' @param scenario `"a"` or `"b"`.
#' @param emap     Edge-index map from [idx_map()] (required).
#' @return Numeric vector of length `d*(d-1)/2`, indexed by edge index.
#' @export
compute_ggm_pvalues <- function(theta_d, var_mat, n, scenario = "a", emap = NULL) {
  if (is.null(emap)) stop("`emap` must be provided explicitly.")
  d <- nrow(theta_d)
  n_edges  <- d * (d - 1L) / 2L
  p_values <- rep(1.0, n_edges)

  for (u in seq_len(d - 1L)) {
    for (v in (u + 1L):d) {
      idx <- emap[u, v]
      se  <- sqrt(var_mat[u, v])
      if (se < 1e-15) {
        p_values[idx] <- 1.0
        next
      }
      stat <- sqrt(n) * theta_d[u, v] / se
      p_values[idx] <- if (scenario == "a") {
        2 - 2 * stats::pnorm(abs(stat))
      } else {
        1 - stats::pnorm(stat)
      }
    }
  }
  p_values
}

# ---------------------------------------------------------------------------
# Dominating error term xi_i(e) (equation 4.1) -- dependence diagnostics
# ---------------------------------------------------------------------------

#' Compute the i.i.d. summands xi_i(e) for every edge.
#'
#' For the GGM, `xi_i(e = (u, v)) = (X_i^\top Theta*_u)(X_i^\top Theta*_v) -
#' I(u == v)`.  Uses the precomputed product `TX = X %*% theta_true`.
#'
#' @param X          n-by-d data matrix.
#' @param theta_true True precision matrix.
#' @param emap       Edge-index map.
#' @return Matrix of dimension `n` by `d*(d-1)/2`; column `e` holds
#'   `xi_1(e), ..., xi_n(e)`.
#' @export
compute_xi_summands <- function(X, theta_true, emap) {
  n <- nrow(X)
  d <- ncol(X)
  n_edges <- d * (d - 1L) / 2L
  xi <- matrix(0, nrow = n, ncol = n_edges)

  TX <- X %*% theta_true
  for (u in seq_len(d - 1L)) {
    for (v in (u + 1L):d) {
      idx <- emap[u, v]
      xi[, idx] <- TX[, u] * TX[, v]
    }
  }
  xi
}
