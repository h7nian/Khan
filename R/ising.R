# =============================================================================
# ising.R -- Ferromagnetic Ising model: data generation, estimation, p-values.
#
# Implements the Ising model component of Section 2.4 (equation 2.9). The
# distribution is P(x) = (1/Z) exp(sum_{(u,v)} w*_{uv} x_u x_v), w*_{uv} >= 0,
# with edge weight W*_{uv} = E[X_u X_v] - tanh(theta).
# =============================================================================

# ---------------------------------------------------------------------------
# Graph generation for Ising tree simulation
# ---------------------------------------------------------------------------

#' Generate a tree-structured Ising graph.
#'
#' Creates disjoint trees, each with a random number of nodes drawn from
#' `m_range`.  Edge weights are sampled uniformly from `w_range`.  Remaining
#' nodes (if any) are disconnected.
#'
#' @param d               Number of nodes.
#' @param m_range         Integer length-2 vector: range of tree sizes.
#' @param w_range         Numeric length-2 vector: range of edge weights.
#' @param theta_threshold Lower bound theta from equation (2.10); retained for
#'   interface symmetry with the estimation routines.
#' @return List with `omega`, `theta`, `d`, and `tree_info`.
#' @export
generate_ising_trees <- function(d, m_range = c(6L, 10L),
                                 w_range = c(0.9, 1.0),
                                 theta_threshold = 0.45) {
  omega <- matrix(0, nrow = d, ncol = d)
  theta <- matrix(0L, nrow = d, ncol = d)

  offset     <- 0L
  tree_sizes <- integer(0)

  while (offset + m_range[1L] <= d) {
    m <- sample(m_range[1L]:m_range[2L], 1L)
    m <- min(m, d - offset)
    if (m < m_range[1L]) break

    nodes <- offset + seq_len(m)
    tree_edges <- .random_tree_edges(m)

    for (k in seq_len(nrow(tree_edges))) {
      u <- nodes[tree_edges[k, 1L]]
      v <- nodes[tree_edges[k, 2L]]
      w <- stats::runif(1L, w_range[1L], w_range[2L])
      omega[u, v] <- w; omega[v, u] <- w
      theta[u, v] <- 1L; theta[v, u] <- 1L
    }

    tree_sizes <- c(tree_sizes, m)
    offset <- offset + m
  }

  list(omega = omega, theta = theta, d = d,
       tree_info = list(tree_sizes = tree_sizes, n_trees = length(tree_sizes)))
}

# Random tree on m nodes via a random Prufer sequence.
.random_tree_edges <- function(m) {
  if (m <= 1L) return(matrix(integer(0), ncol = 2L))
  if (m == 2L) return(matrix(c(1L, 2L), ncol = 2L))

  prufer <- sample(seq_len(m), m - 2L, replace = TRUE)
  deg    <- rep(1L, m)
  for (p in prufer) deg[p] <- deg[p] + 1L

  edges <- matrix(0L, nrow = m - 1L, ncol = 2L)
  row   <- 1L
  for (p in prufer) {
    for (leaf in seq_len(m)) {
      if (deg[leaf] == 1L) {
        edges[row, ] <- sort(c(leaf, p))
        deg[leaf] <- deg[leaf] - 1L
        deg[p]    <- deg[p] - 1L
        row <- row + 1L
        break
      }
    }
  }
  edges[row, ] <- sort(which(deg == 1L))
  edges
}

# ---------------------------------------------------------------------------
# Data sampling
# ---------------------------------------------------------------------------

#' Sample n i.i.d. observations from a ferromagnetic Ising model.
#'
#' Wraps [IsingSampler::IsingSampler()].  The internal reparameterisation
#' converts between the `{0, 1}` convention used by that package and the
#' `{-1, +1}` convention used in the paper so that zero external field and a
#' pure pairwise coupling `w* x_u x_v` are produced.
#'
#' @param n               Sample size.
#' @param omega           Symmetric non-negative weight matrix.
#' @param theta_threshold Kept for interface symmetry (unused here).
#' @param method          Sampling method forwarded to `IsingSampler`,
#'   e.g. `"MH"` or `"CFTP"`.
#' @return Integer matrix of dimension `n` by `d` with entries in `{-1, +1}`.
#' @export
generate_ising_data <- function(n, omega, theta_threshold = 0.45, method = "MH") {
  graph      <- 4 * omega
  thresholds <- -2 * rowSums(omega)

  X <- IsingSampler::IsingSampler(n, graph = graph, thresholds = thresholds,
                                  method = method)
  if (min(X) >= 0) X <- 2 * X - 1
  X
}

# ---------------------------------------------------------------------------
# Estimation (Section 2.4, equation 2.12)
# ---------------------------------------------------------------------------

#' Estimate Ising edge weights and edge-variance estimates.
#'
#' `w_hat_{uv}   = mean(X[, u] * X[, v]) - tanh(theta_threshold)`.
#' `var_hat_{uv} = 1 - mean(X[, u] * X[, v])^2`.
#'
#' @param X               n-by-d matrix with entries in `{-1, +1}`.
#' @param theta_threshold Value of theta from equation (2.10).
#' @param emap            Edge-index map.
#' @return List with `w_hat` (numeric vector) and `var_hat` (numeric vector),
#'   both of length `d*(d-1)/2`.
#' @export
estimate_ising_weights <- function(X, theta_threshold = 0.45, emap) {
  d       <- ncol(X)
  n_edges <- d * (d - 1L) / 2L
  w_hat   <- numeric(n_edges)
  var_hat <- numeric(n_edges)
  tanh_th <- tanh(theta_threshold)

  for (u in seq_len(d - 1L)) {
    for (v in (u + 1L):d) {
      idx <- emap[u, v]
      m_uv <- mean(X[, u] * X[, v])
      w_hat[idx]   <- m_uv - tanh_th
      var_hat[idx] <- 1 - m_uv^2
    }
  }
  list(w_hat = w_hat, var_hat = var_hat)
}

# ---------------------------------------------------------------------------
# p-value computation (equation 2.5, scenario b)
# ---------------------------------------------------------------------------

#' One-sided p-values for Ising edges.
#'
#' `p_e = 1 - Phi(sqrt(n) * w_hat_e / sigma_e)`.
#'
#' @param w_hat   Edge-weight estimates.
#' @param var_hat Edge-variance estimates.
#' @param n       Sample size.
#' @return Numeric vector of p-values.
#' @export
compute_ising_pvalues <- function(w_hat, var_hat, n) {
  se   <- sqrt(pmax(var_hat, 1e-15))
  stat <- sqrt(n) * w_hat / se
  1 - stats::pnorm(stat)
}

# ---------------------------------------------------------------------------
# Ising correlation verification
# ---------------------------------------------------------------------------

#' Verify the Neykov-Liu characterisation for Ising models.
#'
#' For every node pair, compares the empirical correlation
#' `mean(X_u X_v) - tanh(theta)` to the true weight `omega[u, v]`.
#'
#' @param X               n-by-d Ising samples (`{-1, +1}`).
#' @param omega           True weight matrix.
#' @param theta_threshold Theta parameter.
#' @return Data frame with columns `u`, `v`, `is_edge`, `w_true`,
#'   `corr_empirical`, `w_hat`, `correctly_classified`.
#' @export
verify_ising_correlation <- function(X, omega, theta_threshold = 0.45) {
  d       <- ncol(X)
  tanh_th <- tanh(theta_threshold)

  rows <- list()
  for (u in seq_len(d - 1L)) {
    for (v in (u + 1L):d) {
      m_uv    <- mean(X[, u] * X[, v])
      w_hat   <- m_uv - tanh_th
      w_true  <- omega[u, v]
      is_edge <- w_true > 0
      rows[[length(rows) + 1L]] <- data.frame(
        u = u, v = v,
        is_edge = is_edge,
        w_true = w_true,
        corr_empirical = m_uv,
        w_hat = w_hat,
        correctly_classified = (w_hat > 0) == is_edge,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Dominating error term for Ising (equation 4.2)
# ---------------------------------------------------------------------------

#' Compute xi summands for Ising model.
#'
#' `xi_i(e) = X_{iu} X_{iv} - E[X_u X_v]` with `E[X_u X_v] = tanh(theta) + w*`.
#'
#' @param X               n-by-d Ising data.
#' @param omega           True weight matrix.
#' @param theta_threshold Theta parameter.
#' @param emap            Edge-index map.
#' @return Matrix `n` by `d*(d-1)/2`.
#' @export
compute_xi_summands_ising <- function(X, omega, theta_threshold, emap) {
  n       <- nrow(X)
  d       <- ncol(X)
  n_edges <- d * (d - 1L) / 2L
  xi      <- matrix(0, nrow = n, ncol = n_edges)
  tanh_th <- tanh(theta_threshold)

  for (u in seq_len(d - 1L)) {
    for (v in (u + 1L):d) {
      idx <- emap[u, v]
      e_xu_xv <- tanh_th + omega[u, v]
      xi[, idx] <- X[, u] * X[, v] - e_xu_xv
    }
  }
  xi
}
