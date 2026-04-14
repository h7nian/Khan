# =============================================================================
# dgs.R -- Algorithm 2: Discrete Gram-Schmidt for homology generators.
#
# Selects linearly independent homological generators at a fixed filtration
# level mu, combining p-value screening with incremental rank computation.
# k = 1 uses a Union-Find fast path; k >= 2 falls back to rank_increase().
# =============================================================================

#' Discrete Gram-Schmidt algorithm for homology selection.
#'
#' At a fixed filtration level, edges are processed in increasing p-value
#' order; each edge that increases the cycle-group rank contributes `ell`
#' copies of its p-value to the DGS list, to which BHq is then applied.
#'
#' @param p_edges  Numeric vector of edge p-values.
#' @param emap     Edge-index map.
#' @param d        Number of nodes.
#' @param q        FDR level.
#' @param k        Homology dimension (1 or 2).
#' @param rev_emap Optional reverse edge map (built internally if `NULL`).
#' @return List with `selected_edges`, `homology_rank`, `alpha_hat`,
#'   `pvalue_list`.
#' @export
dgs <- function(p_edges, emap, d, q, k = 1L, rev_emap = NULL) {
  if (k == 1L) return(.dgs_k1_fast(p_edges, emap, d, q, rev_emap))
  .dgs_generic(p_edges, emap, d, q, k)
}

# Fast path for H1 using Union-Find.
.dgs_k1_fast <- function(p_edges, emap, d, q, rev_emap = NULL) {
  e_candidates <- which(p_edges < q)
  if (length(e_candidates) == 0L) {
    return(list(selected_edges = integer(0), homology_rank = 0L,
                alpha_hat = 0, pvalue_list = numeric(0)))
  }

  if (is.null(rev_emap)) rev_emap <- build_reverse_emap(emap)

  cand_order <- e_candidates[order(p_edges[e_candidates])]

  parent  <- seq_len(d)
  uf_rank <- rep(0L, d)

  uf_find <- function(x) {
    root <- x
    while (parent[root] != root) root <- parent[root]
    while (parent[x] != root) {
      next_x <- parent[x]
      parent[x] <<- root
      x <- next_x
    }
    root
  }

  uf_union <- function(x, y) {
    rx <- uf_find(x)
    ry <- uf_find(y)
    if (rx == ry) return(FALSE)
    if (uf_rank[rx] < uf_rank[ry]) {
      parent[rx] <<- ry
    } else if (uf_rank[rx] > uf_rank[ry]) {
      parent[ry] <<- rx
    } else {
      parent[ry] <<- rx
      uf_rank[rx] <<- uf_rank[rx] + 1L
    }
    TRUE
  }

  pvalue_list <- numeric(0)
  for (e_star in cand_order) {
    uv <- rev_emap[[e_star]]
    if (!uf_union(uv[1L], uv[2L])) {
      pvalue_list <- c(pvalue_list, p_edges[e_star])
    }
  }

  bar_j <- length(pvalue_list)
  if (bar_j == 0L) {
    return(list(selected_edges = integer(0), homology_rank = 0L,
                alpha_hat = 0, pvalue_list = numeric(0)))
  }

  bh <- .dgs_bhq(pvalue_list, q, bar_j)
  selected <- if (bh$alpha_hat > 0) {
    e_candidates[p_edges[e_candidates] < bh$alpha_hat]
  } else {
    integer(0)
  }

  list(
    selected_edges = selected,
    homology_rank  = bh$j_max,
    alpha_hat      = bh$alpha_hat,
    pvalue_list    = pvalue_list
  )
}

# Generic DGS for k >= 2, using rank_increase().
.dgs_generic <- function(p_edges, emap, d, q, k) {
  e_candidates <- which(p_edges < q)
  if (length(e_candidates) == 0L) {
    return(list(selected_edges = integer(0), homology_rank = 0L,
                alpha_hat = 0, pvalue_list = numeric(0)))
  }

  cand_order  <- e_candidates[order(p_edges[e_candidates])]
  tilde_e     <- integer(0)
  pvalue_list <- numeric(0)

  for (e_star in cand_order) {
    ell <- rank_increase(tilde_e, e_star, d, emap, k)
    if (ell > 0L) {
      pvalue_list <- c(pvalue_list, rep(p_edges[e_star], ell))
    }
    tilde_e <- c(tilde_e, e_star)
  }

  bar_j <- length(pvalue_list)
  if (bar_j == 0L) {
    return(list(selected_edges = integer(0), homology_rank = 0L,
                alpha_hat = 0, pvalue_list = numeric(0)))
  }

  bh <- .dgs_bhq(pvalue_list, q, bar_j)
  selected <- if (bh$alpha_hat > 0) {
    e_candidates[p_edges[e_candidates] < bh$alpha_hat]
  } else {
    integer(0)
  }

  list(
    selected_edges = selected,
    homology_rank  = bh$j_max,
    alpha_hat      = bh$alpha_hat,
    pvalue_list    = pvalue_list
  )
}

# BHq on a DGS p-value list: uses the non-strict inequality required by DGS.
.dgs_bhq <- function(alphas, q, bar_j) {
  sorted <- sort(alphas)
  j_max  <- 0L
  for (i in seq_along(sorted)) {
    if (sorted[i] <= q * i / bar_j) j_max <- i
  }
  alpha_hat <- if (j_max > 0L) q * j_max / bar_j else 0
  list(j_max = as.integer(j_max), alpha_hat = alpha_hat)
}

#' Edge p-values at a given filtration level mu.
#'
#' Scenario `"a"`: `stat = sqrt(n) * (|w_hat_e| - mu) / sigma_e`,
#' `p_e = 2 - 2 Phi(max(stat, 0))`.
#' Scenario `"b"`: `stat = sqrt(n) * (w_hat_e - mu) / sigma_e`,
#' `p_e = 1 - Phi(stat)`.
#'
#' @param w_hat    Edge-weight estimates.
#' @param sigma_e  Edge standard deviations.
#' @param n        Sample size.
#' @param mu       Filtration level.
#' @param scenario `"a"` (two-sided) or `"b"` (one-sided).
#' @return Numeric vector of p-values, same length as `w_hat`.
#' @export
pvalues_at_mu <- function(w_hat, sigma_e, n, mu, scenario = "a") {
  se <- pmax(sigma_e, 1e-15)
  if (scenario == "a") {
    stat <- sqrt(n) * (abs(w_hat) - mu) / se
    2 - 2 * stats::pnorm(pmax(stat, 0))
  } else {
    stat <- sqrt(n) * (w_hat - mu) / se
    1 - stats::pnorm(stat)
  }
}
