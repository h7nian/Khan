# =============================================================================
# edge_index.R -- Edge indexing utilities.
#
# A consistent mapping between node pairs (u, v) with u < v and a linear edge
# index in 1..d*(d-1)/2 for the upper triangle of a d x d matrix.  All other
# code uses these helpers rather than recomputing indices from scratch.
# =============================================================================

#' Build the d-by-d edge-index map.
#'
#' Entry `(u, v)` with `u < v` holds the linear index of edge `{u, v}` in
#' `1..d*(d-1)/2`; the lower triangle is filled symmetrically and the
#' diagonal is zero.
#'
#' @param d Number of nodes (positive integer).
#' @return Integer matrix of dimension `d` by `d`.
#' @export
#' @examples
#' emap <- idx_map(4)
#' emap[1, 2]  # edge index of (1, 2)
idx_map <- function(d) {
  d <- as.integer(d)
  stopifnot(d >= 2L)
  emap <- matrix(0L, nrow = d, ncol = d)
  idx  <- 1L
  for (u in seq_len(d - 1L)) {
    for (v in (u + 1L):d) {
      emap[u, v] <- idx
      emap[v, u] <- idx
      idx <- idx + 1L
    }
  }
  emap
}

#' Convert linear edge indices to node pairs.
#'
#' @param edge_idx Integer vector of edge indices in `1..d*(d-1)/2`.
#' @param d        Number of nodes.
#' @return Integer matrix with columns `u` and `v` (`u < v`).
#' @export
edge_to_nodes <- function(edge_idx, d) {
  d <- as.integer(d)
  n_edges <- d * (d - 1L) / 2L
  stopifnot(all(edge_idx >= 1L & edge_idx <= n_edges))

  u <- integer(length(edge_idx))
  v <- integer(length(edge_idx))
  for (k in seq_along(edge_idx)) {
    idx <- edge_idx[k]
    cum <- 0L
    for (i in seq_len(d - 1L)) {
      next_cum <- cum + (d - i)
      if (idx <= next_cum) {
        u[k] <- i
        v[k] <- i + (idx - cum)
        break
      }
      cum <- next_cum
    }
  }
  cbind(u = u, v = v)
}

#' Look up edge indices for given node pairs.
#'
#' @param u    Integer vector of first endpoints.
#' @param v    Integer vector of second endpoints (same length as `u`).
#' @param emap Edge-index map from [idx_map()].
#' @return Integer vector of edge indices.
#' @export
nodes_to_edge <- function(u, v, emap) {
  d <- nrow(emap)
  stopifnot(all(u >= 1L & u <= d), all(v >= 1L & v <= d), all(u != v))
  emap[cbind(u, v)]
}

#' Total number of undirected edges in a complete graph.
#'
#' @param d Number of nodes.
#' @return Integer.
#' @export
n_edges_total <- function(d) {
  as.integer(d * (d - 1L) / 2L)
}

#' Build a fast reverse lookup from edge index to node pair.
#'
#' Returns a list indexed by edge index whose elements are the two-element
#' integer vector `c(u, v)` with `u < v`.  Much faster than repeated
#' `which(emap == idx, arr.ind = TRUE)` calls.
#'
#' @param emap Edge-index map from [idx_map()].
#' @return List of length `n_edges_total(d)`; each element is `integer(2)`.
#' @export
build_reverse_emap <- function(emap) {
  d <- nrow(emap)
  n_edges <- d * (d - 1L) / 2L
  rev_map <- vector("list", n_edges)
  for (u in seq_len(d - 1L)) {
    for (v in (u + 1L):d) {
      rev_map[[emap[u, v]]] <- c(u, v)
    }
  }
  rev_map
}
