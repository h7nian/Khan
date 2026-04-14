# =============================================================================
# homology.R -- Homology group rank computation.
#
# Implements rank(Z_k(E)) for k = 1 (fast Union-Find) and k = 2 (boundary
# operator via QR on the clique complex).  Used inside dgs() and the generic
# rank_increase() fallback for k >= 2.
# =============================================================================

#' Rank of the k-th cycle group Z_k(E).
#'
#' For `k = 1` this is `|E| - |V_E| + c(G)`, computed in near-linear time by
#' Union-Find.  For `k = 2` the rank of the boundary operator on the
#' 2-skeleton of the clique complex is computed by QR.
#'
#' @param edge_set Integer vector of edge indices.
#' @param d        Number of nodes.
#' @param emap     Edge-index map.
#' @param k        Homology dimension (1 or 2).
#' @param rev_emap Optional reverse edge map from [build_reverse_emap()];
#'   rebuilt if `NULL`.
#' @return Integer rank.
#' @export
compute_homology_rank <- function(edge_set, d, emap, k = 1L, rev_emap = NULL) {
  if (length(edge_set) == 0L) return(0L)

  if (k == 1L) {
    return(.homology_rank_k1(edge_set, d, emap, rev_emap))
  }
  if (k == 2L) {
    return(.homology_rank_k2(edge_set, d, emap, rev_emap))
  }
  stop("Homology dimension k > 2 not implemented.")
}

# Fast H1 rank via Union-Find: beta_1 = |E| - |V_E| + c(G).
.homology_rank_k1 <- function(edge_set, d, emap, rev_emap = NULL) {
  if (is.null(rev_emap)) rev_emap <- build_reverse_emap(emap)

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

  n_cycles <- 0L
  for (idx in edge_set) {
    uv <- rev_emap[[idx]]
    ru <- uf_find(uv[1L])
    rv <- uf_find(uv[2L])
    if (ru == rv) {
      n_cycles <- n_cycles + 1L
    } else {
      if (uf_rank[ru] < uf_rank[rv]) {
        parent[ru] <- rv
      } else if (uf_rank[ru] > uf_rank[rv]) {
        parent[rv] <- ru
      } else {
        parent[rv] <- ru
        uf_rank[ru] <- uf_rank[ru] + 1L
      }
    }
  }
  n_cycles
}

# H2 rank via boundary operator on the clique complex.
.homology_rank_k2 <- function(edge_set, d, emap, rev_emap = NULL) {
  if (is.null(rev_emap)) rev_emap <- build_reverse_emap(emap)

  adj <- matrix(0L, nrow = d, ncol = d)
  for (idx in edge_set) {
    uv <- rev_emap[[idx]]
    adj[uv[1L], uv[2L]] <- 1L
    adj[uv[2L], uv[1L]] <- 1L
  }

  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  triangles <- igraph::cliques(g, min = 3L, max = 3L)
  if (length(triangles) == 0L) return(0L)

  n_tri <- length(triangles)
  n_edg <- length(edge_set)
  boundary <- matrix(0, nrow = n_edg, ncol = n_tri)

  edge_lookup <- integer(max(edge_set))
  edge_lookup[] <- NA_integer_
  for (i in seq_along(edge_set)) edge_lookup[edge_set[i]] <- i

  for (j in seq_len(n_tri)) {
    tri_nodes <- sort(as.integer(triangles[[j]]))
    e1 <- emap[tri_nodes[1L], tri_nodes[2L]]
    e2 <- emap[tri_nodes[1L], tri_nodes[3L]]
    e3 <- emap[tri_nodes[2L], tri_nodes[3L]]
    if (!is.na(edge_lookup[e1])) boundary[edge_lookup[e1], j] <-  1
    if (!is.na(edge_lookup[e2])) boundary[edge_lookup[e2], j] <- -1
    if (!is.na(edge_lookup[e3])) boundary[edge_lookup[e3], j] <-  1
  }

  as.integer(n_tri - qr(boundary)$rank)
}

#' Increase in homology rank from adding a single edge.
#'
#' @param current_edges Integer vector of current edge indices.
#' @param new_edge      Single edge index to add.
#' @param d             Number of nodes.
#' @param emap          Edge-index map.
#' @param k             Homology dimension.
#' @return Non-negative integer.
#' @export
rank_increase <- function(current_edges, new_edge, d, emap, k = 1L) {
  rank_before <- compute_homology_rank(current_edges, d, emap, k)
  rank_after  <- compute_homology_rank(c(current_edges, new_edge), d, emap, k)
  as.integer(rank_after - rank_before)
}
