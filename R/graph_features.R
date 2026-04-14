# =============================================================================
# graph_features.R -- Subgraph enumeration (cycles, cliques, trees, paths).
#
# Every function takes an explicit `emap` (from idx_map()) argument instead of
# reading it from the enclosing environment.  Feature edges are always
# returned as integer vectors of edge indices so downstream code (bhq_test,
# dgs, ...) can index edge-level vectors directly.
# =============================================================================

# ---------------------------------------------------------------------------
# Cycle finding
# ---------------------------------------------------------------------------

#' Find all simple cycles up to a given size.
#'
#' @param g        An undirected `igraph` object.
#' @param max_size Maximum number of nodes in a cycle.
#' @return A named list of integer matrices (one entry per cycle length in
#'   `3..max_size`); each row is the ordered node IDs of one cycle.
#' @export
find_cycles <- function(g, max_size = 5L) {
  adj <- as.matrix(igraph::as_adjacency_matrix(g))
  d   <- nrow(adj)
  result <- list()

  for (len in 3:max_size) {
    cycles <- .enumerate_cycles(adj, d, len)
    if (length(cycles) > 0L) {
      result[[as.character(len)]] <- do.call(rbind, cycles)
    }
  }
  result
}

.enumerate_cycles <- function(adj, d, len) {
  cycles <- list()
  .dfs_cycle <- function(path, start) {
    cur <- path[length(path)]
    if (length(path) == len) {
      if (adj[cur, start] == 1L) {
        canon <- .canonicalize_cycle(path)
        key   <- paste(canon, collapse = ",")
        if (is.null(cycles[[key]])) {
          cycles[[key]] <<- canon
        }
      }
      return()
    }
    for (v in seq_len(d)) {
      if (v > start && !(v %in% path) && adj[cur, v] == 1L) {
        .dfs_cycle(c(path, v), start)
      }
    }
  }
  for (s in seq_len(d)) {
    .dfs_cycle(s, s)
  }
  unname(cycles)
}

.canonicalize_cycle <- function(cyc) {
  n   <- length(cyc)
  idx <- which.min(cyc)
  cyc <- cyc[c(idx:n, if (idx > 1L) 1:(idx - 1L))]
  if (n >= 3L && cyc[2L] > cyc[n]) {
    cyc <- c(cyc[1L], rev(cyc[2L:n]))
  }
  cyc
}

# ---------------------------------------------------------------------------
# Clique finding
# ---------------------------------------------------------------------------

#' Find all cliques of specified sizes.
#'
#' @param g        An `igraph` object.
#' @param min_size Minimum clique size.
#' @param max_size Maximum clique size.
#' @return Named list keyed by clique size; each element is a matrix whose
#'   rows are the node IDs of one clique.
#' @export
find_cliques_by_size <- function(g, min_size = 3L, max_size = 5L) {
  all_cliques <- igraph::cliques(g, min = min_size, max = max_size)
  result <- list()
  for (cl in all_cliques) {
    nodes <- sort(as.integer(cl))
    key   <- as.character(length(nodes))
    if (is.null(result[[key]])) result[[key]] <- list()
    result[[key]] <- c(result[[key]], list(nodes))
  }
  for (key in names(result)) {
    result[[key]] <- do.call(rbind, result[[key]])
  }
  result
}

# ---------------------------------------------------------------------------
# Five-node tree shapes (Ising model, Figure 7b)
# ---------------------------------------------------------------------------

#' Find the three types of 5-node trees in a graph.
#'
#' Type a: path `1-2-3-4-5` (max degree 2).
#' Type b: one node of degree 3 plus leaves.
#' Type c: star with center of degree 4.
#'
#' @param g An `igraph` object.
#' @return Named list with components `type_a`, `type_b`, `type_c`; each is
#'   an integer matrix with five columns whose rows are sorted node IDs.
#' @export
find_trees_5node <- function(g) {
  adj <- as.matrix(igraph::as_adjacency_matrix(g))
  d   <- nrow(adj)

  type_a <- list()
  type_b <- list()
  type_c <- list()

  combs <- utils::combn(d, 5L)
  for (i in seq_len(ncol(combs))) {
    nodes  <- combs[, i]
    sub_adj <- adj[nodes, nodes]
    n_edges <- sum(sub_adj) / 2L
    if (n_edges != 4L) next

    sub_g <- igraph::graph_from_adjacency_matrix(sub_adj, mode = "undirected")
    if (!igraph::is_connected(sub_g)) next

    degs    <- igraph::degree(sub_g)
    max_deg <- max(degs)

    if (max_deg == 2L) {
      type_a <- c(type_a, list(sort(nodes)))
    } else if (max_deg == 3L) {
      type_b <- c(type_b, list(sort(nodes)))
    } else if (max_deg == 4L) {
      type_c <- c(type_c, list(sort(nodes)))
    }
  }

  list(
    type_a = if (length(type_a) > 0L) do.call(rbind, type_a) else matrix(integer(0), ncol = 5L),
    type_b = if (length(type_b) > 0L) do.call(rbind, type_b) else matrix(integer(0), ncol = 5L),
    type_c = if (length(type_c) > 0L) do.call(rbind, type_c) else matrix(integer(0), ncol = 5L)
  )
}

# ---------------------------------------------------------------------------
# Edge-set extraction from node-level feature descriptions
# ---------------------------------------------------------------------------

#' Edge indices of a cycle.
#'
#' @param nodes Integer vector of node IDs, in cyclic order.
#' @param emap  Edge-index map.
#' @return Integer vector of edge indices.
#' @export
edges_from_cycle <- function(nodes, emap) {
  n <- length(nodes)
  edges <- integer(n)
  for (k in seq_len(n)) {
    u <- nodes[k]
    v <- nodes[if (k < n) k + 1L else 1L]
    edges[k] <- emap[u, v]
  }
  sort(unique(edges))
}

#' Edge indices of an open path.
#'
#' @param nodes Integer vector of node IDs in path order.
#' @param emap  Edge-index map.
#' @return Integer vector of edge indices.
#' @export
edges_from_path <- function(nodes, emap) {
  n <- length(nodes)
  edges <- integer(n - 1L)
  for (k in seq_len(n - 1L)) {
    edges[k] <- emap[nodes[k], nodes[k + 1L]]
  }
  sort(unique(edges))
}

#' Edge indices of a clique.
#'
#' @param nodes Integer vector of node IDs.
#' @param emap  Edge-index map.
#' @return Integer vector of edge indices.
#' @export
edges_from_clique <- function(nodes, emap) {
  n <- length(nodes)
  edges <- integer(0)
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      edges <- c(edges, emap[nodes[i], nodes[j]])
    }
  }
  sort(unique(edges))
}

#' Edge indices of a 5-node tree.
#'
#' @param nodes Integer vector of 5 node IDs.
#' @param adj   Full adjacency matrix.
#' @param emap  Edge-index map.
#' @return Integer vector of edge indices.
#' @export
edges_from_tree <- function(nodes, adj, emap) {
  edges <- integer(0)
  for (i in seq_len(4L)) {
    for (j in (i + 1L):5L) {
      if (adj[nodes[i], nodes[j]] == 1L) {
        edges <- c(edges, emap[nodes[i], nodes[j]])
      }
    }
  }
  sort(unique(edges))
}

# ---------------------------------------------------------------------------
# Master shape finder
# ---------------------------------------------------------------------------

#' Enumerate all subgraph features of the requested types.
#'
#' @param adj            Adjacency matrix (0/1, symmetric, zero diagonal).
#' @param emap           Edge-index map from [idx_map()].
#' @param feature_types  Character vector, a subset of `c("triangle",
#'   "four_cycle", "five_cycle", "four_clique", "five_clique", "tree_a",
#'   "tree_b", "tree_c")`.
#' @return Named list; each component is a list of integer vectors of edge
#'   indices, one per feature instance.
#' @export
find_all_features <- function(adj, emap,
                              feature_types = c("triangle", "four_cycle", "five_cycle")) {
  g   <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  out <- list()

  cycle_types  <- intersect(feature_types, c("triangle", "four_cycle", "five_cycle"))
  clique_types <- intersect(feature_types, c("four_clique", "five_clique"))
  tree_types   <- intersect(feature_types, c("tree_a", "tree_b", "tree_c"))

  # Cycles
  if (length(cycle_types) > 0L) {
    max_cyc <- 3L
    if ("four_cycle" %in% cycle_types) max_cyc <- max(max_cyc, 4L)
    if ("five_cycle" %in% cycle_types) max_cyc <- max(max_cyc, 5L)
    cycs <- find_cycles(g, max_size = max_cyc)
    type_map <- c(triangle = "3", four_cycle = "4", five_cycle = "5")
    for (ft in cycle_types) {
      key <- type_map[ft]
      if (!is.null(cycs[[key]])) {
        mat <- cycs[[key]]
        out[[ft]] <- lapply(seq_len(nrow(mat)), function(r) edges_from_cycle(mat[r, ], emap))
      } else {
        out[[ft]] <- list()
      }
    }
  }

  # Cliques
  if (length(clique_types) > 0L) {
    min_cl <- if ("four_clique" %in% clique_types) 4L else 5L
    max_cl <- if ("five_clique" %in% clique_types) 5L else 4L
    cls <- find_cliques_by_size(g, min_size = min_cl, max_size = max_cl)
    type_map <- c(four_clique = "4", five_clique = "5")
    for (ft in clique_types) {
      key <- type_map[ft]
      if (!is.null(cls[[key]])) {
        mat <- cls[[key]]
        out[[ft]] <- lapply(seq_len(nrow(mat)), function(r) edges_from_clique(mat[r, ], emap))
      } else {
        out[[ft]] <- list()
      }
    }
  }

  # Trees
  if (length(tree_types) > 0L) {
    trees <- find_trees_5node(g)
    type_map <- c(tree_a = "type_a", tree_b = "type_b", tree_c = "type_c")
    for (ft in tree_types) {
      key <- type_map[ft]
      mat <- trees[[key]]
      if (nrow(mat) > 0L) {
        out[[ft]] <- lapply(seq_len(nrow(mat)), function(r) edges_from_tree(mat[r, ], adj, emap))
      } else {
        out[[ft]] <- list()
      }
    }
  }

  out
}
