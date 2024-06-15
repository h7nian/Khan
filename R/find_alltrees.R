#' Title
#'
#' @param SigEdge
#' @param p
#'
#' @return
#' @export
#'
#' @examples
find_alltrees <- function(SigEdge, p) {
  d <- nrow(SigEdge)
  size <- d / p

  E5a <- E5b <- E5c <- NULL
  for (j in 1:p)
  {
    subID <- ((j - 1) * size + 1):(j * size)
    SigEdgeSub <- SigEdge[subID, subID]
    if (sum(SigEdgeSub) < 3) {
      next
    }

    g.hat.sub <- igraph::graph.adjacency(SigEdge, mode = "undirected", diag = FALSE)
    Trees <- find_trees(g.hat.sub, nodesize = 4)
    if (is.null(Trees)) {
      next
    }

    Treesa <- t(sapply(Trees$Treesa, matrix))
    Treesb <- t(sapply(Trees$Treesb, matrix))
    Treesc <- t(sapply(Trees$Treesc, matrix))

    if (length(Treesa) > 0) {
      Treesa <- Treesa + (j - 1) * size
      E5a <- rbind(E5a, find_edge_tree(Treesa, type = "a"))
    }
    if (length(Treesb) > 0) {
      Treesb <- Treesb + (j - 1) * size
      E5b <- rbind(E5b, find_edge_tree(Treesb, type = "b"))
    }

    if (length(Treesc) > 0) {
      Treesc <- Treesc + (j - 1) * size
      E5c <- rbind(E5c, find_edge_tree(Treesc, type = "c"))
    }
  }
  list("E5a" = E5a, "E5b" = E5b, "E5c" = E5c)
}
