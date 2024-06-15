#' Title
#'
#' @param SigEdge
#' @param p
#'
#' @return
#' @export
#'
#' @examples
find_allpaths <- function(SigEdge, p) {
  d <- nrow(SigEdge)
  size <- d / p

  E3 <- E4 <- E5 <- E4.cli <- E5.cli <- NULL
  for (j in 1:p)
  {
    subID <- ((j - 1) * size + 1):(j * size)
    SigEdgeSub <- SigEdge[subID, subID]
    if (sum(SigEdgeSub) < 3) {
      next
    }

    g.hat.sub <- igraph::graph.adjacency(SigEdgeSub, mode = "undirected", diag = FALSE)
    Paths <- find_paths(g.hat.sub, nodesize = 4)
    if (is.null(Paths)) {
      next
    }
    paths.len <- sapply(Paths, length)

    ThreeNode <- do.call(rbind, Paths[paths.len == 3])
    if (length(ThreeNode) > 0) {
      E3 <- rbind(E3, find_edge(ThreeNode))
    }

    FourNode <- do.call(rbind, Paths[paths.len == 4])
    if (length(FourNode) > 0) {
      E4 <- rbind(E4, find_edge(FourNode))
    }

    FiveNode <- do.call(rbind, Paths[paths.len == 5])
    if (length(FiveNode) > 0) {
      E5 <- rbind(E5, find_edge(FiveNode))
    }
  }

  return(list("E3" = E3, "E4" = E4, "E5" = E5))
}
