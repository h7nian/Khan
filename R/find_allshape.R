# Find triangle, 4-cycle, 5-cycle, 4-clique, 5-clique
#' Title
#'
#' @param SigEdge
#' @param p
#'
#' @return
#' @export
#'
#' @examples
find_allshape <- function(SigEdge, p) {
  d <- nrow(SigEdge)
  size <- d / p

  E3 <- E4 <- E5 <- E4.cli <- E5.cli <- NULL
  for (j in 1:p) {
    subID <- ((j - 1) * size + 1):(j * size)
    SigEdgeSub <- SigEdge[subID, subID]
    if (sum(SigEdgeSub) < 3) {
      next
    }

    g.hat.sub <- igraph::graph.adjacency(SigEdgeSub, mode = "undirected", diag = FALSE)
    Cycles <- find_cycles(g.hat.sub, nodesize = 4)
    if (is.null(Cycles)) {
      next
    }

    cyc.len <- sapply(Cycles, length)

    Triangle <- do.call(rbind, Cycles[cyc.len == 3])
    if (length(Triangle) > 0) {
      Triangle <- t(apply(Triangle, 1, sort))
      Triangle <- unique(Triangle) + (j - 1) * size
      E3 <- rbind(E3, find_edge(Triangle))
    }

    Fourcyc <- do.call(rbind, Cycles[cyc.len == 4])
    if (length(Fourcyc) > 0) {
      Fourcyc <- Fourcyc + (j - 1) * size
      E4 <- rbind(E4, find_edge(Fourcyc))
    }

    Fivecyc <- do.call(rbind, Cycles[cyc.len == 5])
    if (length(Fivecyc) > 0) {
      Fivecyc <- Fivecyc + (j - 1) * size
      E5 <- rbind(E5, find_edge(Fivecyc))
    }


    Cliques <- igraph::cliques(g.hat.sub, min = 4, max = 5)
    if (length(Cliques) == 0) {
      next
    }

    cli.len <- sapply(Cliques, length)


    Fourcli <- do.call(rbind, Cliques[cli.len == 4])
    if (length(Fourcli) > 0) {
      Fourcli <- Fourcli + (j - 1) * size
      E4.cli <- rbind(E4.cli, find_edge_cli(Fourcli))
    }

    Fivecli <- do.call(rbind, Cliques[cli.len == 5])
    if (length(Fivecli) > 0) {
      Fivecli <- Fivecli + (j - 1) * size
      E5.cli <- rbind(E5.cli, find_edge_cli(Fivecli))
    }
  }

  return(list("E3" = E3, "E4" = E4, "E5" = E5, "E4.cli" = E4.cli, "E5.cli" = E5.cli))
}
