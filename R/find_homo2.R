#' Title
#'
#' @param single.edge
#' @param nodesize
#' @param Emap
#' @param d
#'
#' @return
#' @export
#'
#' @examples
find_homo2 <- function(single.edge, d, Emap, nodesize = 4) {
  ## Find homo
  thetatmp <- matrix(0, d, d)
  homo2 <- 0
  E <- NULL ## the edge and pvalue

  if (length(single.edge) <= 2) {
    return(list(homo.num = homo2, E = E))
  }

  FourCliquesNum <- FiveCliquesNum <- 0
  for (j in 1:nrow(single.edge))
  {
    index <- which(Emap == single.edge[j, 1], arr.ind = T)
    thetatmp[index] <- 1

    node.1st <- index[1, 2]
    node.2nd <- index[1, 1]

    g <- igraph::graph.adjacency(thetatmp, mode = "undirected", diag = FALSE)

    Cliques <- igraph::cliques(g)
    FourCliques <- Cliques[which(sapply(Cliques, length) == 4)]
    FiveCliques <- Cliques[which(sapply(Cliques, length) == 5)]

    Extra4cliques <- length(FourCliques) - FourCliquesNum
    Extra5cliques <- length(FiveCliques) - FiveCliquesNum

    if (Extra4cliques > 0) {
      homo2 <- homo2 + Extra4cliques - Extra5cliques
      E <- rbind(E, c(single.edge[j, ], Extra4cliques - Extra5cliques))
    }

    FourCliquesNum <- length(FourCliques)
    FiveCliquesNum <- length(FiveCliques)
  }
  return(list(homo.num = homo2, E = E))
}
