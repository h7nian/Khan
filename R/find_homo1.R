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
find_homo1 <- function(single.edge, d, Emap, nodesize = 4) {
  ## Find homo
  thetatmp <- matrix(0, d, d)
  homo1 <- 0
  E <- NULL ## the edge and pvalue

  if (length(single.edge) <= 2) {
    return(list(homo.num = homo1, E = E))
  }

  for (j in 1:nrow(single.edge))
  {
    print(j)
    index <- which(Emap == single.edge[j, 1], arr.ind = T)
    thetatmp[index] <- 1

    node.1st <- index[1, 2]
    node.2nd <- index[1, 1]

    g <- igraph::graph.adjacency(thetatmp, mode = "undirected", diag = FALSE)

    TempCyc <- NULL
    TempCyc <- igraph::all_simple_paths(g, node.1st, node.2nd, mode = "out", cutoff = nodesize)
    TempCyc <- TempCyc[which(sapply(TempCyc, length) > 2)]

    if (length(TempCyc) > 0) {
      homo1 <- homo1 + 1
      E <- rbind(E, c(single.edge[j, ], 1))
    }
  }
  return(list(homo.num = homo1, E = E))
}
