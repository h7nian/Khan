#' Title
#'
#' @param single.edge
#' @param dOmega
#' @param Emap
#' @param nodesize
#' @param d
#'
#' @return
#' @export
#'
#' @examples
find_homo1_realdata <- function(single.edge, dOmega, d, Emap, nodesize = 4) {
  thetatmp <- matrix(0, d, d)
  homo1 <- 0
  E <- NULL ## the edge and pvalue
  loop <- NULL

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
      distancenode.num <- sapply(TempCyc, function(l) sum(l %in% 1:173))
      most.distancenode <- TempCyc[distancenode.num == max(distancenode.num)]

      if (length(most.distancenode) > 1) {
        ## Calculate the mean weight
        weights <- sapply(most.distancenode, function(l) {
          temp <- dOmega[utils::head(l, -1), utils::tail(l, -1)]
          sum(c(dOmega[node.1st, node.2nd], diag(temp))) / length(l)
        })

        filtered_loop <- most.distancenode[which.max(weights)]
      } else {
        filtered_loop <- most.distancenode
      }

      homo1 <- homo1 + 1
      E <- rbind(E, c(single.edge[j, ], 1))
      loop <- c(loop, filtered_loop)
    }
  }
  return(list(homo.num = homo1, E = E, loop = loop))
}
