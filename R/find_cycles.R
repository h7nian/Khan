#' Title
#'
#' @param g
#' @param nodesize
#'
#' @return
#' @export
#'
#' @examples
find_cycles <- function(g, nodesize = 4) {
  Cycles <- NULL
  for (v1 in igraph::V(g)) {
    if (igraph::degree(g, v1, mode = "in") == 0) {
      next
    }
    GoodNeighbors <- igraph::neighbors(g, v1, mode = "out")
    GoodNeighbors <- GoodNeighbors[GoodNeighbors > v1]
    for (v2 in GoodNeighbors) {
      TempCyc <- igraph::all_simple_paths(g, v2, v1, mode = "out", cutoff = nodesize)
      TempCyc <- TempCyc[which(sapply(TempCyc, length) > 2)]
      TempCyc <- TempCyc[sapply(TempCyc, min) == v1]
      Cycles <- c(Cycles, TempCyc)
    }
  }
  Cycles
}
