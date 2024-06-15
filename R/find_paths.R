#' Title
#'
#' @param g
#' @param nodesize
#'
#' @return
#' @export
#'
#' @examples
find_paths <- function(g, nodesize = 4) {
  Paths <- NULL
  for (v1 in igraph::V(g))
  {
    if (igraph::degree(g, v1, mode = "out") == 0) {
      next
    }
    TempPath <- igraph::all_simple_paths(g, v1, cutoff = nodesize)
    TempPath <- TempPath[which(sapply(TempPath, length) > 2)]
    Paths <- c(Paths, TempPath)
  }
  return(Paths)
}
