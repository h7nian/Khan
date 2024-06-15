#' Title
#'
#' @param Node
#' @param type
#' @param Emap
#'
#' @return
#' @export
#'
#' @examples
find_edge_tree <- function(Node, type, Emap) {
  Edges <- NULL
  if (type == "a") {
    for (i in 1:(ncol(Node) - 1)) {
      Edges <- cbind(Edges, Emap[cbind(Node[, i], Node[, i + 1])])
    }
  } else if (type == "b") {
    Edges <- cbind(Edges, Emap[cbind(Node[, 1], Node[, 2])])
    Edges <- cbind(Edges, Emap[cbind(Node[, 1], Node[, 4])])
    Edges <- cbind(Edges, Emap[cbind(Node[, 1], Node[, 5])])
    Edges <- cbind(Edges, Emap[cbind(Node[, 2], Node[, 3])])
  } else {
    Edges <- cbind(Edges, Emap[cbind(Node[, 1], Node[, 2])])
    Edges <- cbind(Edges, Emap[cbind(Node[, 1], Node[, 3])])
    Edges <- cbind(Edges, Emap[cbind(Node[, 1], Node[, 4])])
    Edges <- cbind(Edges, Emap[cbind(Node[, 1], Node[, 5])])
  }
  Edges <- t(apply(Edges, 1, sort))
  Edges <- unique(Edges)
  return(Edges)
}
