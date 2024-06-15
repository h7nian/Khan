#' Title
#'
#' @param Node
#' @param Emap
#'
#' @return
#' @export
#'
#' @examples
find_edge <- function(Node, Emap) {
  Edges <- NULL
  for (i in 1:(ncol(Node) - 1)) {
    Edges <- cbind(Edges, Emap[cbind(Node[, i], Node[, i + 1])])
  }
  Edges <- cbind(Edges, Emap[cbind(Node[, ncol(Node)], Node[, 1])])
  Edges <- t(apply(Edges, 1, sort))
  Edges <- unique(Edges)
}
