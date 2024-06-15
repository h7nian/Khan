#' Title
#'
#' @param Node
#' @param Emap
#'
#' @return
#' @export
#'
#' @examples
find_edge_cli <- function(Node, Emap) {
  Edges <- t(apply(Node, 1, function(node) {
    tmp <- Emap[node, node]
    sort(tmp[upper.tri(tmp)])
  }))
  Edges <- unique(Edges)
  Edges
}
