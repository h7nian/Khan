#' Title
#'
#' @param g
#' @param nodesize
#'
#' @return
#' @export
#'
#' @examples
find_trees <- function(g, nodesize = 4) {
  Treesa <- Treesb <- Treesc <- NULL
  for (v1 in igraph::V(g))
  {
    if (igraph::degree(g, v1, mode = "out") == 0) {
      next
    }
    if (igraph::degree(g, v1, mode = "out") >= 4) {
      neigh <- igraph::neighbors(g, v1, mode = "out")
      Tree <- utils::combn(neigh, 4)
      Tree <- split(Tree, col(Tree))
      Tree <- lapply(Tree, function(x) {
        append(v1, x)
      })

      Treesc <- c(Treesc, Tree)
    }
    if (igraph::degree(g, v1, mode = "out") >= 3) {
      neigh <- igraph::neighbors(g, v1, mode = "out")
      goodneigh <- as.numeric(neigh[which(sapply(neigh, FUN = function(x) {
        igraph::degree(g, x, mode = "out")
      }) >= 2)])
      if (length(goodneigh) != 0) {
        for (v in goodneigh)
        {
          v_neibor <- igraph::neighbors(g, v, mode = "out")
          v_neibor <- setdiff(v_neibor, v1)
          Tree <- utils::combn(setdiff(neigh, v), 2)
          Tree <- split(Tree, col(Tree))
          subTree <- v_neibor
          Tree2 <- apply(expand.grid(subTree, Tree), 1, FUN = function(x) {
            c(x$Var1, x$Var2)
          })
          Tree2 <- t(Tree2)
          Tree2 <- Tree2[which(apply(Tree2, 1, FUN = function(x) {
            length(unique(x)) == 3
          })), ]
          if (length(Tree2) == 0) {
            next
          }
          if (!is.matrix(Tree2)) {
            Treesb <- c(Treesb, list(c(v1, v, Tree2[1], Tree2[2], Tree2[3])))
            next
          }
          Treesb <- c(Treesb, lapply(split(Tree2, 1:nrow(Tree2)), FUN = function(x) {
            c(v1, v, x)
          }))
        }
      }
    }
    TempTree <- igraph::all_simple_paths(g, v1, cutoff = nodesize)
    Tree <- TempTree[which(sapply(TempTree, length) == 5)]
    Treesa <- c(Treesa, Tree)
  }
  list(Treesa = Treesa, Treesb = Treesb, Treesc = Treesc)
}
