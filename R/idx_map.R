#' Title
#'
#' @param d
#'
#' @return
#' @export
#'
#' @examples
idx_map <- function(d) {
  mat <- matrix(0, d, d)
  mat[upper.tri(mat)] <- 1:(d * (d - 1) / 2)
  mat <- mat + t(mat)
  return(mat)
}
