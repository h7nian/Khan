#' Title
#'
#' @param mat1
#' @param mat2
#'
#' @return
#' @export
#'
#' @examples
normalize_mat <- function(mat1, mat2) {
  vec <- diag(mat2)
  mat <- mat1 / sqrt(vec)
  mat <- t(mat) / sqrt(vec)
  return(t(mat))
}
