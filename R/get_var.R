#' Title
#'
#' @param data
#' @param Theta.hat
#'
#' @return
#' @export
#'
#' @examples
get_var <- function(data, Theta.hat) {
  d <- ncol(Theta.hat)
  n <- nrow(data)
  M <- data %*% Theta.hat

  idx <- upper.tri(matrix(0, d, d), diag = FALSE)
  M.dot <- apply(M, 1, function(xi) {
    (xi %*% t(xi))[idx]
  })
  Theta.hatUp <- Theta.hat[idx]
  vec <- diag(Theta.hat)
  mat <- matrix(sqrt(vec)) %*% t(sqrt(vec))
  scale_term <- as.vector(mat[idx])

  Var.e <- rowSums((M.dot - Theta.hatUp)^2) / scale_term^2 / n^2
  return(Var.e)
}
