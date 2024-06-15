#' Title
#'
#' @param data
#' @param Theta.hat
#' @param num.B
#'
#' @return
#' @export
#'
#' @examples
set_bootstrap_normalize_test <- function(data, Theta.hat, num.B) {
  d <- ncol(Theta.hat)
  n <- nrow(data)
  M <- data %*% Theta.hat

  mat.norm <- matrix(stats::rnorm(n * num.B), ncol = n)
  idx <- upper.tri(matrix(0, d, d), diag = FALSE)

  vec <- diag(Theta.hat)
  mat <- matrix(sqrt(vec)) %*% t(sqrt(vec))
  scale_term <- as.vector(mat[idx])

  M.dot <- apply(M, 1, function(xi) {
    (xi %*% t(xi))[idx]
  })
  Theta.hatUp <- Theta.hat[idx]

  bootstrapped.replicas <- mat.norm %*% t(M.dot) - rowSums(mat.norm) %*% t(Theta.hatUp)
  bootstrapped.replicas <- t(bootstrapped.replicas)
  bootstrapped.replicas <- bootstrapped.replicas / scale_term

  bootstrapped.replicas <- t(bootstrapped.replicas) / n
  return(bootstrapped.replicas)
}
