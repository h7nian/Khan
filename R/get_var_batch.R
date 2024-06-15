#' Title
#'
#' @param data
#' @param Theta.hat
#' @param batch.size
#'
#' @return
#' @export
#'
#' @examples
get_var_batch <- function(data, Theta.hat, batch.size = 32) {
  d <- ncol(Theta.hat)
  n <- nrow(data)

  # calculate number of batches
  num.batches <- ceiling(n / batch.size)

  # initialize variance vector
  Var.e <- rep(0, d)

  # iterate over each batch
  for (i in 1:num.batches)
  {
    print(i)
    # select batch indices
    start.idx <- (i - 1) * batch.size + 1
    end.idx <- min(i * batch.size, n)
    batch <- data[start.idx:end.idx, ]

    # calculate M.dot
    M <- batch %*% Theta.hat
    idx <- upper.tri(matrix(0, d, d), diag = FALSE)
    M.dot <- apply(M, 1, function(xi) {
      (xi %*% t(xi))[idx]
    })

    # calculate scale term
    Theta.hatUp <- Theta.hat[idx]
    vec <- diag(Theta.hat)
    mat <- matrix(sqrt(vec)) %*% t(sqrt(vec))
    scale_term <- as.vector(mat[idx])

    # update variance vector
    Var.e <- Var.e + rowSums((M.dot - Theta.hatUp)^2) / scale_term^2 / n^2

    # remove temporary variables
    rm(batch, M, M.dot, Theta.hatUp, vec, mat, scale_term)
  }
  return(Var.e)
}
