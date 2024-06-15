#' Title
#'
#' @param Sigmahat
#' @param Omegahat
#'
#' @return
#' @export
#'
#' @examples
graph_debias <- function(Sigmahat, Omegahat) {
  dOmega <- matrix(0, nrow = nrow(Sigmahat), ncol = ncol(Sigmahat))
  for (j in 1:nrow(Sigmahat)) {
    for (k in 1:ncol(Sigmahat)) {
      dOmega[j, k] <- Omegahat[j, k] - (t(Omegahat[j, ]) %*% Sigmahat %*% Omegahat[k, ] - Omegahat[j, k]) / t(Omegahat[j, ]) %*% Sigmahat[j, ]
    }
  }
  dOmega
}
