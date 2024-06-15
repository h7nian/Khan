#' Title
#'
#' @param alphas
#' @param q
#' @param ID
#'
#' @return
#' @export
#'
#' @examples
SumRes <- function(alphas, q, ID) {
  num.rej <- FDP <- Power <- c()
  for (q0 in q) {
    num.rej <- c(num.rej, sum(alphas < q0))
    FDP <- c(FDP, sum(alphas[-ID] < q0) / max(sum(alphas < q0), 1))
    if (length(ID) > 0) {
      Power <- c(Power, sum(alphas[ID] < q0) / length(ID))
    } else {
      Power <- c(Power, 0.99)
    }
  }
  Tab <- cbind(q, num.rej, FDP, Power)
  Tab <- t(Tab)
  res <- c(length(ID), length(ID) / length(alphas), as.vector(Tab))
  names(res) <- c("numH1", "ratioH1", rep(c("q", "numRej", "FDP", "power"), length(q)))
  res
}
