#' Title
#'
#' @param E.rej
#' @param E.star
#' @param pvals
#' @param J
#' @param qlist
#'
#' @return
#' @export
#'
#' @examples
MyStat <- function(E.rej, E.star, pvals, J, qlist) {
  E.rej <- apply(E.rej, 1, function(e) {
    paste(e, collapse = ",")
  })
  E.star <- apply(E.star, 1, function(e) {
    paste(e, collapse = ",")
  })

  n1 <- length(E.star)
  o.match <- E.rej %in% E.star

  fdps <- powers <- c()
  for (q in qlist) {
    jmax <- max(which(sort(pvals) < (1:length(pvals)) * q / J))
    rej.id <- rank(pvals) < jmax

    n.rej <- sum(rej.id)
    n.rej1 <- sum(rej.id & o.match)

    fdp <- 1 - n.rej1 / n.rej
    power <- n.rej1 / n1

    fdps <- c(fdps, fdp)
    powers <- c(powers, power)
  }
  return(list("fdp" = fdps, "power" = powers))
}
