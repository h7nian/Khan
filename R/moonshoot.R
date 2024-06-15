#' Title
#'
#' @param X
#' @param SigmaHat
#' @param ThetaHat
#' @param q
#' @param numB
#' @param V0
#' @param possibleSet
#' @param prescreen
#'
#' @return
#' @export
#'
#' @examples
moonshoot <- function(X, SigmaHat, ThetaHat, q, numB = 1000, V0, possibleSet, prescreen = T) {
  d <- ncol(SigmaHat)
  dOmega <- graph_debias(SigmaHat, ThetaHat)
  dOmega <- normalize_mat(dOmega, dOmega)
  cat("Generating bootstrap samples with numB =", numB, "...")
  boot_stats <- set_bootstrap_normalize_test(X, ThetaHat, numB)
  cat("...done", "\n")

  dOmega.up <- dOmega[upper.tri(dOmega, diag = F)]
  ncorrect <- 0
  # The p-value of each edge: to do precreening.
  Single_edge.p.value <- sapply(1:length(dOmega.up), function(i) {
    (ncorrect + sum(abs(boot_stats[, i]) >= abs(dOmega.up[i]))) / (ncorrect + numB)
  })

  J <- nrow(V0)
  m <- ncol(V0)
  idxmap_mat <- idx_map(d)
  # Initialization:
  alphas <- matrix(q + 0.01 + numeric(J), nrow = J, ncol = length(possibleSet))

  # Assign the edge id to each vertex set.
  E.all <- NULL
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      ge <- idxmap_mat[cbind(V0[, i], V0[, j])]
      E.all <- cbind(E.all, ge)
    }
  }

  # Pre-screening
  if (prescreen) {
    E1.id <- Single_edge.p.value < q
    SubH1 <- NULL
    for (j in 1:length(possibleSet)) {
      subS <- 1
      for (i in 1:length(possibleSet[[j]])) {
        subS <- subS * E1.id[E.all[, possibleSet[[j]][i]]]
      }
      SubH1 <- cbind(SubH1, subS)
    }

    SubH1.id <- which(rowSums(SubH1) > 0)
    cat("The number of input subgraph is ", J, " and after pre-screen is ", length(SubH1.id), "\n")

    if (length(SubH1.id) == 0) {
      return(list("pval" = alphas, "ID" = SubH1))
    }
    E.all <- matrix(E.all[SubH1.id, ], nrow = length(SubH1.id))
  } else {
    SubH1.id <- 1:nrow(E.all)
  }

  alpha1 <- t(apply(E.all, 1, function(edges) {
    obs_stats <- dOmega.up[edges]
    stats_ord <- order(abs(obs_stats), decreasing = TRUE)

    max1 <- 0
    alpha0 <- rep(q + 0.01, length(edges))
    for (i in length(edges):1) {
      max1 <- pmax(abs(boot_stats[, edges[stats_ord[i]]]), max1)
      a <- (ncorrect + sum(max1 >= abs(obs_stats[stats_ord[i]]))) / (ncorrect + numB)
      alpha0[stats_ord[i]] <- a
    }
    unlist(lapply(possibleSet, function(setj) {
      # max(alpha0[stats_ord[1:max(match(setj,stats_ord))]])
      max(alpha0[setj])
    }))
  }))

  alphas[SubH1.id, ] <- alpha1

  return(list("pval" = alphas, "ID" = SubH1, "alpha1" = alpha1))
}
