# =============================================================================
# bhq.R -- Algorithm 1: General graph feature selection with FDR control.
#
# Implements the BHq procedure from Section 2.2.  Given edge p-values and a
# candidate feature list, feature-level p-values are aggregated by
# `alpha(F_j) = max_{e in E(F_j)} p_e` and BH thresholding is applied.
# =============================================================================

#' BHq test for graph feature selection (Algorithm 1).
#'
#' @param p_edges       Numeric vector of edge p-values (length `n_edges`).
#' @param feature_edges List of integer vectors; `feature_edges[[j]]` contains
#'   the edge indices of feature `F_j`.
#' @param q             FDR target level in `(0, 1)`.
#' @return List with:
#'   \describe{
#'     \item{`rejected`}{Integer vector of rejected feature indices.}
#'     \item{`alpha_hat`}{Data-driven threshold `hat{alpha}`.}
#'     \item{`psi`}{Binary decision vector of length `J`.}
#'     \item{`alpha_values`}{Feature-level p-values `alpha(F_j)`.}
#'   }
#' @export
bhq_test <- function(p_edges, feature_edges, q) {
  J <- length(feature_edges)
  if (J == 0L) {
    return(list(rejected = integer(0), alpha_hat = 0,
                psi = integer(0), alpha_values = numeric(0)))
  }

  # Line 3: edges with individual p-value below q.
  sig_edges <- which(p_edges < q)

  # Line 4: feature-level p-values, only for features whose edges all pass.
  alpha_values <- rep(1.0, J)
  for (j in seq_len(J)) {
    edges_j <- feature_edges[[j]]
    if (all(edges_j %in% sig_edges)) {
      alpha_values[j] <- max(p_edges[edges_j])
    }
  }

  bh <- .benjamini_hochberg(alpha_values, q, J)

  list(
    rejected     = bh$rejected,
    alpha_hat    = bh$alpha_hat,
    psi          = bh$psi,
    alpha_values = alpha_values
  )
}

# Benjamini-Hochberg thresholding on a vector of p-values.
# Returns j_max = max{0 <= j <= J : alpha_(j) < q * j / J} and the
# corresponding reject set.  The initial j_max = 0 correctly handles the
# case where no index satisfies the strict inequality.
.benjamini_hochberg <- function(pvals, q, J) {
  ord    <- order(pvals)
  sorted <- pvals[ord]

  j_max <- 0L
  for (j in seq_len(J)) {
    if (sorted[j] < q * j / J) j_max <- j
  }

  if (j_max == 0L) {
    return(list(rejected = integer(0), alpha_hat = 0, psi = rep(0L, J)))
  }

  alpha_hat <- q * j_max / J
  psi       <- as.integer(pvals < alpha_hat)
  list(rejected = which(psi == 1L), alpha_hat = alpha_hat, psi = psi)
}

#' Bonferroni correction for graph feature selection.
#'
#' Rejects feature `j` whenever `alpha(F_j) < q / J`.
#'
#' @inheritParams bhq_test
#' @return Same structure as [bhq_test()].
#' @export
bonferroni_test <- function(p_edges, feature_edges, q) {
  J <- length(feature_edges)
  if (J == 0L) {
    return(list(rejected = integer(0), alpha_hat = 0,
                psi = integer(0), alpha_values = numeric(0)))
  }

  alpha_values <- vapply(feature_edges,
                         function(edges_j) max(p_edges[edges_j]),
                         numeric(1L))
  threshold <- q / J
  psi       <- as.integer(alpha_values < threshold)

  list(
    rejected     = which(psi == 1L),
    alpha_hat    = threshold,
    psi          = psi,
    alpha_values = alpha_values
  )
}
