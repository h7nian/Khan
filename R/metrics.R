# =============================================================================
# metrics.R -- FDP, Power, uFDP evaluation metrics.
#
# Operates on the decision vectors (`psi`) and/or DGS results returned by
# bhq_test(), dgs(), khan(); no raw p-values are (re-)thresholded here.
# =============================================================================

#' FDP and power for subgraph feature selection.
#'
#' A feature `F_j` is declared true (`H_1`) iff every edge in `E(F_j)`
#' belongs to the true adjacency `true_adj`.
#'
#' @param psi                 Binary decision vector.
#' @param true_features       Unused (kept for backwards compatibility); the
#'   ground-truth is computed from `true_adj` and `candidate_features`.
#' @param candidate_features  List of edge-index vectors (one per candidate).
#' @param true_adj            True adjacency matrix (0/1).
#' @param emap                Edge-index map.
#' @return List with `fdp`, `power`, `n_rejected`, `n_true_positive`,
#'   `n_false_positive`, `n_true_features`, `n_null_features`.
#' @export
compute_fdp_power <- function(psi, true_features, candidate_features,
                              true_adj = NULL, emap = NULL) {
  J <- length(psi)
  rejected   <- which(psi == 1L)
  n_rejected <- length(rejected)

  if (n_rejected == 0L) {
    return(list(fdp = 0, power = 0, n_rejected = 0L,
                n_true_positive = 0L, n_false_positive = 0L,
                n_true_features = 0L, n_null_features = 0L))
  }
  if (is.null(true_adj) || is.null(emap)) {
    stop("`true_adj` and `emap` are required to classify features.")
  }

  rev_emap <- build_reverse_emap(emap)

  is_true <- logical(J)
  for (j in seq_len(J)) {
    edges_j <- candidate_features[[j]]
    all_present <- TRUE
    for (e in edges_j) {
      uv <- rev_emap[[e]]
      if (is.null(uv) || true_adj[uv[1L], uv[2L]] == 0L) {
        all_present <- FALSE
        break
      }
    }
    is_true[j] <- all_present
  }

  false_positives <- sum(psi == 1L & !is_true)
  fdp             <- false_positives / max(1L, n_rejected)

  n_true         <- sum(is_true)
  true_positives <- sum(psi == 1L & is_true)
  power          <- if (n_true > 0L) true_positives / n_true else 0

  list(
    fdp              = fdp,
    power            = power,
    n_rejected       = n_rejected,
    n_true_positive  = true_positives,
    n_false_positive = false_positives,
    n_true_features  = n_true,
    n_null_features  = sum(!is_true)
  )
}

#' FDP and power from pre-classified feature indices.
#'
#' @param rejected   Integer vector of rejected feature indices.
#' @param h0_indices Integer vector of null hypothesis indices.
#' @param h1_indices Integer vector of alternative hypothesis indices.
#' @return List with `fdp` and `power`.
#' @export
compute_fdp_power_from_indices <- function(rejected, h0_indices, h1_indices) {
  n_rejected <- length(rejected)
  if (n_rejected == 0L) {
    return(list(fdp = 0, power = 0))
  }
  false_positives <- length(intersect(rejected, h0_indices))
  fdp             <- false_positives / max(1L, n_rejected)

  true_positives <- length(intersect(rejected, h1_indices))
  power <- if (length(h1_indices) > 0L) {
    true_positives / length(h1_indices)
  } else {
    0
  }

  list(fdp = fdp, power = power)
}

#' uFDP over a filtration interval.
#'
#' Evaluates
#' \deqn{\text{uFDP} = \sup_{\mu \in [\mu_0, \mu_1]}
#'   \frac{\operatorname{rank}(\hat Z(\mu)) -
#'         \operatorname{rank}(Z(\mu) \cap \hat Z(\mu))}
#'        {\max(1, \operatorname{rank}(\hat Z(\mu)))}}
#' on a discrete grid.
#'
#' @param khan_result  Output of [khan()].
#' @param true_rank_fn Function `mu -> rank(Z(mu))`.
#' @param overlap_fn   Function `mu -> rank(Z(mu) cap hat Z(mu))`.
#' @param mu_grid      Numeric grid of `mu` values.
#' @return List with `ufdp` (scalar) and `fdp_by_mu` (numeric vector).
#' @export
compute_ufdp <- function(khan_result, true_rank_fn, overlap_fn, mu_grid) {
  fdp_by_mu <- numeric(length(mu_grid))
  for (i in seq_along(mu_grid)) {
    mu  <- mu_grid[i]
    res <- khan_result$selected_by_mu(mu)
    if (is.null(res) || res$homology_rank == 0L) {
      fdp_by_mu[i] <- 0
      next
    }
    fdp_by_mu[i] <-
      (res$homology_rank - overlap_fn(mu)) / max(1L, res$homology_rank)
  }
  list(ufdp = max(fdp_by_mu), fdp_by_mu = fdp_by_mu)
}

#' Uniform-power-style power curve for persistent homology.
#'
#' @param khan_result  Output of [khan()].
#' @param true_rank_fn Function `mu -> rank(Z(mu))`.
#' @param d            Number of nodes.
#' @param n            Sample size.
#' @param C_const      Signal strength constant `C` in the shift offset.
#' @param mu_grid      Numeric grid of `mu` values.
#' @return List with `power_by_mu` (numeric vector).
#' @export
compute_uniform_power <- function(khan_result, true_rank_fn, d, n,
                                  C_const = 1, mu_grid) {
  power_by_mu <- numeric(length(mu_grid))
  offset <- C_const * sqrt(log(d) / n)

  for (i in seq_along(mu_grid)) {
    mu  <- mu_grid[i]
    res <- khan_result$selected_by_mu(mu)
    if (is.null(res)) {
      power_by_mu[i] <- 0
      next
    }
    true_rank_shifted <- true_rank_fn(mu + offset)
    if (true_rank_shifted == 0L) {
      power_by_mu[i] <- 1
      next
    }
    power_by_mu[i] <- min(res$homology_rank, true_rank_shifted) / true_rank_shifted
  }
  list(power_by_mu = power_by_mu)
}

#' Summarize FDP and Power across Monte Carlo repetitions.
#'
#' @param results_list List of lists, each with `fdp` and `power` components.
#' @return List with `mean_fdp`, `se_fdp`, `mean_power`, `se_power`, `n_rep`.
#' @export
summarize_repetitions <- function(results_list) {
  fdps   <- vapply(results_list, function(r) r$fdp,   numeric(1L))
  powers <- vapply(results_list, function(r) r$power, numeric(1L))
  n_rep  <- length(results_list)

  list(
    mean_fdp   = mean(fdps),
    se_fdp     = stats::sd(fdps)   / sqrt(n_rep),
    mean_power = mean(powers),
    se_power   = stats::sd(powers) / sqrt(n_rep),
    n_rep      = n_rep
  )
}
