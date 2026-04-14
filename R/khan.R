# =============================================================================
# khan.R -- Algorithm 3: KHAN (k-dimensional persistent Homology Adaptive
#           selectioN).
#
# Iteratively locates change points on the filtration line and applies DGS at
# each change point, providing uFDR-controlled persistent homology selection
# (Section 3.2 of the paper).
# =============================================================================

#' KHAN algorithm for persistent homology selection (Algorithm 3).
#'
#' @param w_hat    Edge-weight estimates (numeric vector, length `n_edges`).
#' @param sigma_e  Edge standard deviations (numeric vector).
#' @param n        Sample size.
#' @param emap     Edge-index map.
#' @param d        Number of nodes.
#' @param q        FDR level.
#' @param mu_range Numeric length-2 vector `c(mu_0, mu_1)`; use `c(0, Inf)`
#'   for unbounded tails.
#' @param scenario `"a"` (two-sided / GGM) or `"b"` (one-sided / Ising).
#' @param k        Homology dimension (1 or 2).
#' @return List with:
#'   \describe{
#'     \item{`barcode`}{Data frame with columns `generator_id`, `birth_mu`,
#'       `death_mu`.}
#'     \item{`change_points`}{Numeric vector of change points `mu^{(t)}`.}
#'     \item{`selected_by_mu`}{Function `mu -> dgs()` result at the last
#'       change point `<= mu`.}
#'     \item{`results_at_changes`}{List of DGS results, one per change point.}
#'   }
#' @export
khan <- function(w_hat, sigma_e, n, emap, d, q,
                 mu_range = c(0, Inf), scenario = "a", k = 1L) {
  mu_0 <- mu_range[1L]
  mu_1 <- mu_range[2L]

  # bar{J} = rank(Z_k(K_d)).  For k = 1 this is (d-1)(d-2)/2 by Proposition E.1;
  # the same formula with k substituted is used as a uniform upper bound.
  bar_J <- as.integer((d - k) * (d - k - 1L) / 2L)
  if (bar_J <= 0L) bar_J <- 1L

  rev_emap <- build_reverse_emap(emap)

  p_mu0 <- pvalues_at_mu(w_hat, sigma_e, n, mu_0, scenario)
  dgs_0 <- dgs(p_mu0, emap, d, q, k, rev_emap = rev_emap)

  change_points      <- mu_0
  results_at_changes <- list(dgs_0)
  barcode_entries    <- list()

  current_edges <- dgs_0$selected_edges
  current_rank  <- dgs_0$homology_rank

  t        <- 0L
  max_iter <- 500L
  prev_mu  <- -Inf

  while (current_rank > 0L && t < max_iter) {
    alpha_dynamic <- if (scenario == "a") {
      q * current_rank / (2 * bar_J)
    } else {
      q * current_rank / bar_J
    }

    mu_next <- .find_next_change_point(
      w_hat, sigma_e, n, current_edges, scenario, alpha_dynamic
    )

    if (is.infinite(mu_next) || mu_next > mu_1) {
      for (i in seq_len(current_rank)) {
        barcode_entries <- c(barcode_entries, list(data.frame(
          generator_id = length(barcode_entries) + 1L,
          birth_mu     = change_points[length(change_points)],
          death_mu     = mu_1
        )))
      }
      break
    }

    # Guard against non-advancing mu.
    if (mu_next <= prev_mu + 1e-12) mu_next <- prev_mu + 1e-6
    prev_mu <- mu_next

    p_new   <- pvalues_at_mu(w_hat, sigma_e, n, mu_next, scenario)
    dgs_new <- dgs(p_new, emap, d, q, k, rev_emap = rev_emap)

    rank_diff <- current_rank - dgs_new$homology_rank
    if (rank_diff > 0L) {
      for (i in seq_len(rank_diff)) {
        barcode_entries <- c(barcode_entries, list(data.frame(
          generator_id = length(barcode_entries) + 1L,
          birth_mu     = change_points[length(change_points)],
          death_mu     = mu_next
        )))
      }
    }

    change_points      <- c(change_points, mu_next)
    results_at_changes <- c(results_at_changes, list(dgs_new))
    current_edges      <- dgs_new$selected_edges
    current_rank       <- dgs_new$homology_rank
    t <- t + 1L
  }

  barcode <- if (length(barcode_entries) > 0L) {
    do.call(rbind, barcode_entries)
  } else {
    data.frame(generator_id = integer(0),
               birth_mu     = numeric(0),
               death_mu     = numeric(0))
  }

  cp_vec   <- change_points
  res_list <- results_at_changes
  selected_by_mu <- function(mu) {
    if (mu < mu_0) return(NULL)
    idx <- max(which(cp_vec <= mu))
    res_list[[idx]]
  }

  list(
    barcode            = barcode,
    change_points      = change_points,
    selected_by_mu     = selected_by_mu,
    results_at_changes = results_at_changes
  )
}

# Internal: next filtration change point from equation (3.1).
.find_next_change_point <- function(w_hat, sigma_e, n, current_edges,
                                    scenario, alpha) {
  if (length(current_edges) == 0L) return(Inf)

  z_alpha <- stats::qnorm(1 - alpha)
  offset  <- z_alpha * sigma_e[current_edges] / sqrt(n)

  candidates <- if (scenario == "a") {
    abs(w_hat[current_edges]) - offset
  } else {
    w_hat[current_edges] - offset
  }

  valid <- candidates[candidates > 0]
  if (length(valid) == 0L) return(Inf)
  min(valid)
}
