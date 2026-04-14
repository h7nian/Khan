# =============================================================================
# bootstrap.R -- Gaussian multiplier bootstrap baselines.
#
# These functions provide the bootstrap-calibrated max-statistic comparison
# baselines used in the paper.  KHAN itself does not require a bootstrap; the
# routines here exist to benchmark against BHq / DGS.
# =============================================================================

#' Gaussian multiplier bootstrap test statistics.
#'
#' For each bootstrap replicate `b = 1..num_B` and each edge `e`:
#' `T^{(b)}_e = (1 / sqrt(n)) sum_i g_i * xi_i(e) / sigma_e`
#' where `g_i ~ N(0, 1)`, `xi_i(e)` are the GGM influence-function summands,
#' and `sigma_e` is the empirical standard deviation of `xi_.(e)`.
#'
#' @param X         n-by-d data matrix.
#' @param theta_hat Estimated precision matrix.
#' @param sigma_hat Estimated (or sample) covariance matrix; currently unused
#'   but retained for interface symmetry with other estimators that may need
#'   it.
#' @param emap      Edge-index map.
#' @param num_B     Number of bootstrap replicates.
#' @return Numeric matrix of dimension `num_B` by `d*(d-1)/2`.
#' @export
bootstrap_test_statistics <- function(X, theta_hat, sigma_hat, emap,
                                      num_B = 1000L) {
  n       <- nrow(X)
  d       <- ncol(X)
  n_edges <- d * (d - 1L) / 2L

  xi <- matrix(0, nrow = n, ncol = n_edges)
  for (u in seq_len(d - 1L)) {
    for (v in (u + 1L):d) {
      idx <- emap[u, v]
      theta_u_x <- X %*% theta_hat[u, ]
      theta_v_x <- X %*% theta_hat[, v]
      xi[, idx] <- theta_u_x * theta_v_x - theta_hat[u, v]
    }
  }

  sigma_vec <- apply(xi, 2L, stats::sd)
  sigma_vec[sigma_vec < 1e-15] <- 1e-15

  boot_stats <- matrix(0, nrow = num_B, ncol = n_edges)
  for (b in seq_len(num_B)) {
    g <- stats::rnorm(n)
    boot_stats[b, ] <- crossprod(g, xi) / (sqrt(n) * sigma_vec)
  }
  boot_stats
}

#' Bootstrap-calibrated feature test via per-feature max statistics.
#'
#' For every feature `F_j` the observed statistic `max_{e in E(F_j)} |T_e|`
#' is compared to its bootstrap distribution to form a p-value, followed by
#' BH or Bonferroni correction across features.
#'
#' @param X              n-by-d data matrix.
#' @param theta_hat      Estimated precision.
#' @param theta_d        Debiased precision.
#' @param var_mat        Variance matrix from [estimate_variance_ggm()].
#' @param sigma_hat      Estimated covariance.
#' @param n              Sample size.
#' @param feature_edges  List of edge-index vectors per feature.
#' @param emap           Edge-index map.
#' @param q              Error-rate target.
#' @param num_B          Number of bootstrap replicates.
#' @param correction     `"bh"` or `"bonferroni"`.
#' @return List with `rejected`, `psi`, `p_values_boot`.
#' @export
bootstrap_feature_test <- function(X, theta_hat, theta_d, var_mat, sigma_hat,
                                   n, feature_edges, emap, q,
                                   num_B = 1000L, correction = "bh") {
  J <- length(feature_edges)
  if (J == 0L) {
    return(list(rejected = integer(0), psi = integer(0),
                p_values_boot = numeric(0)))
  }

  d       <- ncol(X)
  n_edges <- d * (d - 1L) / 2L

  obs_stat <- numeric(n_edges)
  for (u in seq_len(d - 1L)) {
    for (v in (u + 1L):d) {
      idx <- emap[u, v]
      se  <- sqrt(max(var_mat[u, v], 1e-15))
      obs_stat[idx] <- sqrt(n) * abs(theta_d[u, v]) / se
    }
  }

  boot_stats <- bootstrap_test_statistics(X, theta_hat, sigma_hat, emap, num_B)

  p_values_boot <- numeric(J)
  for (j in seq_len(J)) {
    edges_j  <- feature_edges[[j]]
    obs_max  <- max(obs_stat[edges_j])
    boot_max <- apply(abs(boot_stats[, edges_j, drop = FALSE]), 1L, max)
    p_values_boot[j] <- mean(boot_max >= obs_max)
  }

  if (correction == "bonferroni") {
    threshold <- q / J
    psi <- as.integer(p_values_boot < threshold)
  } else {
    bh  <- .benjamini_hochberg(p_values_boot, q, J)
    psi <- bh$psi
  }

  list(
    rejected      = which(psi == 1L),
    psi           = psi,
    p_values_boot = p_values_boot
  )
}

#' Fully simultaneous bootstrap test.
#'
#' Constructs a single global critical value from the bootstrap distribution
#' of `max_j max_{e in E(F_j)} |T_e|`.  This is the most conservative
#' baseline.
#'
#' @inheritParams bootstrap_feature_test
#' @return List with `rejected`, `psi`, `global_critical_value`,
#'   `obs_max_per_feature`.
#' @export
bootstrap_simultaneous_test <- function(X, theta_hat, theta_d, var_mat, sigma_hat,
                                        n, feature_edges, emap, q,
                                        num_B = 1000L) {
  J <- length(feature_edges)
  if (J == 0L) {
    return(list(rejected = integer(0), psi = integer(0),
                global_critical_value = NA_real_,
                obs_max_per_feature = numeric(0)))
  }

  d       <- ncol(X)
  n_edges <- d * (d - 1L) / 2L

  obs_stat <- numeric(n_edges)
  for (u in seq_len(d - 1L)) {
    for (v in (u + 1L):d) {
      idx <- emap[u, v]
      se  <- sqrt(max(var_mat[u, v], 1e-15))
      obs_stat[idx] <- sqrt(n) * abs(theta_d[u, v]) / se
    }
  }

  boot_stats <- bootstrap_test_statistics(X, theta_hat, sigma_hat, emap, num_B)

  obs_max_j <- vapply(feature_edges,
                      function(edges_j) max(obs_stat[edges_j]),
                      numeric(1L))

  boot_global_max <- numeric(num_B)
  for (b in seq_len(num_B)) {
    max_across <- 0
    for (j in seq_len(J)) {
      feat_max <- max(abs(boot_stats[b, feature_edges[[j]]]))
      if (feat_max > max_across) max_across <- feat_max
    }
    boot_global_max[b] <- max_across
  }

  critical_value <- stats::quantile(boot_global_max, probs = 1 - q)
  psi            <- as.integer(obs_max_j > critical_value)

  list(
    rejected              = which(psi == 1L),
    psi                   = psi,
    global_critical_value = critical_value,
    obs_max_per_feature   = obs_max_j
  )
}
