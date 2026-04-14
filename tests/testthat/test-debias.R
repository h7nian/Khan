# Regression test for graph_debias.R (#C3): the legacy implementation had
# an operator-precedence bug AND subtracted Omega[j, k] instead of I(j == k).

test_that("debias_glasso recovers Theta* when Sigma = Theta^{-1}", {
  set.seed(1L)
  d     <- 6L
  A     <- matrix(stats::runif(d * d, -0.2, 0.2), d, d)
  theta <- A + t(A) + diag(d) * (d + 1)
  sigma <- solve(theta)

  theta_d <- debias_glasso(sigma, theta)
  # When Theta Sigma = I exactly, the correction term is zero and
  # theta_d should equal theta.
  expect_equal(theta_d, theta, tolerance = 1e-8)
})

test_that("debias_glasso matches a naive element-wise formula", {
  set.seed(2L)
  d         <- 5L
  X         <- matrix(stats::rnorm(200 * d), 200, d)
  sigma_hat <- crossprod(X) / 200
  theta_hat <- solve(sigma_hat + diag(d) * 0.1)

  theta_d_vec <- debias_glasso(sigma_hat, theta_hat)

  # Naive element-wise reference: direct transcription of equation (2.7).
  # denom = theta_hat[u, ] %*% sigma_hat[, u] (scalar).
  # numer = theta_hat[u, ] %*% (sigma_hat %*% theta_hat[, v] - e_v).
  naive <- matrix(0, d, d)
  for (u in seq_len(d)) {
    denom <- sum(theta_hat[u, ] * sigma_hat[, u])
    for (v in seq_len(d)) {
      resid <- sigma_hat %*% theta_hat[, v]
      resid[v] <- resid[v] - 1
      numer <- sum(theta_hat[u, ] * resid)
      naive[u, v] <- theta_hat[u, v] - numer / denom
    }
  }
  expect_equal(theta_d_vec, naive, tolerance = 1e-10)
})

test_that("estimate_variance_ggm is symmetric and non-negative off-diagonal", {
  d   <- 4L
  m   <- matrix(c(4, 1, 0, 0.5,
                  1, 3, 0.2, 0,
                  0, 0.2, 5, 0.1,
                  0.5, 0, 0.1, 2), d, d)
  v   <- estimate_variance_ggm(m)
  expect_equal(v, t(v))
  expect_true(all(v[upper.tri(v)] >= 0))
})
