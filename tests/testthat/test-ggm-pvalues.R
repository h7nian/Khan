# Regression tests for compute_ggm_pvalues (#C4: legacy variance normalisation
# divided by n^2 instead of n; #C9: emap was pulled from lexical scope).

test_that("compute_ggm_pvalues requires explicit emap", {
  theta_d <- diag(3)
  var_mat <- matrix(1, 3, 3)
  expect_error(compute_ggm_pvalues(theta_d, var_mat, n = 10L,
                                   scenario = "a", emap = NULL),
               "emap")
})

test_that("compute_ggm_pvalues returns values in [0, 1]", {
  set.seed(4L)
  d       <- 4L
  emap    <- idx_map(d)
  theta_d <- matrix(stats::rnorm(d * d, sd = 0.1), d, d)
  theta_d <- (theta_d + t(theta_d)) / 2
  diag(theta_d) <- abs(diag(theta_d)) + 1
  var_mat <- estimate_variance_ggm(theta_d)

  p <- compute_ggm_pvalues(theta_d, var_mat, n = 100L,
                           scenario = "a", emap = emap)
  expect_length(p, d * (d - 1L) / 2L)
  expect_true(all(p >= 0 & p <= 1))
})
