# End-to-end smoke test for the KHAN pipeline on a small GGM.

test_that("khan runs on a triangle and produces a non-empty barcode", {
  d       <- 3L
  emap    <- idx_map(d)
  w_hat   <- c(1.2, 1.1, 1.0)
  sigma_e <- c(1.0, 1.0, 1.0)
  n       <- 400L
  res <- khan(w_hat, sigma_e, n, emap, d, q = 0.1,
              mu_range = c(0, 2), scenario = "a", k = 1L)
  expect_true(is.list(res))
  expect_true(is.function(res$selected_by_mu))
  expect_true(length(res$change_points) >= 1L)
})

test_that("khan returns an empty barcode on noise", {
  d       <- 4L
  emap    <- idx_map(d)
  w_hat   <- stats::rnorm(6L, sd = 0.01)
  sigma_e <- rep(1, 6L)
  res <- khan(w_hat, sigma_e, n = 50L, emap, d, q = 0.1,
              mu_range = c(0, 1), scenario = "a", k = 1L)
  expect_equal(nrow(res$barcode), 0L)
})
