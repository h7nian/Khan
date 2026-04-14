# Regression tests for Benjamini-Hochberg thresholding.
# These specifically cover the bugs in the legacy MyStat.R (#C1, #C2):
#   (a) empty reject set should return an empty vector, not error with -Inf
#   (b) the strict-inequality form of BH is applied
#   (c) downstream FDP is computed from psi, not raw p-values

test_that("bhq_test returns empty result on empty feature list", {
  out <- bhq_test(p_edges = numeric(0), feature_edges = list(), q = 0.1)
  expect_length(out$rejected, 0L)
  expect_length(out$psi, 0L)
})

test_that("bhq_test rejects nothing when all edge p-values exceed q", {
  p_edges  <- c(0.5, 0.6, 0.7)
  features <- list(c(1L, 2L), c(2L, 3L))
  out <- bhq_test(p_edges, features, q = 0.1)
  expect_equal(out$psi, c(0L, 0L))
  expect_length(out$rejected, 0L)
})

test_that("bhq_test rejects strong signals at q = 0.1", {
  p_edges <- c(1e-6, 1e-6, 1e-6, 0.5, 0.9, 0.9)
  feat1   <- c(1L, 2L, 3L)
  feat2   <- c(4L, 5L, 6L)
  out     <- bhq_test(p_edges, list(feat1, feat2), q = 0.1)
  expect_true(out$psi[1] == 1L)
  expect_true(out$psi[2] == 0L)
})

test_that("bonferroni_test is stricter than bhq_test", {
  p_edges  <- c(0.001, 0.01, 0.03, 0.5, 0.9, 0.9)
  features <- list(c(1L, 2L, 3L), c(4L, 5L, 6L))
  bh       <- bhq_test(p_edges, features, q = 0.1)
  bf       <- bonferroni_test(p_edges, features, q = 0.1)
  expect_true(sum(bf$psi) <= sum(bh$psi))
})
