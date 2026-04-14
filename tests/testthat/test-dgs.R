# Regression tests for DGS (#C5: the legacy find_homo1 iterated in input
# order and never assigned ell copies).

test_that("dgs selects a triangle's closing edge at q = 0.1", {
  d    <- 3L
  emap <- idx_map(d)
  p_edges <- c(1e-6, 1e-6, 1e-6)
  out     <- dgs(p_edges, emap, d, q = 0.1, k = 1L)
  expect_equal(out$homology_rank, 1L)
  expect_equal(length(out$selected_edges), 3L)
})

test_that("dgs fast path matches generic path for k = 1", {
  skip_if_not_installed("igraph")
  set.seed(3L)
  d    <- 5L
  emap <- idx_map(d)
  p_edges <- stats::runif(10L)
  p_edges[c(1, 2, 3)] <- 1e-4
  fast    <- dgs(p_edges, emap, d, q = 0.1, k = 1L)
  generic <- Khan:::.dgs_generic(p_edges, emap, d, q = 0.1, k = 1L)
  expect_equal(fast$homology_rank, generic$homology_rank)
  expect_equal(sort(fast$selected_edges), sort(generic$selected_edges))
})

test_that("dgs returns empty result when no edge passes the screening", {
  d    <- 4L
  emap <- idx_map(d)
  out  <- dgs(rep(0.9, 6L), emap, d, q = 0.1, k = 1L)
  expect_equal(out$homology_rank, 0L)
  expect_length(out$selected_edges, 0L)
})
