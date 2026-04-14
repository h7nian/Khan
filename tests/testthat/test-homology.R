# Regression tests for compute_homology_rank (#C5, #C6).

test_that("H1 rank of a triangle is 1", {
  d    <- 3L
  emap <- idx_map(d)
  expect_equal(compute_homology_rank(c(1L, 2L, 3L), d, emap, k = 1L), 1L)
})

test_that("H1 rank of a tree is 0", {
  d    <- 4L
  emap <- idx_map(d)
  tree_edges <- c(emap[1, 2], emap[2, 3], emap[3, 4])
  expect_equal(compute_homology_rank(tree_edges, d, emap, k = 1L), 0L)
})

test_that("H1 rank of two disjoint triangles is 2", {
  d    <- 6L
  emap <- idx_map(d)
  edges <- c(emap[1, 2], emap[2, 3], emap[1, 3],
             emap[4, 5], emap[5, 6], emap[4, 6])
  expect_equal(compute_homology_rank(edges, d, emap, k = 1L), 2L)
})

test_that("rank_increase reflects Union-Find structure", {
  d    <- 4L
  emap <- idx_map(d)
  e12  <- emap[1, 2]
  e23  <- emap[2, 3]
  e13  <- emap[1, 3]  # closes triangle with 1-2-3
  e34  <- emap[3, 4]  # bridge

  expect_equal(rank_increase(c(e12, e23), e13, d, emap, k = 1L), 1L)
  expect_equal(rank_increase(c(e12, e23, e13), e34, d, emap, k = 1L), 0L)
})
