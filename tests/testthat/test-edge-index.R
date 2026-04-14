test_that("idx_map is symmetric and covers every edge", {
  d <- 5L
  emap <- idx_map(d)
  expect_equal(dim(emap), c(d, d))
  expect_true(all(emap == t(emap)))
  expect_equal(diag(emap), rep(0L, d))
  expect_equal(sort(unique(emap[upper.tri(emap)])), 1:(d * (d - 1L) / 2L))
})

test_that("edge_to_nodes inverts nodes_to_edge", {
  d    <- 7L
  emap <- idx_map(d)
  all_edges <- 1:(d * (d - 1L) / 2L)
  uv   <- edge_to_nodes(all_edges, d)
  back <- nodes_to_edge(uv[, 1], uv[, 2], emap)
  expect_equal(back, all_edges)
})

test_that("build_reverse_emap agrees with idx_map", {
  d    <- 4L
  emap <- idx_map(d)
  rev_emap <- build_reverse_emap(emap)
  expect_length(rev_emap, d * (d - 1L) / 2L)
  for (e in seq_along(rev_emap)) {
    uv <- rev_emap[[e]]
    expect_equal(emap[uv[1], uv[2]], e)
  }
})
