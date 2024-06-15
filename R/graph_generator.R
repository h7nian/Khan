#' Title
#'
#' @param n
#' @param d
#' @param graph
#' @param v
#' @param u
#' @param p
#' @param prob
#' @param vis
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
graph_generator <- function(n = 200, d = 50, graph = "SBM", v = NULL, u = NULL,
                            p = 10, prob = NULL, vis = FALSE, verbose = TRUE) {
  gcinfo(FALSE)
  if (verbose) {
    cat(
      "Generating data from the multivariate normal distribution with the",
      graph, "graph structure...."
    )
  }

  theta <- omega <- matrix(0, d, d)
  if (d %% p > 0) {
    cat("p should be a factor of d", "\n")
    return(0)
  }
  g <- d / p
  grp.idx <- rep(1:g, rep(p, g))
  if (graph == "SBM") {
    O <- matrix(0.0, g, g)
    diag(O) <- prob

    theta <- 1 * (matrix(stats::runif(d * d), nrow = d) < O[grp.idx, grp.idx])
    diag(theta) <- 0
    theta[lower.tri(theta)] <- t(theta)[lower.tri(theta)]
  }

  e <- v * stats::runif(d * (d - 1) / 2, min = 0.9, max = 1) * sign(stats::runif(d * (d - 1) / 2, min = -1, max = 1))
  omega[upper.tri(omega)] <- e
  omega <- omega + t(omega)
  omega[theta == 0] <- 0

  diag(omega) <- diag(omega) + abs(min(eigen(omega)$values)) + u
  omega <- omega * 10

  min(eigen(omega)$values)
  max(eigen(omega)$values)

  sigma <- solve(omega)
  x <- MASS::mvrnorm(n, rep(0, d), sigma)
  sigmahat <- stats::cov(x)
  if (vis == TRUE) {
    fullfig <- graphics::par(mfrow = c(2, 2), pty = "s", omi = c(
      0.3,
      0.3, 0.3, 0.3
    ), mai = c(0.3, 0.3, 0.3, 0.3))
    fullfig[1] <- graphics::image(theta, col = grDevices::gray.colors(256), main = "Adjacency Matrix")
    fullfig[2] <- graphics::image(sigma, col = grDevices::gray.colors(256), main = "Covariance Matrix")
    g <- igraph::graph.adjacency(theta, mode = "undirected", diag = FALSE)
    layout.grid <- igraph::layout.fruchterman.reingold(g)
    fullfig[3] <- plot(g,
      layout = layout.grid, edge.color = "gray50",
      vertex.color = "red", vertex.size = 3, vertex.label = NA,
      main = "Graph Pattern"
    )
    fullfig[4] <- graphics::image(sigmahat,
      col = grDevices::gray.colors(256),
      main = "Empirical Matrix"
    )
    rm(fullfig, g, layout.grid)
    gc()
  }
  if (verbose) {
    cat("done.\n")
  }
  rm(vis, verbose)
  gc()
  # sigma: the true covariance matrix
  # omega: the ture percision matrix
  # theta: the support set of omega
  sim <- list(
    data = x, sigma = sigma, sigmahat = sigmahat,
    omega = omega, theta = theta,
    sparsity = sum(theta) / (d * (d - 1)), graph.type = graph
  )
  class(sim) <- "sim"
  return(sim)
}
