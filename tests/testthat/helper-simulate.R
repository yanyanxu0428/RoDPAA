# Simulated two-view data with k shared latent factors, used by several test files.
simulate_two_views <- function(n = 60, k = 2, p1 = 30, p2 = 40, sd = 0.5) {
  z <- matrix(rnorm(n * k), n, k)
  X1 <- z %*% matrix(rnorm(k * p1), k, p1) + matrix(rnorm(n * p1, sd = sd), n, p1)
  X2 <- z %*% matrix(rnorm(k * p2), k, p2) + matrix(rnorm(n * p2, sd = sd), n, p2)
  list(X1 = X1, X2 = X2, z = z)
}
