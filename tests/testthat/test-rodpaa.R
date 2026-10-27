test_that("rodpaa recovers a known shared dimension", {
  set.seed(10)
  sim <- simulate_two_views()
  fit <- rodpaa(dist(sim$X1), dist(sim$X2), n_perm = 200, verbose = FALSE)
  expect_s3_class(fit, "rodpaa")
  expect_equal(fit$k_hat, 2)
  expect_equal(dim(fit$shared_factors), c(60, 2))
  # the shared factors span the true latent factors
  r2 <- summary(lm(fit$shared_factors ~ sim$z))
  expect_true(all(vapply(r2, function(s) s$r.squared, numeric(1)) > 0.9))
})

test_that("unrelated views give no shared dimension", {
  set.seed(11)
  X1 <- matrix(rnorm(50 * 20), 50, 20)
  X2 <- matrix(rnorm(50 * 20), 50, 20)
  fit <- rodpaa(dist(X1), dist(X2), n_perm = 200, verbose = FALSE)
  expect_equal(fit$k_hat, 0)
  expect_equal(ncol(fit$shared_factors), 0)
})

test_that("results are reproducible with set.seed", {
  sim <- simulate_two_views()
  set.seed(3); a <- rodpaa(dist(sim$X1), dist(sim$X2), n_perm = 100, verbose = FALSE)
  set.seed(3); b <- rodpaa(dist(sim$X1), dist(sim$X2), n_perm = 100, verbose = FALSE)
  expect_identical(a$permutation$table, b$permutation$table)
})

test_that("kernel input gives the same result as dissimilarity input", {
  sim <- simulate_two_views()
  set.seed(4); a <- rodpaa(dist(sim$X1), dist(sim$X2), input = "dissimilarity", n_perm = 100, verbose = FALSE)
  K1 <- dist_to_kernel(dist(sim$X1)); K2 <- dist_to_kernel(dist(sim$X2))
  set.seed(4); b <- rodpaa(K1, K2, input = "kernel", n_perm = 100, verbose = FALSE)
  expect_equal(a$shared_factors, b$shared_factors, ignore_attr = TRUE)
  expect_equal(a$k_hat, b$k_hat)
})

test_that("user-supplied k skips the permutation test", {
  sim <- simulate_two_views()
  fit <- rodpaa(dist(sim$X1), dist(sim$X2), k = 1, verbose = FALSE)
  expect_null(fit$permutation)
  expect_equal(fit$k_hat, 1)
  expect_equal(ncol(fit$shared_factors), 1)
  r <- min(fit$dimensions$r1, fit$dimensions$r2)
  expect_error(rodpaa(dist(sim$X1), dist(sim$X2), k = r + 1, verbose = FALSE), "signal ranks")
})

test_that("k_extract extracts directions without changing k_hat", {
  set.seed(12)
  X1 <- matrix(rnorm(40 * 10), 40, 10)
  X2 <- matrix(rnorm(40 * 10), 40, 10)
  fit <- rodpaa(dist(X1), dist(X2), n_perm = 100, k_extract = 2, verbose = FALSE)
  expect_equal(fit$k_hat, 0)
  expect_equal(fit$dimensions$k_extracted, 2)
  expect_equal(ncol(fit$shared_factors), 2)
})

test_that("mismatched participants are rejected", {
  X <- matrix(rnorm(20 * 3), 20, 3, dimnames = list(paste0("p", 1:20), NULL))
  Y <- X[20:1, ]
  expect_error(rodpaa(dist(X), dist(Y), verbose = FALSE), "row names")
})

test_that("the example data give two shared factors", {
  set.seed(1)
  fit <- rodpaa(dist(rodpaa_example$view1), dist(rodpaa_example$view2), verbose = FALSE)
  expect_equal(fit$k_hat, 2)
  r2 <- summary(lm(fit$shared_factors ~ rodpaa_example$shared))
  expect_true(all(vapply(r2, function(s) s$r.squared, numeric(1)) > 0.9))
})

test_that("kernel input is centered automatically and keeps participant names", {
  set.seed(6)
  sim <- simulate_two_views()
  ids <- sprintf("p%02d", 1:60)
  rownames(sim$X1) <- rownames(sim$X2) <- ids
  K1 <- dist_to_kernel(dist(sim$X1)); K2 <- dist_to_kernel(dist(sim$X2))
  set.seed(7); a <- rodpaa(K1, K2, input = "kernel", n_perm = 100, verbose = FALSE)
  # an already-centered kernel gives the same result after centering again
  set.seed(7); b <- rodpaa(center_kernel(K1), center_kernel(K2), input = "kernel", n_perm = 100, verbose = FALSE)
  expect_equal(a$dimensions, b$dimensions)
  expect_equal(a$shared_factors, b$shared_factors, tolerance = 1e-8)
  # an uncentered kernel (constant added) gives the same result as the centered one
  set.seed(7); u <- rodpaa(K1 + 5, K2 + 5, input = "kernel", n_perm = 100, verbose = FALSE)
  expect_equal(u$dimensions, a$dimensions)
  expect_equal(u$shared_factors, a$shared_factors, tolerance = 1e-8)
  expect_equal(rownames(a$shared_factors), ids)
})

test_that("non-symmetric or non-square kernels are rejected", {
  K <- dist_to_kernel(dist(matrix(rnorm(40), 10)))
  K_bad <- K; K_bad[1, 2] <- K_bad[1, 2] + 1
  expect_error(rodpaa(K_bad, K, input = "kernel", verbose = FALSE), "`x1` is not symmetric")
  expect_error(rodpaa(K, K[, 1:9], input = "kernel", verbose = FALSE), "`x2` must be a square")
})
