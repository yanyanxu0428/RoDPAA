test_that("dist_to_kernel gives a centred, symmetric, PSD matrix", {
  set.seed(1)
  X <- matrix(rnorm(20 * 5), 20, 5)
  K <- dist_to_kernel(dist(X))
  expect_equal(dim(K), c(20, 20))
  expect_equal(K, t(K))
  expect_equal(unname(rowSums(K)), rep(0, 20), tolerance = 1e-8)
  expect_true(min(eigen(K, symmetric = TRUE, only.values = TRUE)$values) > -1e-8)
})

test_that("for Euclidean distances the kernel equals the centred linear kernel", {
  set.seed(2)
  X <- matrix(rnorm(15 * 4), 15, 4)
  expect_equal(dist_to_kernel(dist(X)), linear_kernel(X), tolerance = 1e-8,
               ignore_attr = TRUE)
})

test_that("negative eigenvalues of non-Euclidean distances are removed", {
  D <- matrix(c(0, 1, 1, 3,
                1, 0, 1, 1,
                1, 1, 0, 1,
                3, 1, 1, 0), 4, 4)   # violates the triangle inequality
  K <- dist_to_kernel(D)
  expect_true(min(eigen(K, symmetric = TRUE, only.values = TRUE)$values) > -1e-10)
})
