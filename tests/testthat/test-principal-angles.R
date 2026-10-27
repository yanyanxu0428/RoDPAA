random_basis <- function(n, r) qr.Q(qr(matrix(rnorm(n * r), n, r)))

test_that("identical subspaces have shared intensities of 1", {
  set.seed(1)
  O <- random_basis(30, 4)
  R <- random_basis(4, 4)                      # a different basis of the same subspace
  a <- principal_angles(O, O %*% R)
  expect_equal(a$sigma, rep(1, 4), tolerance = 1e-10)
})

test_that("orthogonal subspaces have shared intensities of 0", {
  set.seed(2)
  Q <- random_basis(30, 6)
  a <- principal_angles(Q[, 1:3], Q[, 4:6])
  expect_equal(a$sigma, rep(0, 3), tolerance = 1e-10)
})

test_that("a known angle is recovered as its cosine", {
  theta <- pi / 3
  e <- diag(10)
  O1 <- e[, 1, drop = FALSE]
  O2 <- cos(theta) * e[, 1, drop = FALSE] + sin(theta) * e[, 2, drop = FALSE]
  expect_equal(principal_angles(O1, O2)$sigma, cos(theta))
})

test_that("shared intensities do not depend on the choice of basis", {
  set.seed(3)
  O1 <- random_basis(40, 5)
  O2 <- random_basis(40, 3)
  a <- principal_angles(O1, O2)
  b <- principal_angles(O1 %*% random_basis(5, 5), O2 %*% random_basis(3, 3))
  expect_equal(a$sigma, b$sigma, tolerance = 1e-10)
  expect_true(all(a$sigma >= 0 & a$sigma <= 1 + 1e-12))
  expect_length(a$sigma, 3)                    # min(r1, r2)
})

test_that("paired shared directions have inner product equal to the shared intensity", {
  set.seed(4)
  n <- 40
  O1 <- random_basis(n, 4)
  O2 <- random_basis(n, 5)
  a <- principal_angles(O1, O2)
  S1 <- list(r = 4, O = O1)
  S2 <- list(r = 5, O = O2)
  K <- diag(n)
  sf <- construct_shared_factors(S1, S2, a, K, K, k = 3)
  expect_equal(colSums(sf$O_J1 * sf$O_J2), a$sigma[1:3], tolerance = 1e-10)
  expect_equal(sf$O_J, (sf$O_J1 + sf$O_J2) / 2)
  expect_equal(ncol(sf$O_I1), 1)               # r1 - k view-specific directions
  expect_equal(ncol(sf$O_I2), 2)               # r2 - k
})
