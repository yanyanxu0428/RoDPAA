set.seed(21)
sim <- simulate_two_views(n = 40)
ids <- sprintf("p%02d", 1:40)
X1 <- sim$X1; X2 <- sim$X2
rownames(X1) <- rownames(X2) <- ids
A1 <- exp(X1); A2 <- exp(X2)                     # positive data, for Bray-Curtis

test_that("data input gives the same fit as precomputed dissimilarities (Euclidean)", {
  set.seed(1); a <- rodpaa(X1, X2, input = "data", n_perm = 100, verbose = FALSE)
  set.seed(1); b <- rodpaa(dist(X1), dist(X2), input = "dissimilarity", n_perm = 100, verbose = FALSE)
  expect_identical(unclass(a), unclass(b))
  expect_equal(rownames(a$shared_factors), ids)  # participants kept, in order
})

test_that("data input gives the same fit as precomputed dissimilarities (Bray-Curtis)", {
  skip_if_not_installed("vegan")
  set.seed(1); a <- rodpaa(A1, A2, input = "data", dissimilarity1 = "bray", dissimilarity2 = "bray",
                           n_perm = 100, verbose = FALSE)
  set.seed(1); b <- rodpaa(vegan::vegdist(A1, "bray"), vegan::vegdist(A2, "bray"),
                           input = "dissimilarity", n_perm = 100, verbose = FALSE)
  expect_equal(a$dimensions, b$dimensions)
  expect_equal(a$p_values, b$p_values)
  expect_equal(a$shared_factors, b$shared_factors, tolerance = 1e-8)
})

test_that("the built-in Bray-Curtis equals vegan::vegdist", {
  skip_if_not_installed("vegan")
  expect_equal(as.matrix(bray_curtis(A1)), as.matrix(vegan::vegdist(A1, "bray")), tolerance = 1e-12)
})

test_that("each view can use its own dissimilarity, and data frames are accepted", {
  set.seed(1); a <- rodpaa(as.data.frame(A1), X2, input = "data", dissimilarity1 = "bray", n_perm = 50, verbose = FALSE)
  set.seed(1); b <- rodpaa(bray_curtis(A1), dist(X2), input = "dissimilarity", n_perm = 50, verbose = FALSE)
  expect_equal(a$shared_factors, b$shared_factors, tolerance = 1e-8)
})

test_that("invalid data are rejected with clear messages", {
  fit_data <- function(...) rodpaa(..., input = "data", verbose = FALSE)
  bad <- X1; bad[1, 1] <- NA
  expect_error(fit_data(bad, X2), "missing or infinite")
  bad[1, 1] <- Inf
  expect_error(fit_data(bad, X2), "missing or infinite")
  expect_error(fit_data(data.frame(id = ids, X1), X2), "non-numeric columns")
  expect_error(fit_data(X1, X2[-1, ]), "different numbers of participants")
  Y2 <- X2; rownames(Y2) <- rev(ids)
  expect_error(fit_data(X1, Y2), "different row names")
  expect_error(fit_data(X1, X2, dissimilarity1 = "bray"), "nonnegative")
  Z <- A1; Z[3, ] <- 0
  expect_error(fit_data(Z, A2, dissimilarity1 = "bray"), "all zero")
  expect_error(fit_data(X1, X2, dissimilarity1 = "jaccard"))
})

test_that("inputs given with the wrong `input` are caught", {
  D1 <- dist(X1); D2 <- dist(X2)
  K1 <- dist_to_kernel(D1); K2 <- dist_to_kernel(D2)
  # dissimilarities passed as data
  expect_error(rodpaa(D1, D2, input = "data", verbose = FALSE), "Use `input = \"dissimilarity\"`", fixed = TRUE)
  w <- capture_warnings(rodpaa(as.matrix(D1), as.matrix(D2), input = "data", n_perm = 20, verbose = FALSE))
  expect_length(w, 2)                                        # one warning per view
  expect_match(w, "looks like a dissimilarity matrix")
  # data or kernels passed as dissimilarities (the default)
  expect_error(rodpaa(X1, X2, verbose = FALSE), "use `input = \"data\"`", fixed = TRUE)
  expect_error(rodpaa(K1, K2, verbose = FALSE), "looks like a kernel")
  # dissimilarities passed as kernels
  expect_error(rodpaa(as.matrix(D1), as.matrix(D2), input = "kernel", verbose = FALSE), "looks like a dissimilarity")
})

test_that("malformed dissimilarities are rejected, including dist objects", {
  D <- as.matrix(dist(X1))
  D_na <- D; D_na[1, 2] <- D_na[2, 1] <- NA
  expect_error(rodpaa(as.dist(D_na), dist(X2), verbose = FALSE), "no missing or infinite")
  D_neg <- D; D_neg[1, 2] <- D_neg[2, 1] <- -1
  expect_error(rodpaa(as.dist(D_neg), dist(X2), verbose = FALSE), "nonnegative")
  D_asym <- D; D_asym[1, 2] <- D_asym[1, 2] + 1
  expect_error(rodpaa(D_asym, dist(X2), verbose = FALSE), "symmetric")
  D_diag <- D; D_diag[1, 1] <- 1
  expect_error(rodpaa(D_diag, dist(X2), verbose = FALSE), "zero diagonal")
  expect_error(rodpaa(matrix("a", 40, 40), dist(X2), verbose = FALSE), "numeric")
  expect_error(rodpaa(dist(X1[1, , drop = FALSE]), dist(X2[1, , drop = FALSE]), verbose = FALSE),
               "at least two participants")
})

test_that("one-participant data input is rejected", {
  expect_error(rodpaa(X1[1, , drop = FALSE], X2[1, , drop = FALSE], input = "data", verbose = FALSE),
               "at least two participants")
})

test_that("legitimate square data matrices are accepted", {
  set.seed(31)
  S1 <- matrix(rnorm(30 * 30), 30, 30)            # square, but not a dissimilarity
  S2 <- matrix(rnorm(30 * 30), 30, 30)
  expect_no_warning(fit <- rodpaa(S1, S2, input = "data", n_perm = 50, verbose = FALSE))
  expect_s3_class(fit, "rodpaa")
  # square, symmetric, zero diagonal data: a warning, but the analysis runs as requested
  Z1 <- (S1 + t(S1)) / 2; diag(Z1) <- 0; Z1 <- abs(Z1)
  set.seed(1)
  expect_warning(a <- rodpaa(Z1, S2, input = "data", n_perm = 50, verbose = FALSE), "looks like a dissimilarity")
  set.seed(1)
  b <- rodpaa(dist(Z1), dist(S2), input = "dissimilarity", n_perm = 50, verbose = FALSE)
  expect_identical(unclass(a), unclass(b))
})

test_that("looks_like_dissimilarity() never errors on odd input", {
  expect_false(looks_like_dissimilarity(matrix("a", 3, 3)))
  expect_false(looks_like_dissimilarity(matrix(NA_real_, 3, 3)))
  expect_false(looks_like_dissimilarity(matrix(numeric(0), 0, 0)))
  expect_false(looks_like_dissimilarity(matrix(0, 2, 3)))
  expect_false(looks_like_dissimilarity(matrix(0, 1, 1)))
  expect_false(looks_like_dissimilarity(list(1, 2)))
  expect_false(looks_like_dissimilarity(NULL))
  expect_true(looks_like_dissimilarity(as.matrix(dist(X1))))
})
