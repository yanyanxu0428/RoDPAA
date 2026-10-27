# all plots are drawn to a null device
plot_ok <- function(...) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  plot(...)
}

set.seed(5)
sim <- simulate_two_views()
fit <- rodpaa(dist(sim$X1), dist(sim$X2), n_perm = 100, verbose = FALSE)   # k_hat = 2

test_that("all three plot types work and return the fit unchanged", {
  for (tp in c("intensity", "scores", "pvalue")) {
    expect_invisible(out <- plot_ok(fit, type = tp))
    expect_identical(out, fit)
  }
  expect_invisible(plot_ok(fit))                                   # default = intensity
  expect_invisible(plot_ok(fit, type = "scores", group = rep(c("a", "b"), 30)))
})

test_that("plotting does not change the stored P-values, including zeros", {
  expect_true(any(fit$p_adj == 0))
  before <- fit$permutation$table
  plot_ok(fit, type = "pvalue")
  expect_identical(fit$permutation$table, before)
})

test_that("plots work when k_hat = 0", {
  set.seed(11)
  X1 <- matrix(rnorm(50 * 20), 50, 20)
  X2 <- matrix(rnorm(50 * 20), 50, 20)
  f0 <- rodpaa(dist(X1), dist(X2), n_perm = 100, verbose = FALSE)
  expect_equal(f0$k_hat, 0)
  expect_invisible(plot_ok(f0, type = "intensity"))
  expect_invisible(plot_ok(f0, type = "pvalue"))
  expect_error(plot_ok(f0, type = "scores"), "k_extract = 2")
  # exploratory directions can be plotted after k_extract
  f0x <- rodpaa(dist(X1), dist(X2), n_perm = 100, k_extract = 2, verbose = FALSE)
  expect_invisible(plot_ok(f0x, type = "scores"))
})

test_that("plots without permutation results", {
  fk <- rodpaa(dist(sim$X1), dist(sim$X2), k = 2, verbose = FALSE)
  expect_null(fk$permutation)
  expect_invisible(plot_ok(fk, type = "intensity"))
  expect_invisible(plot_ok(fk, type = "scores"))
  expect_error(plot_ok(fk, type = "pvalue"), "no permutation results")
})

test_that("only one displayed component", {
  expect_invisible(plot_ok(fit, type = "intensity", n_show = 1))
  expect_invisible(plot_ok(fit, type = "pvalue", n_show = 1))
})

test_that("invalid plotting arguments give clear errors", {
  expect_error(plot_ok(fit, type = "heatmap"))
  expect_error(plot_ok(fit, n_show = 0), "positive whole number")
  expect_error(plot_ok(fit, n_show = 1.5), "positive whole number")
  expect_error(plot_ok(fit, type = "scores", factors = c(1, 1)), "two different")
  expect_error(plot_ok(fit, type = "scores", factors = c(1, 5)), "between 1 and 2")
  expect_error(plot_ok(fit, type = "scores", group = c("a", "b")), "one non-missing label")
  expect_error(plot_ok(fit, type = "scores", group = c(NA, rep("a", 59))), "one non-missing label")
})
