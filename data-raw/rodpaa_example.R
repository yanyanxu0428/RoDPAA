## Code to create the `rodpaa_example` dataset (simulated, not real participants).
## 80 participants, two views, two latent factors shared by both views,
## plus one factor specific to each view and independent noise.

set.seed(2026)
n  <- 80
p1 <- 50   # number of features in view 1
p2 <- 40   # number of features in view 2

ids <- sprintf("P%02d", seq_len(n))
shared   <- matrix(rnorm(n * 2), n, 2, dimnames = list(ids, c("Z1", "Z2")))  # shared factors
spec1    <- rnorm(n)                                                         # view-1-specific factor
spec2    <- rnorm(n)                                                         # view-2-specific factor

view1 <- shared %*% matrix(rnorm(2 * p1), 2, p1) + 1.5 * spec1 %o% rnorm(p1) +
  matrix(rnorm(n * p1), n, p1)
view2 <- shared %*% matrix(rnorm(2 * p2), 2, p2) + 1.5 * spec2 %o% rnorm(p2) +
  matrix(rnorm(n * p2), n, p2)
dimnames(view1) <- list(ids, sprintf("feature1_%02d", seq_len(p1)))
dimnames(view2) <- list(ids, sprintf("feature2_%02d", seq_len(p2)))

rodpaa_example <- list(view1 = view1, view2 = view2, shared = shared)
usethis::use_data(rodpaa_example, overwrite = TRUE)
