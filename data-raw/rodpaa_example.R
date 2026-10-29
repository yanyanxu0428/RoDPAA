## Code to create the `rodpaa_example` dataset (simulated, not real participants).
## Latent-factor model of the RoDPAA manuscript,
##   X_v = F_J W_v + F_{I,v} U_v + E_v,  v = 1, 2,
## with n = 80, p1 = 50, p2 = 40, two shared factors (k0 = 2), one additional
## factor per view (q1 = q2 = 1), sigma_indiv = 1.5 and sigma_noise = 1.

set.seed(2026)
n  <- 80
p1 <- 50   # number of features in view 1
p2 <- 40   # number of features in view 2
sigma_indiv <- 1.5
sigma_noise <- 1

ids  <- sprintf("P%02d", seq_len(n))
F_J  <- matrix(rnorm(n * 2), n, 2, dimnames = list(ids, c("F_J1", "F_J2")))  # shared factors
F_J  <- sweep(F_J, 2, colMeans(F_J))
F_I1 <- rnorm(n)                                                            # additional factor, view 1
F_I2 <- rnorm(n)                                                            # additional factor, view 2

W1 <- matrix(rnorm(2 * p1), 2, p1)
U1 <- matrix(rnorm(p1, sd = sigma_indiv), 1, p1)
E1 <- matrix(rnorm(n * p1, sd = sigma_noise), n, p1)
W2 <- matrix(rnorm(2 * p2), 2, p2)
U2 <- matrix(rnorm(p2, sd = sigma_indiv), 1, p2)
E2 <- matrix(rnorm(n * p2, sd = sigma_noise), n, p2)

view1 <- F_J %*% W1 + F_I1 %*% U1 + E1
view2 <- F_J %*% W2 + F_I2 %*% U2 + E2
dimnames(view1) <- list(ids, sprintf("feature1_%02d", seq_len(p1)))
dimnames(view2) <- list(ids, sprintf("feature2_%02d", seq_len(p2)))

rodpaa_example <- list(view1 = view1, view2 = view2, shared = F_J)
usethis::use_data(rodpaa_example, overwrite = TRUE)
