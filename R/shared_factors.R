#' Construct the shared factors (Step 5)
#'
#' For the first \eqn{k} principal-angle pairs, forms the paired shared
#' directions \eqn{f_{1j} = O_1 p_j} and \eqn{f_{2j} = O_2 q_j} and averages
#' them into the shared factor \eqn{\mathrm{SF}_j = (f_{1j} + f_{2j}) / 2}.
#' The remaining directions of each view are its view-specific directions.
#' Because \eqn{f_{1j}^\top f_{2j} = \sigma_j > 0}, each pair is consistently
#' oriented before averaging.
#'
#' @param S1,S2 Signal subspaces from `estimate_signal_subspace()`.
#' @param angles Output of `principal_angles()`.
#' @param K1,K2 The two centred similarity matrices (for the reconstructions).
#' @param k Number of shared directions to construct.
#' @return A list with the shared directions `O_J1`, `O_J2`, the shared factors
#'   `O_J`, the view-specific directions `O_I1`, `O_I2`, their eigenvalue
#'   matrices, and the joint and individual reconstructions of `K1` and `K2`.
#' @keywords internal
#' @noRd
construct_shared_factors <- function(S1, S2, angles, K1, K2, k) {
  n <- nrow(K1)
  r1 <- S1$r
  r2 <- S2$r

  # paired shared directions and their average
  if (k > 0) {
    O_J1 <- S1$O %*% angles$P[, 1:k, drop = FALSE]
    O_J2 <- S2$O %*% angles$Q[, 1:k, drop = FALSE]
    O_J <- (O_J1 + O_J2) / 2
  } else {
    O_J1 <- matrix(0, n, 0)
    O_J2 <- matrix(0, n, 0)
    O_J <- matrix(0, n, 0)
  }

  # view-specific directions
  O_I1 <- if (r1 > k) S1$O %*% angles$P[, (k + 1):r1, drop = FALSE] else matrix(0, n, 0)
  O_I2 <- if (r2 > k) S2$O %*% angles$Q[, (k + 1):r2, drop = FALSE] else matrix(0, n, 0)

  # eigenvalue matrices and reconstructions
  Lambda_1J <- if (k > 0) t(O_J1) %*% K1 %*% O_J1 else matrix(0, 0, 0)
  Lambda_2J <- if (k > 0) t(O_J2) %*% K2 %*% O_J2 else matrix(0, 0, 0)
  Lambda_1I <- if (r1 > k) t(O_I1) %*% K1 %*% O_I1 else matrix(0, 0, 0)
  Lambda_2I <- if (r2 > k) t(O_I2) %*% K2 %*% O_I2 else matrix(0, 0, 0)

  list(
    O_J1 = O_J1, O_J2 = O_J2, O_J = O_J,
    O_I1 = O_I1, O_I2 = O_I2,
    Lambda_1J = Lambda_1J, Lambda_2J = Lambda_2J,
    Lambda_1I = Lambda_1I, Lambda_2I = Lambda_2I,
    H1_joint      = if (k > 0) O_J1 %*% Lambda_1J %*% t(O_J1) else matrix(0, n, n),
    H1_individual = if (r1 > k) O_I1 %*% Lambda_1I %*% t(O_I1) else matrix(0, n, n),
    H2_joint      = if (k > 0) O_J2 %*% Lambda_2J %*% t(O_J2) else matrix(0, n, n),
    H2_individual = if (r2 > k) O_I2 %*% Lambda_2I %*% t(O_I2) else matrix(0, n, n)
  )
}
