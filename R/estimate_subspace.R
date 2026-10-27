#' Estimate the signal subspace of one view (Step 2)
#'
#' Eigendecomposes the centered similarity matrix \eqn{K_v}, chooses the signal
#' rank \eqn{\hat{r}_v} (unless supplied) and keeps the leading \eqn{\hat{r}_v}
#' eigenvectors as the signal subspace \eqn{\hat{O}_v}. The remaining
#' eigenvectors form the noise part.
#'
#' @param K Centered \eqn{n \times n} similarity matrix of one view.
#' @param r Optional signal rank. If `NULL`, chosen by `select_signal_rank()`.
#' @param rank_method,var_thresh Passed to `select_signal_rank()`.
#' @return A list with the signal rank `r`, the signal basis `O` (\eqn{n \times
#'   r}), its eigenvalue matrix `Lambda`, the signal part `K_signal`
#'   \eqn{= O \Lambda O^\top}, the noise part `E` and all eigenvalues `values`
#'   (negative values set to zero).
#' @keywords internal
#' @noRd
estimate_signal_subspace <- function(K, r = NULL, rank_method = "variance", var_thresh = 0.80) {
  n <- nrow(K)
  eig <- eigen(K, symmetric = TRUE)
  eig$values <- pmax(eig$values, 0)

  if (is.null(r)) {
    r <- select_signal_rank(K, method = rank_method,
                            var_thresh = var_thresh, tol = 1e-16,
                            cap_centered_n = n - 1L, min_rank = 0L)
  }

  O      <- if (r > 0) eig$vectors[, 1:r, drop = FALSE] else matrix(0, n, 0)
  Lambda <- if (r > 0) diag(eig$values[1:r], nrow = r) else matrix(0, 0, 0)
  K_signal <- if (r > 0) O %*% Lambda %*% t(O) else matrix(0, n, n)

  # noise part: the remaining eigenvectors
  if (r < n) {
    O_noise <- eig$vectors[, (r + 1):n, drop = FALSE]
    E <- O_noise %*% diag(eig$values[(r + 1):n], nrow = n - r) %*% t(O_noise)
  } else {
    E <- matrix(0, n, n)
  }

  list(r = r, O = O, Lambda = Lambda, K_signal = K_signal, E = E, values = eig$values)
}
