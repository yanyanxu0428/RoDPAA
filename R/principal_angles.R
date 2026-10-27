#' Principal-angle alignment of two signal subspaces (Step 3)
#'
#' Computes \eqn{G = O_1^\top O_2} and its singular value decomposition
#' \eqn{G = P \Sigma Q^\top}. The singular values (shared intensities) are the
#' cosines of the principal angles between \eqn{\mathrm{span}(O_1)} and
#' \eqn{\mathrm{span}(O_2)}.
#'
#' @param O1,O2 Orthonormal \eqn{n \times r_1} and \eqn{n \times r_2} bases of
#'   the two signal subspaces.
#' @return A list with `G`, the singular values `sigma` (\eqn{\min(r_1, r_2)}
#'   of them, decreasing), and the full singular vector matrices `P`
#'   (\eqn{r_1 \times r_1}) and `Q` (\eqn{r_2 \times r_2}). The first columns
#'   of `P` and `Q` give the paired shared directions; the remaining columns
#'   give the view-specific directions.
#' @keywords internal
#' @noRd
principal_angles <- function(O1, O2) {
  G <- t(O1) %*% O2
  s <- svd(G, nu = ncol(O1), nv = ncol(O2))
  list(G = G, sigma = s$d, P = s$u, Q = s$v)
}

#' Elbow estimate of the shared dimension
#'
#' Largest relative drop between consecutive shared intensities. Used only
#' when `use_permutation = FALSE`.
#'
#' @param sigma Shared intensities (decreasing).
#' @param max_k Upper bound, \eqn{\min(r_1, r_2)}.
#' @return An integer.
#' @keywords internal
#' @noRd
elbow_shared_dimension <- function(sigma, max_k) {
  s <- sigma[sigma > 1e-16]
  L <- length(s)
  if (L <= 1) return(L)
  drops <- (s[-L] - s[-1]) / s[-L]
  min(which.max(drops), max_k)
}
