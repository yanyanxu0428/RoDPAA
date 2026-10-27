#' Convert a dissimilarity matrix to a centered similarity kernel
#'
#' Constructs a centered sample-by-sample similarity matrix from a
#' user-specified dissimilarity matrix using Gower double-centering:
#'
#' \deqn{
#' K = -\frac{1}{2} C D^{(2)} C,
#' \qquad
#' C = I_n - \frac{1}{n}\mathbf{1}\mathbf{1}^{\top}.
#' }
#'
#' Here, \eqn{D^{(2)}} denotes the element-wise squared dissimilarity
#' matrix, and \eqn{n} is the number of samples.
#'
#' @details
#' For Euclidean distances, Gower double-centering produces the
#' centered linear kernel, which is positive semidefinite.
#'
#' For non-Euclidean dissimilarities, Gower double-centering may
#' produce a matrix with negative eigenvalues. RoDPAA sets these
#' eigenvalues to zero, retaining the positive spectral component
#' used to estimate the signal subspace. Consequently, the
#' resulting kernel may differ from the original double-centered
#' matrix.
#'
#' @param D A symmetric, nonnegative \eqn{n \times n}
#'   dissimilarity matrix with zero diagonal, or a `dist` object.
#'   Rows and columns must represent the same samples in the
#'   same order.
#'
#' @return A centered, positive semidefinite \eqn{n \times n}
#'   similarity matrix. Sample names are preserved when available.
#'
#' @export
#'
#' @examples
#' # Euclidean dissimilarity
#' set.seed(123)
#' X <- matrix(rnorm(50), nrow = 10)
#' D <- dist(X)
#' K <- dist_to_kernel(D)
#'
#' dim(K)
#' # eigenvalues are nonnegative (up to rounding error of about 1e-15)
#' round(eigen(K, symmetric = TRUE)$values, 10)
dist_to_kernel <- function(D) {
  D <- as.matrix(D)
  if (nrow(D) != ncol(D)) stop("`D` must be a square dissimilarity matrix.")
  n <- nrow(D)
  C <- diag(n) - 1 / n           # centring matrix C = I_n - (1/n) 1 1^T
  D2 <- D * D                    # D^(2): element-wise squared dissimilarities
  K <- -0.5 * C %*% D2 %*% C     # K = -1/2 C D^(2) C
  # set negative eigenvalues of K to zero (positive semi-definite)
  eK <- eigen(K, symmetric = TRUE)
  K <- eK$vectors %*% diag(pmax(0, eK$values), nrow = n) %*% t(eK$vectors)
  dimnames(K) <- dimnames(D)
  K
}

#' Double-center a user-supplied similarity matrix
#'
#' Checks that `K` is a square, symmetric matrix and returns \eqn{C K C}.
#' A matrix that is already centered is unchanged, up to rounding error.
#'
#' @param K An \eqn{n \times n} similarity matrix.
#' @param name Name of the argument, for error messages.
#' @return \eqn{C K C} with the centering matrix
#'   \eqn{C = I_n - \frac{1}{n} 1 1^\top}, with the row and column names of `K`.
#' @keywords internal
#' @noRd
center_kernel <- function(K, name = "K") {
  K <- as.matrix(K)
  if (nrow(K) != ncol(K)) stop("`", name, "` must be a square similarity matrix.")
  if (!isSymmetric(unname(K), tol = sqrt(.Machine$double.eps))) {
    stop("`", name, "` is not symmetric. With `input = \"kernel\"`, each view must be a ",
         "symmetric similarity matrix. For a dissimilarity matrix, use `input = \"dissimilarity\"`.")
  }
  if (looks_like_dissimilarity(K)) {
    stop("`", name, "` looks like a dissimilarity matrix (zero diagonal), not a kernel. ",
         "Use `input = \"dissimilarity\"`.")
  }
  n <- nrow(K)
  C <- diag(n) - matrix(1 / n, n, n)
  KC <- C %*% K %*% C
  dimnames(KC) <- dimnames(K)
  KC
}

#' Linear kernel of a data matrix
#'
#' @param X An \eqn{n \times p} data matrix (participants in rows).
#' @return The double-centred linear kernel \eqn{C X X^\top C}, with the
#'   centring matrix \eqn{C = I_n - \frac{1}{n} 1 1^\top}. It equals
#'   `dist_to_kernel(dist(X))`, so with this kernel RoDPAA uses Euclidean
#'   geometry.
#' @keywords internal
#' @noRd
linear_kernel <- function(X) {
  center_kernel(X %*% t(X))
}
