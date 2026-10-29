#' Select the signal rank of a similarity matrix
#'
#' Chooses how many leading eigenvectors of a symmetric similarity matrix form
#' its signal subspace.
#'
#' @param H A symmetric \eqn{n \times n} similarity matrix.
#' @param method `"variance"` keeps the smallest number of eigenvectors whose
#'   eigenvalues account for at least `var_thresh` of the total.
#'   `"relative_diff"` places the cut at the largest relative drop between
#'   consecutive eigenvalues (elbow).
#' @param var_thresh Proportion of variance for `method = "variance"`.
#' @param tol Eigenvalues at or below `tol` are treated as zero.
#' @param cap_centered_n Optional upper bound on the rank, for example
#'   \eqn{n - 1} for a double-centered matrix.
#' @param min_rank Lower bound on the returned rank.
#' @return An integer rank.
#' @keywords internal
#' @noRd
select_signal_rank <- function(H,
                               method = c("variance", "relative_diff"),
                               var_thresh = 0.80,
                               tol = 1e-16,
                               cap_centered_n = NULL,
                               min_rank = 0L) {
  method <- match.arg(method)

  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  ev <- pmax(ev, 0)          # clip tiny negatives
  ev <- ev[ev > tol]         # drop near-zero eigenvalues
  if (length(ev) == 0L) return(0L)

  max_possible <- length(ev)
  if (!is.null(cap_centered_n)) max_possible <- min(max_possible, cap_centered_n)

  if (method == "variance") {
    cumvar <- cumsum(ev) / sum(ev)
    k <- which(cumvar >= var_thresh)[1]
    if (is.na(k)) k <- max_possible
  } else {
    # Eigenvalues that are zero up to rounding error (for example beyond the
    # number of features) would otherwise produce a spurious largest drop.
    ev <- ev[ev > max(ev) * sqrt(.Machine$double.eps)]
    L <- length(ev)
    if (L <= 1L) {
      k <- L
    } else {
      denom <- pmax(ev[-L], 1e-12)
      drops <- (ev[-L] - ev[-1]) / denom
      k <- which.max(drops)
    }
  }

  as.integer(max(min(k, max_possible), min_rank))
}
