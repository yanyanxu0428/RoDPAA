#' Permutation test for the shared dimension
#'
#' Compares the singular values of \eqn{G = O_1^\top O_2} (the cosines of the
#' principal angles between the two signal subspaces) with their distribution
#' when the participant correspondence between the views is broken by permuting
#' the rows of \eqn{O_2}. The estimated shared dimension is the number of
#' leading singular values that are significant, counted up to the first
#' non-significant one.
#'
#' Set a seed with [set.seed()] before calling this function for reproducible
#' results.
#'
#' @param O1_signal,O2_signal Orthonormal \eqn{n \times r_1} and
#'   \eqn{n \times r_2} bases of the two signal subspaces.
#' @param n_perm Number of permutations.
#' @param alpha Significance level.
#' @param use_multiple_testing Adjust the p-values across components?
#' @param correction_method `"bonferroni"`, `"holm"` or `"fdr"`.
#' @return A list with the estimated shared dimension `k_hat`, the observed
#'   singular values `sigma_orig`, the permutation distribution `sigma_perm`
#'   (one row per permutation), its summaries, raw and adjusted p-values, and
#'   the settings used.
#' @keywords internal
#' @noRd
test_shared_dimension <- function(O1_signal, O2_signal,
                                  n_perm = 1000,
                                  alpha = 0.05,
                                  use_multiple_testing = TRUE,
                                  correction_method = "bonferroni") {
  G_orig <- t(O1_signal) %*% O2_signal
  sigma_orig <- svd(G_orig, nu = 0, nv = 0)$d
  max_rank <- length(sigma_orig)

  n <- nrow(O1_signal)
  sigma_perm <- matrix(NA_real_, nrow = n_perm, ncol = max_rank)
  for (i in seq_len(n_perm)) {
    O2_perm <- O2_signal[sample(n), , drop = FALSE]
    G_perm  <- t(O1_signal) %*% O2_perm
    svd_p   <- svd(G_perm, nu = 0, nv = 0)$d
    len <- min(length(svd_p), max_rank)
    sigma_perm[i, seq_len(len)] <- svd_p[seq_len(len)]
  }

  p_values <- colMeans(sigma_perm >= matrix(sigma_orig, nrow(sigma_perm), ncol(sigma_perm), byrow = TRUE),
                       na.rm = TRUE)

  if (use_multiple_testing) {
    p_adj <- switch(correction_method,
                    bonferroni = pmin(p_values * max_rank, 1),
                    holm       = stats::p.adjust(p_values, method = "holm"),
                    fdr        = stats::p.adjust(p_values, method = "fdr"),
                    { warning("Unknown correction method; using raw p-values."); p_values }
    )
  } else {
    p_adj <- p_values
  }

  significant <- p_adj < alpha
  first_false <- match(FALSE, significant, nomatch = max_rank + 1)
  k_hat <- first_false - 1

  list(
    k_hat = k_hat,
    sigma_orig = sigma_orig,
    sigma_perm = sigma_perm,
    perm_mean = colMeans(sigma_perm, na.rm = TRUE),
    perm_sd   = apply(sigma_perm, 2, stats::sd, na.rm = TRUE),
    perm_q95  = apply(sigma_perm, 2, function(x) stats::quantile(x, 0.95, na.rm = TRUE)),
    perm_q99  = apply(sigma_perm, 2, function(x) stats::quantile(x, 0.99, na.rm = TRUE)),
    p_values = p_values,
    p_adj = p_adj,
    significant = significant,
    alpha = alpha,
    use_multiple_testing = use_multiple_testing,
    correction_method = if (use_multiple_testing) correction_method else "none"
  )
}
