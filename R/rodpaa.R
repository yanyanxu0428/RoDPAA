#' Robust Dissimilarity-based Principal Angle Analysis (RoDPAA)
#'
#' Identifies participant-level variation shared by two data views measured on
#' the same \eqn{n} participants. Each view is represented by its own
#' dissimilarity, so each can keep the geometry that suits it, for example
#' Bray-Curtis for microbiome counts and Euclidean distance for log-transformed
#' metabolomic profiles.
#'
#' Each view can be given in one of three forms (`input`):
#' * `"dissimilarity"` (default): dissimilarity matrices, for any
#'   dissimilarity, for example Bray-Curtis, UniFrac or Aitchison distance.
#' * `"data"`: the data matrices, as a shortcut. The dissimilarities are
#'   computed with `dissimilarity1` and `dissimilarity2`. The data are used as
#'   given: apply any transformation (for example CLR or log) or filtering
#'   first.
#' * `"kernel"` (advanced): similarity (kernel) matrices.
#'
#' The steps are:
#' 1. convert each dissimilarity matrix to a centered similarity matrix
#'    ([dist_to_kernel()]);
#' 2. keep the leading eigenvectors of each similarity matrix as its signal
#'    subspace, with signal ranks \eqn{r_1} and \eqn{r_2} chosen by the
#'    proportion of variance explained (`var_thresh`) unless supplied;
#' 3. take the singular value decomposition of \eqn{G = O_1^\top O_2}. The
#'    singular values (shared intensities) are the cosines of the principal
#'    angles between the two subspaces;
#' 4. estimate the shared dimension \eqn{\hat{k}} with a permutation test that
#'    breaks participant correspondence. \eqn{\hat{k}} is the number of
#'    leading singular values that are significant, counted up to the first
#'    non-significant one;
#' 5. average each significant pair of directions across the views to give the
#'    shared factors.
#'
#' The permutation test is random. Call [set.seed()] first for reproducible
#' results.
#'
#' @param x1,x2 The two views on the same participants, in the same order
#'   (they are never reordered). Depending on `input`:
#'   * `"dissimilarity"`: \eqn{n \times n} dissimilarity matrices or `dist`
#'     objects (symmetric, nonnegative, zero diagonal). A distance (for example
#'     Euclidean) is a particular type of dissimilarity.
#'   * `"data"`: numeric \eqn{n \times p_1} and \eqn{n \times p_2} matrices or
#'     data frames, participants in rows and features in columns, with no
#'     missing values.
#'   * `"kernel"`: \eqn{n \times n} symmetric similarity matrices.
#'
#'   If both views have row names, they must be identical.
#' @param input `"dissimilarity"` (default), `"data"` or `"kernel"`. Kernel
#'   matrices are double-centered as \eqn{C K C} before analysis, so they need
#'   not be centered beforehand. Kernels that are already centered are
#'   unchanged. Positive semidefiniteness is not required: only eigenvectors
#'   with positive eigenvalues are used.
#' @param dissimilarity1,dissimilarity2 Dissimilarity for each view when
#'   `input = "data"`: `"euclidean"` (default) or `"bray"` (Bray-Curtis, for
#'   nonnegative data such as counts or relative abundances). Ignored for
#'   the other inputs.
#' @param k Optional shared dimension. If given, the permutation test is skipped
#'   and `k` directions are treated as shared.
#' @param r1,r2 Optional signal ranks of the two views. If `NULL`, chosen by
#'   `rank_method` and `var_thresh`.
#' @param rank_method Signal-rank rule. `"variance"` (default) keeps the
#'   smallest number of eigenvectors whose eigenvalues account for at least
#'   `var_thresh` of the total. `"relative_diff"` cuts at the largest relative
#'   drop between consecutive eigenvalues.
#' @param var_thresh Variance threshold for `rank_method = "variance"`.
#' @param use_permutation Use the permutation test to choose \eqn{\hat{k}}? If
#'   `FALSE`, the largest relative drop in the singular values (elbow) is used.
#' @param n_perm Number of permutations.
#' @param alpha Significance level of the permutation test.
#' @param use_multiple_testing Adjust the permutation p-values across
#'   components?
#' @param correction_method `"bonferroni"` (default), `"holm"` or `"fdr"`.
#' @param k_extract Optional number of directions to extract, which may differ
#'   from \eqn{\hat{k}}, for example to visualise the leading directions when
#'   none is significant. The test result is still reported in `k_hat`.
#' @param verbose Print the selected shared dimension?
#'
#' @return An object of class `"rodpaa"`, a list including
#' \describe{
#'   \item{shared_factors}{\eqn{n \times k} matrix of shared factors (the
#'     average of the paired directions of the two views). As with principal
#'     components, the sign of each factor is arbitrary: a factor and its
#'     negative describe the same shared variation.}
#'   \item{k_hat}{Estimated shared dimension.}
#'   \item{shared_intensities}{Singular values of \eqn{G}, the cosines of the
#'     principal angles.}
#'   \item{O_J1, O_J2, O_J}{Shared directions of view 1, view 2 and their
#'     average (`O_J` equals `shared_factors`).}
#'   \item{O_I1, O_I2}{View-specific (individual) directions.}
#'   \item{dimensions}{List with `k_hat`, `k_extracted`, `r1` and `r2`.}
#'   \item{permutation}{Permutation results, including a `table` with the
#'     singular value, raw and adjusted p-value of each component, or `NULL`
#'     when the test was not run.}
#' }
#' and the joint and individual reconstructions of the two similarity
#' matrices.
#'
#' @export
#' @examples
#' # two views on the same 80 participants (simulated, two shared factors)
#' D1 <- dist(rodpaa_example$view1)   # Euclidean dissimilarity, view 1
#' D2 <- dist(rodpaa_example$view2)   # Euclidean dissimilarity, view 2
#'
#' set.seed(1)                        # the permutation test is random
#' fit <- rodpaa(D1, D2, input = "dissimilarity")
#' fit
#' summary(fit)
#' plot(fit)
#' head(fit$shared_factors)
#'
#' # the same fit directly from the data matrices
#' set.seed(1)
#' fit_x <- rodpaa(rodpaa_example$view1, rodpaa_example$view2, input = "data",
#'                 dissimilarity1 = "euclidean", dissimilarity2 = "euclidean")
#' all.equal(fit$shared_factors, fit_x$shared_factors)
rodpaa <- function(x1, x2,
                   input = c("dissimilarity", "data", "kernel"),
                   dissimilarity1 = c("euclidean", "bray"),
                   dissimilarity2 = c("euclidean", "bray"),
                   k = NULL, r1 = NULL, r2 = NULL,
                   rank_method = c("variance", "relative_diff"),
                   var_thresh = 0.80,
                   use_permutation = TRUE,
                   n_perm = 1000, alpha = 0.05,
                   use_multiple_testing = TRUE,
                   correction_method = "bonferroni",
                   k_extract = NULL,
                   verbose = TRUE) {
  input <- match.arg(input)
  rank_method <- match.arg(rank_method)
  dissimilarity1 <- match.arg(dissimilarity1)
  dissimilarity2 <- match.arg(dissimilarity2)

  # Step 1: data -> dissimilarities -> centered similarity matrices (Gower double-centering)
  if (input == "data") {
    x1 <- as_data_matrix(x1, "x1")
    x2 <- as_data_matrix(x2, "x2")
    if (nrow(x1) != nrow(x2)) {
      stop("`x1` and `x2` have different numbers of participants (rows): ",
           nrow(x1), " and ", nrow(x2), ".")
    }
    x1 <- compute_dissimilarity(x1, dissimilarity1, "x1")
    x2 <- compute_dissimilarity(x2, dissimilarity2, "x2")
  }
  if (input %in% c("data", "dissimilarity")) {
    H1 <- dist_to_kernel(check_dissimilarity(x1, "x1"))
    H2 <- dist_to_kernel(check_dissimilarity(x2, "x2"))
  } else {
    H1 <- center_kernel(x1, "x1")
    H2 <- center_kernel(x2, "x2")
  }
  n <- nrow(H1)
  if (ncol(H1) != n || nrow(H2) != n || ncol(H2) != n) {
    stop("Both views must be n x n matrices on the same n participants.")
  }
  ids <- rownames(H1)
  if (!is.null(ids) && !is.null(rownames(H2)) && !identical(ids, rownames(H2))) {
    stop("The two views have different row names. Put the participants in the same order.")
  }

  # Step 2: signal ranks and signal subspaces
  S1 <- estimate_signal_subspace(H1, r = r1, rank_method = rank_method, var_thresh = var_thresh)
  S2 <- estimate_signal_subspace(H2, r = r2, rank_method = rank_method, var_thresh = var_thresh)
  r1 <- S1$r
  r2 <- S2$r

  # Step 3: principal-angle alignment, G = O1^T O2 = P Sigma Q^T
  angles <- principal_angles(S1$O, S2$O)

  # Step 4: shared dimension by permutation inference (elbow if use_permutation = FALSE)
  perm_result <- NULL
  k_hat_perm <- NA_integer_
  k_hat_relative_rank <- 0
  if (is.null(k)) {
    k_hat_relative_rank <- elbow_shared_dimension(angles$sigma, min(r1, r2))
    n_positive <- sum(angles$sigma > 1e-16)
    if (use_permutation && n_positive > 0) {
      perm_result <- test_shared_dimension(
        S1$O, S2$O,
        n_perm = n_perm,
        alpha = alpha,
        use_multiple_testing = use_multiple_testing,
        correction_method = correction_method
      )
      k_hat_perm <- perm_result$k_hat
      k_hat <- k_hat_perm
      if (verbose) message("Permutation test selected k = ", k_hat_perm)
    } else {
      k_hat <- if (n_positive == 0) 0 else k_hat_relative_rank
    }
  } else {
    k_hat <- k
  }

  # Number of directions actually extracted
  k_use <- if (!is.null(k_extract)) k_extract else k_hat
  if (k_use > min(r1, r2)) {
    stop("Cannot extract ", k_use, " shared directions: the signal ranks are r1 = ", r1,
         " and r2 = ", r2, ". Use a smaller `k`/`k_extract` or larger `r1`, `r2`.")
  }
  if (!is.null(k_extract) && verbose && k_use != k_hat) {
    message("Extracting k = ", k_use, " directions (shared dimension k_hat = ", k_hat, ").")
  }

  # Step 5: shared factors = average of the paired shared directions
  SF <- construct_shared_factors(S1, S2, angles, H1, H2, k_use)
  O_J <- SF$O_J
  if (k_use > 0) {
    colnames(O_J) <- paste0("SF", seq_len(k_use))
    if (!is.null(ids)) rownames(O_J) <- ids
  }

  # Return results
  result <- list(
    shared_factors = O_J,
    k_hat = k_hat,
    shared_intensities = angles$sigma,
    O_J1 = SF$O_J1, O_J2 = SF$O_J2, O_J = O_J,
    Lambda_1J = SF$Lambda_1J, Lambda_2J = SF$Lambda_2J,
    O_I1 = SF$O_I1, O_I2 = SF$O_I2,
    Lambda_1I = SF$Lambda_1I, Lambda_2I = SF$Lambda_2I,
    E1 = S1$E, E2 = S2$E,
    H1_joint = SF$H1_joint, H1_individual = SF$H1_individual,
    H2_joint = SF$H2_joint, H2_individual = SF$H2_individual,
    H1_signal = S1$K_signal, H2_signal = S2$K_signal,
    dimensions = list(k_hat = k_hat, k_extracted = k_use, r1 = r1, r2 = r2),
    singular_values_G = angles$sigma,
    eigenvalues1 = S1$values, eigenvalues2 = S2$values,
    use_permutation = use_permutation,
    rank_selection = list(
      k_perm  = k_hat_perm,
      k_elbow = k_hat_relative_rank,
      k_hat   = k_hat,
      k_extracted = k_use,
      rule    = if (!is.null(k)) "user" else if (use_permutation && length(angles$sigma) > 0) "perm" else "elbow"
    ),
    permutation = NULL
  )

  if (!is.null(perm_result)) {
    result$permutation <- list(
      k_hat             = perm_result$k_hat,
      alpha             = perm_result$alpha,
      correction_method = perm_result$correction_method,
      table = data.frame(
        component      = seq_along(perm_result$sigma_orig),
        singular_value = perm_result$sigma_orig,
        p_value        = perm_result$p_values,
        p_adj          = perm_result$p_adj,
        significant    = perm_result$significant
      ),
      sigma_orig = perm_result$sigma_orig,
      sigma_perm = perm_result$sigma_perm,
      perm_mean  = perm_result$perm_mean,
      perm_sd    = perm_result$perm_sd,
      perm_q95   = perm_result$perm_q95,
      perm_q99   = perm_result$perm_q99
    )
    result$p_values    <- perm_result$p_values
    result$p_adj       <- perm_result$p_adj
    result$significant <- perm_result$significant
  }

  class(result) <- "rodpaa"
  result
}
