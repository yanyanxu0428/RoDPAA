#' Robust Dissimilarity-based Principal Angle Analysis (RoDPAA)
#'
#' Identifies participant-level variation shared across two paired data views
#' using dissimilarity-based signal-subspace analysis. Each view may use a
#' dissimilarity measure appropriate to its data characteristics, such as
#' Bray-Curtis dissimilarity for microbial abundance profiles and Euclidean
#' distance for metabolomic measurements. Principal-angle analysis identifies
#' corresponding directions across the two views, allowing statistically
#' supported shared variation to be detected even when the directions are only
#' partially aligned.
#'
#' @details
#' ## Input
#' The two views must contain the same \eqn{n}{n} participants in the same
#' order. Each view can be supplied in one of three forms, selected by
#' `input`:
#' * `"dissimilarity"` (default): an \eqn{n \times n}{n x n} dissimilarity
#'   matrix or `dist` object, computed with any dissimilarity appropriate to
#'   the data.
#' * `"data"`: an \eqn{n \times p_v}{n x p_v} data matrix. Euclidean distance
#'   (`"euclidean"`) or Bray-Curtis dissimilarity (`"bray"`) is computed
#'   internally, as selected by `dissimilarity1` and `dissimilarity2`. For any
#'   other dissimilarity, compute the matrix first and use
#'   `input = "dissimilarity"`.
#' * `"kernel"`: an \eqn{n \times n}{n x n} symmetric similarity matrix.
#'
#' The inputs are used as supplied. Any preprocessing, such as a centered
#' log-ratio or logarithmic transformation, must be applied beforehand.
#'
#' ## Method
#' RoDPAA proceeds in five steps:
#' 1. *Similarity matrices.* Each dissimilarity matrix \eqn{D_v}{D_v} is
#'    converted into a centered similarity matrix by Gower double-centering,
#'    \eqn{K_v = -\frac{1}{2} C D_v^{(2)} C}{K_v = -(1/2) C D_v^2 C}, where
#'    \eqn{C = I_n - \frac{1}{n} \mathbf{1}\mathbf{1}^\top}{C = I_n - (1/n) 1 1'}
#'    and \eqn{D_v^{(2)}}{D_v^2} contains the squared dissimilarities (see
#'    [dist_to_kernel()]). Kernel inputs are centered as
#'    \eqn{C K_v C}{C K_v C}.
#' 2. *Signal subspaces.* The leading \eqn{r_v}{r_v} eigenvectors associated
#'    with positive eigenvalues define the signal subspace \eqn{O_v}{O_v} of
#'    each view, so the similarity matrices need not be positive
#'    semidefinite.
#' 3. *Principal-angle alignment.* The singular values
#'    \eqn{\sigma_1 \ge \cdots \ge \sigma_L}{sigma_1 >= ... >= sigma_L} of
#'    \eqn{G = O_1^\top O_2}{G = t(O_1) O_2}, where
#'    \eqn{L = \min(r_1, r_2)}{L = min(r_1, r_2)}, are the shared
#'    intensities: the cosines of the principal angles between the two signal
#'    subspaces.
#' 4. *Permutation inference.* Participant correspondence is disrupted by
#'    permuting the rows of \eqn{O_2}{O_2}, without re-estimating the signal
#'    subspaces. The resulting permutation distributions of the shared
#'    intensities are used to assess their significance and to estimate the
#'    shared dimension \eqn{\hat{k}}{k_hat}.
#' 5. *Shared factors.* For each of the first \eqn{\hat{k}}{k_hat} components,
#'    the paired principal vectors of the two views are averaged to obtain a
#'    shared factor, which assigns one score to each participant.
#'
#' ## Signal ranks
#' The signal ranks \eqn{r_1}{r_1} and \eqn{r_2}{r_2} are selected by
#' `rank_method` unless they are supplied as `r1` and `r2`:
#' * `"variance"` (default) selects the smallest number of eigenvalues that
#'   account for at least `var_thresh` of the sum of the positive eigenvalues.
#' * `"relative_diff"` (elbow) places the cut at the largest relative drop
#'   between consecutive positive eigenvalues. It is suited to eigenvalues of
#'   which a few are clearly separated from the rest.
#'
#' ## Permutation inference
#' RoDPAA uses permutation testing to assess the significance of
#' principal-angle shared intensities and estimate the shared dimension
#' (\eqn{\hat{k}}{k_hat}). Empirical P-values are adjusted for multiple
#' testing using `correction_method`. The estimated shared dimension is the
#' number of consecutive significant leading components.
#'
#' Users may alternatively specify the shared dimension using `k` or disable
#' permutation inference with `use_permutation = FALSE`. The `k_extract`
#' argument allows additional directions to be extracted for exploratory
#' analysis.
#'
#' Use [set.seed()] before calling `rodpaa()` for reproducible permutation
#' results.
#'
#' ## Interpreting the shared factors
#' Each shared factor represents a participant-level pattern of variation
#' identified through principal-angle alignment between the two data views.
#' The paired principal vectors of the two views are consistently oriented
#' relative to each other by the singular value decomposition, although the
#' overall sign of each shared factor is arbitrary. Shared factors can be
#' interpreted by examining associations between their scores and the
#' original features of each view, or with participant characteristics, such
#' as clinical outcomes. Because RoDPAA estimates the shared factors without
#' using participant characteristics, such associations can be investigated
#' in downstream analyses.
#'
#' @param x1,x2 The two paired data views, containing the same participants
#'   in the same row order. Depending on `input`, provide numeric data
#'   matrices (`"data"`), dissimilarity matrices or `dist` objects
#'   (`"dissimilarity"`), or symmetric similarity matrices (`"kernel"`).
#' @param input Type of matrices supplied as `x1` and `x2`:
#'   * `"dissimilarity"` (default): dissimilarity matrices or `dist` objects;
#'   * `"data"`: numeric data matrices, from which the dissimilarities are
#'     computed as specified by `dissimilarity1` and `dissimilarity2`;
#'   * `"kernel"`: symmetric similarity (kernel) matrices.
#'
#'   See "Input" in Details.
#' @param dissimilarity1,dissimilarity2 Dissimilarity computed for each view
#'   when `input = "data"`: `"euclidean"` (default) or `"bray"`. Bray-Curtis
#'   dissimilarity requires nonnegative data, such as counts or relative
#'   abundances. For any other dissimilarity, compute the dissimilarity
#'   matrices first and supply them with `input = "dissimilarity"`.
#' @param k Optional number of shared directions, specified in advance. If
#'   `NULL` (default), the shared dimension is estimated from the data by the
#'   permutation test. If a number is supplied, the first `k` components are
#'   treated as shared and no permutation test is performed.
#' @param r1,r2 Optional signal ranks of the two views, specified in advance.
#'   If `NULL` (default), they are estimated from the data using
#'   `rank_method`.
#' @param rank_method Rule for estimating the signal ranks when `r1` and `r2`
#'   are not supplied:
#'   * `"variance"` (default): the smallest number of eigenvalues that account
#'     for at least `var_thresh` (by default 80%) of the sum of the positive
#'     eigenvalues;
#'   * `"relative_diff"`: the cut at the largest relative drop between
#'     consecutive positive eigenvalues (elbow).
#' @param var_thresh Proportion of the sum of the positive eigenvalues used by
#'   `rank_method = "variance"`.
#' @param use_permutation Select the shared dimension by the permutation test?
#'   If `FALSE`, the largest relative drop in the shared intensities is used.
#' @param n_perm Number of permutations \eqn{B}{B}.
#' @param alpha Significance level, applied to the adjusted P-values.
#' @param use_multiple_testing Logical. If `TRUE` (default), the permutation
#'   P-values are adjusted for multiple testing across the \eqn{L}{L}
#'   components, using `correction_method`. If `FALSE`, the unadjusted
#'   P-values are compared with `alpha`.
#' @param correction_method Multiple-testing correction: `"bonferroni"`
#'   (default), `"holm"` or `"fdr"` (Benjamini-Hochberg).
#' @param k_extract Optional number of shared factors to extract. If `NULL`
#'   (default), the number of shared factors equals the estimated shared
#'   dimension (returned as `k_hat`). A larger number can be supplied to
#'   inspect further directions. Directions beyond the estimated shared
#'   dimension are exploratory: they are not supported by the permutation
#'   test, and `k_hat` is reported unchanged.
#' @param verbose Logical. If `TRUE` (default), a message reports the shared
#'   dimension selected by the permutation test, for example
#'   "Shared dimension selected by the permutation test: k = 2", and notes
#'   when exploratory directions are extracted with `k_extract`. If `FALSE`,
#'   no messages are printed. The results are the same either way.
#'
#' @return An object of class `"rodpaa"`, a list containing the estimated
#'   shared factors, shared dimension, shared intensities, and permutation
#'   results. Its components are listed below. See
#'   `vignette("Getting-started-with-RoDPAA", package = "RoDPAA")` for
#'   examples.
#'
#'   \strong{Main results}
#' \describe{
#'   \item{`shared_factors`}{\eqn{n \times}{n x} `k_extracted` matrix of
#'     shared factors, with columns `SF1`, `SF2`, and so on.}
#'   \item{`k_hat`}{Estimated shared dimension.}
#'   \item{`shared_intensities`}{Shared intensities (cosines of the principal
#'     angles), in decreasing order.}
#'   \item{`dimensions`}{List with `k_hat`, `k_extracted`, and the signal
#'     ranks `r1` and `r2`. `k_extracted` is the number of extracted
#'     directions, which determines the number of columns of the shared
#'     factors and directions. It equals `k_hat` unless `k_extract` is
#'     supplied.}
#'   \item{`permutation`}{Results of the permutation test, or `NULL` if it was
#'     not performed. The element `table` is a data frame with one row per
#'     component and columns `component`, `singular_value` (shared
#'     intensity), `p_value`, `p_adj` and `significant`. The element
#'     `sigma_perm` is the \eqn{B \times L}{B x L} matrix of permuted shared
#'     intensities.}
#'   \item{`p_values`, `p_adj`, `significant`}{Raw and adjusted P-values and
#'     significance of each component. Present only if the permutation test
#'     was performed.}
#' }
#'
#'   \strong{Directions and eigenvalues}
#' \describe{
#'   \item{`O_J1`, `O_J2`, `O_J`}{Paired shared directions of view 1 and
#'     view 2, and their average. `O_J` is identical to `shared_factors`.}
#'   \item{`O_I1`, `O_I2`}{View-specific directions: the remaining
#'     \eqn{r_v}{r_v} \eqn{-}{-} `k_extracted` directions of each signal
#'     subspace.}
#'   \item{`eigenvalues1`, `eigenvalues2`}{Eigenvalues of \eqn{K_1}{K_1} and
#'     \eqn{K_2}{K_2}, with negative values set to zero.}
#'   \item{`singular_values_G`}{Identical to `shared_intensities`.}
#' }
#'
#' @seealso [plot.rodpaa()] to plot the results, and
#'   `vignette("Getting-started-with-RoDPAA", package = "RoDPAA")` for a
#'   tutorial.
#'
#' @export
#' @examples
#' ## Simulated data: two views measured on the same 80 participants
#' X1 <- rodpaa_example$view1      # 80 participants x 50 features
#' X2 <- rodpaa_example$view2      # 80 participants x 40 features
#'
#' ## Dissimilarity matrices (Euclidean distance)
#' D1 <- dist(X1)                  # 80 x 80
#' D2 <- dist(X2)                  # 80 x 80
#'
#' ## Similarity (kernel) matrices, by Gower double-centering
#' K1 <- dist_to_kernel(D1)        # 80 x 80
#' K2 <- dist_to_kernel(D2)        # 80 x 80
#'
#' ## 1. Dissimilarity matrices (default input)
#' set.seed(1)
#' fit <- rodpaa(D1, D2, input = "dissimilarity")
#'
#' ## 2. Data matrices, with the dissimilarities computed internally
#' set.seed(1)
#' fit_data <- rodpaa(X1, X2, input = "data",
#'                    dissimilarity1 = "euclidean", dissimilarity2 = "euclidean")
#'
#' ## 3. Similarity (kernel) matrices
#' set.seed(1)
#' fit_kernel <- rodpaa(K1, K2, input = "kernel")
#'
#' ## The three inputs give the same shared dimension
#' c(fit$k_hat, fit_data$k_hat, fit_kernel$k_hat)    # 2 2 2
#'
#' ## Results
#' fit
#' c(r1 = fit$dimensions$r1, r2 = fit$dimensions$r2) # signal ranks: 3 and 3
#' fit$k_hat                       # shared dimension: 2
#' round(fit$shared_intensities, 3)                  # 0.978 0.967 0.081
#' fit$permutation$table           # adjusted P-values: 0, 0 and 0.489
#' head(fit$shared_factors)        # shared factors: 80 x 2 matrix
#' plot(fit)
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
      if (verbose) message("Shared dimension selected by the permutation test: k = ", k_hat_perm)
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
    message("Extracting ", k_use, " directions (shared dimension k = ", k_hat, "). ",
            "Directions beyond k are exploratory.")
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
