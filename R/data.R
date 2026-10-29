#' Simulated two-view example data
#'
#' A small simulated dataset for trying out [rodpaa()]. Two views are measured
#' on the same 80 participants. Both views contain two shared latent factors,
#' each view also contains one additional factor, and both contain
#' independent noise. RoDPAA should therefore find two shared factors.
#'
#' @details
#' For \eqn{n = 80} participants and views \eqn{v = 1, 2} with
#' \eqn{p_1 = 50} and \eqn{p_2 = 40} features, the data were generated as
#' \deqn{X_v = F_J W_v + F_{I,v} U_v + E_v,}
#' where
#' * \eqn{F_J} (\eqn{n \times 2}) contains the two shared factors, common to
#'   both views (returned as `shared`);
#' * \eqn{W_v} (\eqn{2 \times p_v}) contains the loadings of the shared
#'   factors in view \eqn{v};
#' * \eqn{F_{I,v}} (\eqn{n \times 1}) is an additional factor of view
#'   \eqn{v}, with loadings \eqn{U_v} (\eqn{1 \times p_v});
#' * \eqn{E_v} (\eqn{n \times p_v}) is noise.
#'
#' The entries of \eqn{F_J}, \eqn{F_{I,v}} and \eqn{W_v} are independent
#' \eqn{N(0, 1)}, the entries of \eqn{U_v} are independent
#' \eqn{N(0, \sigma_{\mathrm{indiv}}^2)} with
#' \eqn{\sigma_{\mathrm{indiv}} = 1.5}, and the entries of \eqn{E_v} are
#' independent \eqn{N(0, \sigma_{\mathrm{noise}}^2)} with
#' \eqn{\sigma_{\mathrm{noise}} = 1}. \eqn{F_{I,1}} and \eqn{F_{I,2}} are
#' independent, so the additional factors are not shared. Each view therefore
#' has signal rank 3, and the two views share two dimensions.
#' The data were generated with `set.seed(2026)`. The full code is in
#' `data-raw/rodpaa_example.R` in the package source on GitHub.
#'
#' @format A list with three matrices, all with participants `"P01"` to
#'   `"P80"` in rows, in the same order:
#' \describe{
#'   \item{view1}{80 x 50 data matrix \eqn{X_1} of view 1.}
#'   \item{view2}{80 x 40 data matrix \eqn{X_2} of view 2.}
#'   \item{shared}{80 x 2 matrix of the true shared factors \eqn{F_J}, for
#'     checking the result.}
#' }
#' @source Simulated, see Details.
#' @examples
#' ## Data: two views measured on the same 80 participants
#' dim(rodpaa_example$view1)               # 80 x 50
#' dim(rodpaa_example$view2)               # 80 x 40
#'
#' ## Dissimilarity matrix of each view (Euclidean distance)
#' D1 <- dist(rodpaa_example$view1)
#' D2 <- dist(rodpaa_example$view2)
#'
#' ## Fit RoDPAA, with the signal ranks selected by the variance rule
#' set.seed(1)                               # the permutation test is random
#' fit <- rodpaa(D1, D2,
#'               input = "dissimilarity",
#'               rank_method = "variance",
#'               var_thresh = 0.80,
#'               n_perm = 1000,
#'               alpha = 0.05,
#'               correction_method = "bonferroni")
#'
#' ## Overview of the fit and its components
#' fit
#' names(fit)
#'
#' ## Signal ranks of the two views
#' c(r1 = fit$dimensions$r1, r2 = fit$dimensions$r2)        # 3 and 3
#'
#' ## Shared dimension, selected by the permutation test
#' fit$k_hat                                                # 2
#'
#' ## Shared intensities (cosines of the principal angles)
#' round(fit$shared_intensities, 3)                         # 0.978 0.967 0.081
#'
#' ## Permutation test of each component
#' summary(fit)
#'
#' ## Shared factors: one score per participant and shared direction
#' dim(fit$shared_factors)                                  # 80 x 2
#' head(round(fit$shared_factors, 3))
#'
#' ## Agreement with the true shared subspace: largest principal angle (radians)
#' ## between the estimated shared factors and the true shared factors F_J,
#' ## where 0 means identical subspaces
#' if (requireNamespace("pracma", quietly = TRUE)) {
#'   pracma::subspace(fit$shared_factors, rodpaa_example$shared)  # 0.12
#' }
#'
#' ## Figures
#' plot(fit)                                 # shared intensities and permutation null
#' plot(fit, type = "scores")                # participant scores on SF1 and SF2
#'
#' ## Signal ranks selected by the elbow rule (largest relative drop)
#' set.seed(1)
#' fit_elbow <- rodpaa(D1, D2, input = "dissimilarity",
#'                     rank_method = "relative_diff")
#' c(r1 = fit_elbow$dimensions$r1, r2 = fit_elbow$dimensions$r2,
#'   k_hat = fit_elbow$k_hat)                               # 3, 3 and 2
#'
#' ## The same analysis, starting from the data matrices
#' set.seed(1)
#' fit_data <- rodpaa(rodpaa_example$view1, rodpaa_example$view2,
#'                    input = "data",
#'                    dissimilarity1 = "euclidean", dissimilarity2 = "euclidean")
#' fit_data$k_hat                                           # 2
"rodpaa_example"

#' IBDMDB paired ileal and rectal microbiome profiles
#'
#' Genus-level 16S rRNA gene profiles of paired baseline ileal and rectal
#' biopsies from 68 participants of the Inflammatory Bowel Disease
#' Multi'omics Database (IBDMDB), part of the Integrative Human Microbiome
#' Project. The two biopsy sites form two views of the same participants.
#'
#' @details
#' The data are paired baseline biopsy profiles, with operational taxonomic
#' units aggregated to genus level. Genera were filtered by prevalence within
#' each site, as in the RoDPAA manuscript. The values are observed read counts,
#' without normalization or transformation.
#'
#' @format A list with four elements, with participants in the same order in
#'   every element:
#' \describe{
#'   \item{ileum}{68 x 99 integer matrix of ileal genus counts.}
#'   \item{rectum}{68 x 94 integer matrix of rectal genus counts.}
#'   \item{diagnosis}{Factor with levels `"Control"` and `"IBD"`.}
#'   \item{taxonomy}{List with the full SILVA taxonomy of the `ileum` and
#'     `rectum` genera, named by the short genus names.}
#' }
#' @source Lloyd-Price J, Arze C, Ananthakrishnan AN, et al. (2019).
#'   Multi-omics of the gut microbial ecosystem in inflammatory bowel
#'   diseases. *Nature*, 569, 655-662. Data available from
#'   <https://ibdmdb.org>.
#' @examples
#' str(ibdmdb_micro[c("ileum", "rectum")])
#'
#' # Bray-Curtis dissimilarity for both sites, rectum as view 1
#' set.seed(1)
#' fit <- rodpaa(ibdmdb_micro$rectum, ibdmdb_micro$ileum, input = "data",
#'               dissimilarity1 = "bray", dissimilarity2 = "bray")
#' fit
"ibdmdb_micro"
