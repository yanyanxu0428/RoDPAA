#' Simulated two-view example data
#'
#' A small simulated dataset for trying out [rodpaa()]. Two views are measured
#' on the same 80 participants. Both views contain two shared latent factors,
#' each view also contains one view-specific factor, and both contain
#' independent noise. RoDPAA should therefore find two shared factors.
#'
#' @details
#' For \eqn{n = 80} participants and views \eqn{v = 1, 2} with
#' \eqn{p_1 = 50} and \eqn{p_2 = 40} features, the data were generated as
#' \deqn{X_v = Z B_v + 1.5\, s_v a_v^\top + E_v,}
#' where
#' * \eqn{Z} (\eqn{n \times 2}) contains the two shared factors, common to
#'   both views (returned as `shared`);
#' * \eqn{B_v} (\eqn{2 \times p_v}) contains the loadings of the shared
#'   factors in view \eqn{v};
#' * \eqn{s_v} (\eqn{n \times 1}) is a factor specific to view \eqn{v}, with
#'   loadings \eqn{a_v} (\eqn{p_v \times 1});
#' * \eqn{E_v} (\eqn{n \times p_v}) is noise.
#'
#' All entries of \eqn{Z}, \eqn{B_v}, \eqn{s_v}, \eqn{a_v} and \eqn{E_v} are
#' independent standard normal, and \eqn{s_1} and \eqn{s_2} are independent.
#' Each view therefore has signal rank 3 (two shared factors and one
#' view-specific factor), and the two views share two dimensions.
#' The data were generated with `set.seed(2026)`. The full code is in
#' `data-raw/rodpaa_example.R` in the package source on GitHub.
#'
#' @format A list with three matrices, all with participants `"P01"` to
#'   `"P80"` in rows, in the same order:
#' \describe{
#'   \item{view1}{80 x 50 feature matrix of view 1.}
#'   \item{view2}{80 x 40 feature matrix of view 2.}
#'   \item{shared}{80 x 2 matrix of the true shared factors \eqn{Z}, for
#'     checking the result.}
#' }
#' @source Simulated, see Details.
#' @examples
#' str(rodpaa_example)
#'
#' # RoDPAA recovers the two shared factors
#' set.seed(1)
#' fit <- rodpaa(dist(rodpaa_example$view1), dist(rodpaa_example$view2))
#' fit$k_hat                                                # 2
#' summary(lm(fit$shared_factors[, 1] ~ rodpaa_example$shared))$r.squared
"rodpaa_example"
