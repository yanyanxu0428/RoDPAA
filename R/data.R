#' Simulated two-view example data
#'
#' A small simulated dataset for trying out [rodpaa()]. Two views are measured
#' on the same 80 participants. Both views contain two shared latent factors,
#' each view also contains one view-specific factor, and both contain
#' independent noise. RoDPAA should therefore find two shared factors.
#'
#' @format A list with three matrices, all with participants `"P01"` to
#'   `"P80"` in rows, in the same order:
#' \describe{
#'   \item{view1}{80 x 50 feature matrix of view 1.}
#'   \item{view2}{80 x 40 feature matrix of view 2.}
#'   \item{shared}{80 x 2 matrix of the true shared factors, for checking the
#'     result.}
#' }
#' @source Simulated. The code is in `data-raw/rodpaa_example.R` in the
#'   package source on GitHub.
#' @examples
#' str(rodpaa_example)
"rodpaa_example"
