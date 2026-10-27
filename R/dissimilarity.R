# Step 1 inputs: data matrices -> dissimilarity matrices, and input checks.
# Nothing here transforms, normalizes or reorders the data.

#' Check and convert one data view to a numeric matrix
#'
#' @param x A data matrix or data frame, participants in rows.
#' @param name Argument name, for error messages.
#' @return A numeric matrix.
#' @keywords internal
#' @noRd
as_data_matrix <- function(x, name) {
  if (inherits(x, "dist")) {
    stop("`", name, "` is a dissimilarity (`dist`) object. Use `input = \"dissimilarity\"`.")
  }
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) {
      stop("`", name, "` has non-numeric columns. Keep only the numeric features, ",
           "with participant IDs as row names.")
    }
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x)) stop("`", name, "` must be a numeric matrix or data frame.")
  if (any(!is.finite(x))) stop("`", name, "` contains missing or infinite values.")
  if (nrow(x) < 2) stop("`", name, "` must have at least two participants (rows).")
  # a square, symmetric, zero-diagonal matrix may still be valid data, so only warn
  if (looks_like_dissimilarity(x)) {
    warning("`", name, "` looks like a dissimilarity matrix (square, symmetric, zero diagonal). ",
            "It is analysed as a data matrix, as requested by `input = \"data\"`. ",
            "If it is a dissimilarity matrix, use `input = \"dissimilarity\"`.", call. = FALSE)
  }
  x
}

#' Bray-Curtis dissimilarity between the rows of a data matrix
#'
#' \eqn{d_{ij} = \sum_g |x_{ig} - x_{jg}| / \sum_g (x_{ig} + x_{jg})}, the same
#' as `vegan::vegdist(x, "bray")`.
#'
#' @param x Nonnegative numeric matrix, participants in rows, no all-zero rows.
#' @return A `dist` object.
#' @keywords internal
#' @noRd
bray_curtis <- function(x) {
  n <- nrow(x)
  rs <- rowSums(x)
  D <- matrix(0, n, n, dimnames = list(rownames(x), rownames(x)))
  for (i in seq_len(n - 1)) {
    j <- (i + 1):n
    num <- colSums(abs(t(x[j, , drop = FALSE]) - x[i, ]))
    D[i, j] <- num / (rs[i] + rs[j])
  }
  stats::as.dist(D + t(D))
}

#' Dissimilarity matrix of one data view
#'
#' @param x Numeric data matrix, participants in rows.
#' @param method `"euclidean"` or `"bray"`.
#' @param name Argument name, for error messages.
#' @return A `dist` object.
#' @keywords internal
#' @noRd
compute_dissimilarity <- function(x, method, name) {
  if (method == "euclidean") return(stats::dist(x))
  if (any(x < 0)) {
    stop("Bray-Curtis needs nonnegative values (counts or relative abundances), but `", name,
         "` has negative values, for example after a CLR or log transform. ",
         "Use \"euclidean\" for transformed data.")
  }
  if (any(rowSums(x) == 0)) {
    stop("Bray-Curtis is undefined for participants whose values are all zero: `", name,
         "` has ", sum(rowSums(x) == 0), " such row(s).")
  }
  bray_curtis(x)
}

#' Is a matrix square, symmetric, nonnegative with a zero diagonal?
#'
#' Returns `FALSE` (never an error) for anything that is not a finite numeric
#' square matrix with at least two rows.
#' @keywords internal
#' @noRd
looks_like_dissimilarity <- function(D) {
  D <- tryCatch(as.matrix(D), error = function(e) NULL)
  if (is.null(D) || !is.numeric(D) || length(D) == 0 || length(dim(D)) != 2 ||
      nrow(D) != ncol(D) || nrow(D) < 2 || any(!is.finite(D))) {
    return(FALSE)
  }
  tol <- sqrt(.Machine$double.eps) * max(1, max(abs(D)))
  nrow(D) == ncol(D) && nrow(D) > 1 &&
    isSymmetric(unname(D), tol = sqrt(.Machine$double.eps)) &&
    all(abs(diag(D)) <= tol) && min(D) >= -tol
}

#' Check a dissimilarity input
#'
#' @param D A dissimilarity matrix or `dist` object.
#' @param name Argument name, for error messages.
#' @return `D`, unchanged.
#' @keywords internal
#' @noRd
check_dissimilarity <- function(D, name) {
  # `dist` objects are checked too: they can hold missing or negative values
  M <- tryCatch(as.matrix(D), error = function(e) NULL)
  if (is.null(M) || !is.numeric(M) || length(dim(M)) != 2) {
    stop("`", name, "` must be a numeric dissimilarity matrix or `dist` object.")
  }
  if (any(!is.finite(M))) stop("`", name, "` must be numeric with no missing or infinite values.")
  if (nrow(M) != ncol(M)) {
    stop("`", name, "` is not square, so it is not a dissimilarity matrix. ",
         "If it is a data matrix (participants x features), use `input = \"data\"`.")
  }
  if (nrow(M) < 2) stop("`", name, "` must cover at least two participants.")
  if (!looks_like_dissimilarity(M)) {
    if (isSymmetric(unname(M), tol = sqrt(.Machine$double.eps)) && all(diag(M) > 0)) {
      stop("`", name, "` looks like a kernel (similarity) matrix, not a dissimilarity matrix: ",
           "its diagonal is not zero. Use `input = \"kernel\"`.")
    }
    stop("`", name, "` is not a dissimilarity matrix: it must be symmetric, nonnegative and ",
         "have a zero diagonal. If it is a data matrix, use `input = \"data\"`.")
  }
  D
}
