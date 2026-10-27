# print() and summary() methods for "rodpaa" objects. plot() is in plot.R.

#' @export
print.rodpaa <- function(x, ...) {
  d <- x$dimensions
  cat("RoDPAA fit\n")
  cat("  participants:            ", nrow(x$O_J1), "\n")
  cat("  signal ranks (r1, r2):   ", d$r1, ",", d$r2, "\n")
  cat("  shared dimension (k_hat):", d$k_hat, "\n")
  if (!identical(d$k_extracted, d$k_hat)) {
    cat("  directions extracted:    ", d$k_extracted, "\n")
  }
  if (length(x$shared_intensities) > 0) {
    m <- min(5, length(x$shared_intensities))
    cat("  leading shared intensities:", paste(round(x$shared_intensities[1:m], 3), collapse = ", "), "\n")
  }
  invisible(x)
}

#' @export
summary.rodpaa <- function(object, ...) {
  print(object)
  if (!is.null(object$permutation)) {
    cat("\nPermutation test (", object$permutation$correction_method,
        " correction, alpha = ", object$permutation$alpha, "):\n", sep = "")
    tab <- object$permutation$table
    print(utils::head(tab, max(object$k_hat + 3, 5)), row.names = FALSE)
  }
  invisible(object)
}
