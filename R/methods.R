# print() and summary() methods for "rodpaa" objects. plot() is in plot.R.
# These only format the output; the fitted object is not changed.

#' @export
print.rodpaa <- function(x, ...) {
  print_overview(x)
  invisible(x)
}

# Overview lines shared by print() and summary()
print_overview <- function(x) {
  d <- x$dimensions
  s <- x$shared_intensities
  row <- function(label, value) cat(sprintf("  %-22s %s\n", label, value))

  cat("RoDPAA fit\n")
  row("Participants:", nrow(x$O_J1))
  row("Signal ranks:", sprintf("r1 = %d, r2 = %d", d$r1, d$r2))
  row("Shared dimension:", sprintf("k = %d (%s)", d$k_hat, rule_label(x$rank_selection$rule)))
  if (!identical(d$k_extracted, d$k_hat)) {
    row("Directions extracted:", sprintf("%d (%d exploratory)", d$k_extracted, d$k_extracted - d$k_hat))
  }
  if (length(s) > 0) {
    m <- min(5, length(s))
    row("Shared intensities:", paste0(paste(sprintf("%.3f", s[seq_len(m)]), collapse = ", "),
                                      if (length(s) > m) ", ..."))
  }
}

#' @export
summary.rodpaa <- function(object, ...) {
  print_overview(object)
  perm <- object$permutation
  if (!is.null(perm)) {
    B <- nrow(perm$sigma_perm)
    cat(sprintf("\nPermutation test (%s correction, alpha = %s, B = %d permutations)\n",
                method_label(perm$correction_method), format(perm$alpha), B))
    tab <- utils::head(perm$table, max(object$k_hat + 3, 5))
    out <- data.frame(
      Component            = tab$component,
      `Shared intensity`   = sprintf("%.3f", tab$singular_value),
      `P-value`            = sprintf("%.3f", tab$p_value),
      `Adjusted P-value`   = sprintf("%.3f", tab$p_adj),
      Shared               = ifelse(seq_len(nrow(tab)) <= object$k_hat, "yes", "no"),
      check.names = FALSE
    )
    print(out, row.names = FALSE, right = TRUE)
  }
  invisible(object)
}

# How the shared dimension was chosen
rule_label <- function(rule) {
  switch(rule, perm = "permutation test", user = "set by user", elbow = "largest drop", rule)
}
