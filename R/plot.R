#' Plot results from a RoDPAA fit
#'
#' Three plots are available:
#' * `"intensity"` (default): the shared intensities (cosines of the principal
#'   angles between the two signal subspaces) by component. Components counted
#'   in the shared dimension \eqn{\hat{k}} are filled, the others are open.
#'   When the permutation test was run, the grey band shows the central 95% of
#'   the permutation null, as a descriptive reference.
#' * `"scores"`: participant scores on two shared factors, optionally
#'   colored by group labels such as diagnosis.
#' * `"pvalue"`: permutation P-values (adjusted, if a correction was used) on
#'   the \eqn{-\log_{10}} scale, with the significance level as a dashed line.
#'
#' Only the first \eqn{\hat{k}} components are statistically significant.
#' Directions extracted beyond \eqn{\hat{k}} with `k_extract` are exploratory
#' and are labeled as such in the scores plot. In the P-value plot, components
#' after the first non-significant one are not counted in \eqn{\hat{k}} even
#' if their P-value is below the significance level, and are drawn in grey.
#'
#' The permutation P-value of component \eqn{j} is the proportion of the
#' \eqn{B} permutations whose \eqn{j}th shared intensity is at least the
#' observed one. Because a P-value of 0 cannot be shown on a log scale, it is
#' drawn at \eqn{1/B} with a triangle, for display only. The P-values stored in
#' the fit are not changed.
#'
#' Colors follow the color-blind-safe Okabe-Ito palette. The title and axis
#' labels can be replaced with `main`, `xlab` and `ylab`.
#'
#' @param x A `"rodpaa"` object from [rodpaa()].
#' @param type `"intensity"` (default), `"scores"` or `"pvalue"`.
#' @param n_show Number of components to show in the `"intensity"` and
#'   `"pvalue"` plots. Defaults to all of them, or \eqn{\hat{k} + 10} when
#'   there are more than 30.
#' @param factors Two shared factors to plot against each other in the
#'   `"scores"` plot.
#' @param group Optional group labels, one per participant in the same order
#'   as the input, used to color the `"scores"` plot.
#' @param col Color of the components counted in \eqn{\hat{k}} (and of the
#'   points in the `"scores"` plot without `group`).
#' @param ... Further graphical arguments passed to [plot()], for example
#'   `main`, `xlab`, `ylab` or `cex.lab`.
#' @return `x`, invisibly.
#' @export
#' @examples
#' set.seed(1)
#' fit <- rodpaa(dist(rodpaa_example$view1), dist(rodpaa_example$view2))
#' plot(fit)                                  # shared intensities
#' plot(fit, type = "scores")                 # SF1 vs SF2
#' plot(fit, type = "pvalue")                 # permutation P-values
plot.rodpaa <- function(x, type = c("intensity", "scores", "pvalue"),
                        n_show = NULL, factors = c(1, 2), group = NULL,
                        col = "#0072B2", ...) {
  type <- match.arg(type)
  if (type == "scores") return(plot_scores(x, factors, group, col, ...))

  s <- x$shared_intensities
  L <- length(s)
  if (L == 0) stop("No shared intensities to plot: a signal rank is 0.")
  k_hat <- x$k_hat
  if (is.null(n_show)) n_show <- if (L > 30) min(L, k_hat + 10) else L
  if (!is.numeric(n_show) || length(n_show) != 1 || !is.finite(n_show) ||
      n_show < 1 || n_show != floor(n_show)) {
    stop("`n_show` must be a positive whole number.")
  }
  idx <- seq_len(min(n_show, L))
  in_k <- idx <= k_hat

  if (type == "intensity") plot_intensity(x, idx, in_k, col, ...) else plot_pvalue(x, idx, in_k, col, ...)
  invisible(x)
}

# ---- shared style ------------------------------------------------------------

# Okabe-Ito palette (color-blind safe) and neutral greys
rodpaa_colors <- list(
  other = "grey55",      # components not counted in k_hat
  line  = "grey45",      # connecting line
  band  = "grey88",      # permutation null
  ref   = "grey40",      # reference lines
  groups = c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000")
)

# Empty plot frame in a clean style (open L-shaped axes, horizontal tick labels).
# Defaults can be overridden by the user through `...` (e.g. main, xlab, ylab).
plot_frame <- function(x, y, defaults, ...) {
  args <- utils::modifyList(
    c(list(x = x, y = y, type = "n", bty = "l", las = 1, tcl = -0.3,
           mgp = c(2.4, 0.6, 0), font.main = 1, cex.main = 1.05), defaults),
    list(...)
  )
  do.call(graphics::plot, args, quote = TRUE)   # quote: keep plotmath labels unevaluated
}

# Legend in the strip above the data: one row when it fits within the plot,
# otherwise two columns
top_legend <- function(labels, ...) {
  args <- list(x = "top", legend = labels, bty = "n", cex = 0.85,
               x.intersp = 0.5, seg.len = 1.2, ...)
  one_row <- do.call(graphics::legend, c(args, horiz = TRUE, plot = FALSE), quote = TRUE)
  if (one_row$rect$w <= diff(graphics::par("usr")[1:2])) {
    do.call(graphics::legend, c(args, horiz = TRUE), quote = TRUE)
  } else {
    do.call(graphics::legend, c(args, ncol = 2), quote = TRUE)
  }
}

# ---- intensity ---------------------------------------------------------------

# Shared intensities with the permutation null band
plot_intensity <- function(x, idx, in_k, col, ...) {
  s <- x$shared_intensities
  k_hat <- x$k_hat
  has_perm <- !is.null(x$permutation)
  pal <- rodpaa_colors

  # y-axis runs to 1 (cosines); the strip above 1 holds the legend, clear of the data
  plot_frame(idx, s[idx], list(
    ylim = c(0, 1.24), xaxt = "n", yaxt = "n",
    xlab = "Principal-angle component",
    ylab = expression("Shared intensity (cos " * theta * ")"),
    main = bquote("Shared dimension " * hat(k) == .(k_hat))), ...)
  graphics::axis(1, at = idx, tcl = -0.3, mgp = c(2.4, 0.6, 0))
  graphics::axis(2, at = seq(0, 1, by = 0.2), las = 1, tcl = -0.3, mgp = c(2.4, 0.6, 0))

  # permutation null, central 95% (descriptive)
  if (has_perm) {
    sp <- x$permutation$sigma_perm[, idx, drop = FALSE]
    lo <- apply(sp, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
    hi <- apply(sp, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
    if (length(idx) == 1) {
      graphics::segments(idx, lo, idx, hi, col = pal$band, lwd = 12, lend = 1)
    } else {
      graphics::polygon(c(idx, rev(idx)), c(lo, rev(hi)), col = pal$band, border = NA)
    }
  }

  if (k_hat > 0 && k_hat < length(idx)) {
    graphics::segments(k_hat + 0.5, 0, k_hat + 0.5, 1.02, lty = 2, col = pal$ref)
  }
  graphics::lines(idx, s[idx], col = pal$line, lwd = 1.2)
  graphics::points(idx, s[idx], pch = 21, cex = 1.1,
                   bg = ifelse(in_k, col, "white"), col = ifelse(in_k, col, pal$other))

  top_legend(c("Shared", "Not shared", if (has_perm) "Permutation null (95%)"),
             pch = c(21, 21, if (has_perm) 15),
             pt.bg = c(col, "white", NA),
             col = c(col, pal$other, if (has_perm) pal$band),
             pt.cex = c(1.1, 1.1, if (has_perm) 1.8))
}

# ---- P-values ----------------------------------------------------------------

# Permutation P-values on the -log10 scale
plot_pvalue <- function(x, idx, in_k, col, ...) {
  perm <- x$permutation
  if (is.null(perm)) {
    stop("This fit has no permutation results (it was fitted with `k` or ",
         "`use_permutation = FALSE`). Refit with the permutation test to plot P-values.")
  }
  pal <- rodpaa_colors
  B <- nrow(perm$sigma_perm)
  adjusted <- !identical(perm$correction_method, "none")
  p <- perm$table$p_adj[idx]

  # plotting only: P = 0 is drawn at 1/B with a triangle
  is_zero <- p == 0
  y <- -log10(ifelse(is_zero, 1 / B, p))
  y_alpha <- -log10(perm$alpha)
  y_top <- max(y, y_alpha)

  plot_frame(idx, y, list(
    ylim = c(0, y_top * 1.4), xaxt = "n",        # headroom for the legend
    xlab = "Principal-angle component",
    ylab = if (adjusted) expression(-log[10] * "(adjusted " * italic(P) * ")") else expression(-log[10] * italic(P)),
    main = paste0("Permutation test (", method_label(perm$correction_method), ", B = ", B, ")")), ...)
  graphics::axis(1, at = idx, tcl = -0.3, mgp = c(2.4, 0.6, 0))
  graphics::segments(idx, 0, idx, y, lwd = 3, lend = 1, col = ifelse(in_k, col, pal$other))
  graphics::abline(h = y_alpha, lty = 2, col = pal$ref)
  if (any(is_zero)) {
    graphics::points(idx[is_zero], y[is_zero], pch = 17, cex = 0.9, col = ifelse(in_k[is_zero], col, pal$other))
  }

  has_zero <- any(is_zero)
  # plotmath (not a Unicode alpha) so the label draws on every graphics device
  labels <- as.expression(c(list("Shared", "Not shared", bquote(alpha == .(perm$alpha))),
                            if (has_zero) list("P = 0, shown at 1/B")))
  top_legend(labels,
             lty = c(1, 1, 2, if (has_zero) NA),
             lwd = c(3, 3, 1, if (has_zero) NA),
             pch = c(NA, NA, NA, if (has_zero) 17),
             col = c(col, pal$other, pal$ref, if (has_zero) "grey30"))
}

# "bonferroni" -> "Bonferroni", "fdr" -> "FDR", ...
method_label <- function(m) {
  switch(m, bonferroni = "Bonferroni", holm = "Holm", fdr = "FDR", none = "no correction", m)
}

# ---- scores ------------------------------------------------------------------

# Participant scores on two shared factors
plot_scores <- function(x, factors, group, col, ...) {
  sf <- x$shared_factors
  k_ext <- ncol(sf)
  if (k_ext < 2) {
    stop("The scores plot needs at least two extracted directions, but this fit has ", k_ext,
         " (k_hat = ", x$k_hat, "). Refit with `k_extract = 2` to plot the two leading ",
         "directions as exploratory.")
  }
  if (!is.numeric(factors) || length(factors) != 2 || anyNA(factors) ||
      any(factors != floor(factors)) || any(factors < 1) || any(factors > k_ext) ||
      factors[1] == factors[2]) {
    stop("`factors` must be two different whole numbers between 1 and ", k_ext, ".")
  }

  if (is.null(group)) {
    pt_col <- col
  } else {
    if (length(group) != nrow(sf) || anyNA(group)) {
      stop("`group` must have one non-missing label per participant (", nrow(sf), ").")
    }
    group <- factor(group)
    ng <- nlevels(group)
    pal <- if (ng <= length(rodpaa_colors$groups)) rodpaa_colors$groups[seq_len(ng)] else
      grDevices::hcl.colors(ng, palette = "Dark 3")
    pt_col <- pal[as.integer(group)]
  }

  # directions beyond k_hat are exploratory, not significant
  lab <- function(j) paste0("Shared factor ", j, if (j > x$k_hat) " (exploratory)")
  xs <- sf[, factors[1]]; ys <- sf[, factors[2]]
  # with a group legend, leave headroom above the points so the legend does not cover them
  ylim <- range(ys)
  if (!is.null(group)) ylim[2] <- ylim[2] + 0.25 * diff(ylim)
  plot_frame(xs, ys, list(xlab = lab(factors[1]), ylab = lab(factors[2]), ylim = ylim,
                          main = "Shared-factor scores"), ...)
  graphics::abline(h = 0, v = 0, lty = 3, col = "grey75")
  graphics::points(xs, ys, pch = 21, cex = 1, col = "white", lwd = 0.5,
                   bg = grDevices::adjustcolor(pt_col, alpha.f = 0.85))
  if (!is.null(group)) {
    graphics::legend("topright", legend = levels(group), pch = 21, pt.bg = pal, col = "white",
                     pt.cex = 1.1, bty = "n", cex = 0.85, inset = 0.02)
  }
  invisible(x)
}
