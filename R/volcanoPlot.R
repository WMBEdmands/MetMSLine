#' Volcano plot slightly modified from \code{\link{VolcanoPlot}} from the metabolomics
#' package.
#' 
#' @param folds A vector of fold changes with metabolite names.
#' @param pvals A vector of corresponding p-values with metabolite names.
#' @param cexcutoff	Font size of the cut-off labels.
#' @param cexlab	Font size of the variable labels.
#' @param plimit A numeric indicating the p value cutoff. The default is set to 0.05.
#' @param fclimit	A numeric indicating the lower fold cutoff. The default is set to 2.
#' @param xlab x-axis label.
#' @param ylab y-axis label
#' @param main Plot title.
#' @param ...	Other graphical parameters. See par.
#' @export
volcanoPlot <- function (folds, pvals, cexcutoff = 0.7, cexlab = 0.5, plimit = 0.05, 
                         fclimit = 2, xlab = "log2 Fold Change", ylab = "-log10 t-Test P-value", 
                         main = "Volcano Plot", ...) 
{
  x_min <- (-1.0)
  x_max <- 1.0
  if (min(range(folds, finite = TRUE)) <= x_min) {
    x_min <- min(range(folds, finite = TRUE))
  }
  if (max(range(folds, finite = TRUE)) >= x_max) {
    x_max <- max(range(folds, finite = TRUE))
  }
  x_range <- c(x_min, x_max)
  y_min <- 0
  y_max <- 2
  if (min(range(-log10(pvals), finite = TRUE)) <= y_min) {
    y_min <- min(range(-log10(pvals), finite = TRUE))
  }
  if (max(range(-log10(pvals), finite = TRUE)) >= y_max) {
    y_max <- max(range(-log10(pvals), finite = TRUE))
  }
  y_range <- c(y_min, y_max)
  plot(x_range, y_range, type = "n", xlab = xlab, ylab = ylab, 
       main = main, ...)
  abline(h = -log10(plimit), col = "green", lty = "44")
  mtext(paste("pval =", plimit), side = 2, at = -log10(plimit), 
        cex = cexcutoff, las = 1)
  abline(v = c(-log2(fclimit), log2(fclimit)), col = "violet", 
         lty = "1343")
  mtext(c(paste("-", fclimit, "fold"), paste("+", fclimit, 
                                             "fold")), side = 3, at = c(log2(1/fclimit), log2(fclimit)), 
        cex = cexcutoff, las = 1)
  for (ii in 1:length(pvals)) {
    if (-log10(pvals[ii]) > (-log10(plimit))) {
      if (folds[ii] > (-log2(fclimit))) {
        if (folds[ii] < log2(fclimit)) {
          points(folds[ii], -log10(pvals[ii]), col = "orange", 
                 pch = 20)
        }
        else {
          points(folds[ii], -log10(pvals[ii]), col = "red", 
                 pch = 20)
          text(folds[ii], -log10(pvals[ii]), labels = names(folds)[ii], 
               pos = if (-log10(pvals[ii]) < 0.95 * max(y_range)) {
                 if (folds[ii] < 0.75 * max(x_range)) {
                   4
                 }
                 else {
                   2
                 }
               }
               else {
                 1
               }, cex = cexlab)
        }
      }
      else {
        points(folds[ii], -log10(pvals[ii]), col = "blue", 
               pch = 20)
        text(folds[ii], -log10(pvals[ii]), labels = names(folds)[ii], 
             pos = if (-log10(pvals[ii]) < 0.95 * max(y_range)) {
               if (folds[ii] < 0.75 * max(x_range)) {
                 4
               }
               else {
                 2
               }
             }
             else {
               1
             }, cex = cexlab)
      }
    }
    else {
      points(folds[ii], -log10(pvals[ii]), col = "purple", 
             pch = 20)
    }
  }
} 