#' plot PCA with expanded Hotelling's ellipse adapted from \code{\link{plotPcs}} from 
#' the \code{\link{pcaMethods}} package.
#'  
#' @param object \code{\link{pcaRes}} class object
#' @param exEllipse matrix of expanded Hotelling's Ellipse see \code{\link{pcaOutId}}.
#' @param type character either "scores" or "loadings" for scores or loadings plot respectively
#' @param Further arguments to \code{\link{pairs}} on which this function is based.
#'  
#' @seealso \code{\link{pcaOutId}}, \code{\link{pcaRes}}, \code{\link{pca}}.
#' @export
plotPcsEx <- function(object, exEllipse,  type = c("scores", "loadings"), ...){
  pcs <- 1:nP(object)
  type <- match.arg(type)
  panel <- function(x, y, ...) {
    abline(h = 0, v = 0, col = "black")
      A <- length(pcs)
      el <- hotEllipse(x, y)
      lines(el)
      # expand ellipse
      lines(exEllipse, col="red")
      points(x, y, ...)
    }
  switch(type, scores = {
    labels <- paste("PC", pcs, "\n", "R^2 =", round(object@R2[pcs], 
                                                    2))
    pairs(scores(object)[, pcs], labels = labels, panel = panel, 
          upper.panel = NULL, ...)
  }, loadings = {
    if (method(object) == "nlpca") stop("Loadings plot not applicable for non-linear PCA")
    labels <- paste("PC", pcs, "\n", "R^2 =", round(object@R2[pcs], 
                                                    2))
    pairs(loadings(object)[, pcs], labels = labels, panel = panel, 
          upper.panel = NULL, ...)
  })
}