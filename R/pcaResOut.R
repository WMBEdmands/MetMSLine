#' expand hotellings ellipse from a pcaRes object \code{\link{pcaRes}}
#' and expand hotellings ellipse to identify outliers
#' @param pcaResult a \code{\link{pcaRes}} object. 
#' @param outTol proportional expansion value for Hotelling's ellipse (PC1 and PC2), 
#' any outlying samples beyond this will be removed. 
#' @return list containing pca results file, and named logical vector of possible
#' outliers.
pcaResOut <- function(pcaResult, outTol){
  # expand hotellings ellipse proportionately
  exHotEllipse <- hotEllipse(pcaResult@scores[, 1], pcaResult@scores[, 2])
  # logical negs
  negIndx.tmp <- exHotEllipse < 0
  exHotEllipse <- abs(exHotEllipse) * outTol
##HotE(scores[,2],scores[,3])))*out.tol
  # id poss outliers beyond hotellings ellipse expansion return logical
  possOut.tmp <- apply(pcaResult@scores, 1, function(Samp){
                       outPC1 <- which(exHotEllipse[, 1] <= abs(Samp[1]))
                       outPC1PC2 <- any(exHotEllipse[outPC1, 2] <= abs(Samp[2]))
                       return(outPC1PC2)})
# replace neg values
exHotEllipse[negIndx.tmp[, 1], 1] <- exHotEllipse[negIndx.tmp[, 1], 1]  * -1
exHotEllipse[negIndx.tmp[, 2], 2] <- exHotEllipse[negIndx.tmp[, 2], 2]  * -1
return(list(pcaResult=pcaResult, exHotEllipse=exHotEllipse, possOut=possOut.tmp))
}