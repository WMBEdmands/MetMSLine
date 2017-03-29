#' Identify outliers in a PCA plot based on expansion of Hotelling's ellipse.
#'  
#' @description identify potential analytical/preparative outliers on a pca scores plot based
#' on a proportional expansion of the Hotelling's T2 ellipse using the \code{\link{pcaMethods}} package.
#' 
#' @param peakTable either a data.frame, full file path as a character string to a  .csv file of a peak table in the form observation (samples) in columns and
#' variables (Mass spectral signals) in rows. If argument is not supplied a GUI file selection window will open and a .csv file can be selected.
#' @param obsNames character vector of observation (i.e. sample/ QC/ Blank) names to identify appropriate observation (sample) columns.
#' @param outTol proportional expansion value for Hotelling's ellipse (PC1 and PC2), any outlying samples beyond this will be removed. 
#' default = 1.2
#' @param maxIter number of iterations of pca model calculation/ outlier removal to perform. 
#' The iteration process will stop if no further outliers are detected and the 
#' last PCA model calculated. default = 2
#' @param ... additional arguments to \code{\link{pca}}. 
#' 
#' @details principal components analysis is a commonly used method to identify
#' clustering and potentially outlying samples in multivariate datasets. The
#' Hotellings student's T2 distribution is also commonly used to detect strongly
#' outlying samples on the resulting scores plot. Weakly outlying samples (perhaps
#' biological rather than analytical) in origin will appear close to the Hotelling's
#' ellipse therefore a small proportional expansion of the ellipse can be used to
#' automatically remove strongly outlying samples whilst retaining weakly outlying 
#' samples. The first two principal components (PCs) represent the greatest sources
#' of systematic variation in multivariate datasets therefore considered of the first
#' two PCs should be sufficient to identify strong analytical outliers, representing
#' such things as, strongly contaminated biological samples, failed LC autosampler
#' injections, experimental preparation errors and temporal errors in mass 
#' spectrometer performance etc.
#' 
#' @return a list containing two elements the original peak table with outliers
#' removed and a nested list with a sub list for each iteration of pca calculation,
#' and outlier identification. Each pca result sublist consists of 3 elements:
#' 
#' 1. a \code{\link{pcaRes}} object of the pca calculation.
#' 
#' 2. the coordinates of the expanded hotelling's ellipse.
#' 
#' 3. a named logical vector of each outlier detected.
#' 
#' @seealso \code{\link{pca}}, \code{\link{pcaRes}}
#' @export
pcaOutId <- function(peakTable=NULL, obsNames=NULL, outTol=1.2, maxIter=2, ...){
  
  if(is.null(obsNames)){
    stop('argument obsNames is missing with no default')
  } 
  # error handling or read from csv function
  peakTable <- tableCheckRead(peakTable, stringsAsFactors=F)
 
  # match obsNames to peak table colnames
  obsIndx <- match(obsNames, colnames(peakTable))
  # if less than all matched then stop
  if(length(obsIndx) < length(obsNames)){
    stop(length(obsIndx), " of ", length(obsNames), 
         " observation names were matched in the peakTable column names, check the obsNames and peakTable column names")
  }
  # subset table
  obsTable <- peakTable[, obsIndx]
  # check if any zeros or NAs 
  if(any(is.na(obsTable) | obsTable == 0)){
    stop("peakTable observations contain NAs or zeros use the ?zeroFill function")
  }
  # n iterations counter for while loop
  nIter <- 1
  # create empty pcaResults list
  pcaResults <- vector("list", maxIter + 1)
  # calc first Pca model
  message("Calculating PCA model ", nIter, "...")
  flush.console()
  # calc Pca Model
  pcaResult <- pcaMethods::pca(t(obsTable), nPcs=2, ...)
  # id poss outliers
  possOutliers <- pcaResOut(pcaResult, outTol)
  pcaResults[[nIter]] <- possOutliers
  # if any outliers after first then continue
  if(any(possOutliers$possOut)){
    # remove outliers
    outRem <- obsTable[, possOutliers$possOut == F]
    message(sum(possOutliers$possOut), " outliers identified PCA model ", nIter)
    flush.console()
    # while any poss outliers or n iterations continue to remove outliers 
    while(any(possOutliers$possOut) & nIter < maxIter){
      nIter <- nIter + 1
      message("Calculating PCA model ", nIter, "...")
      flush.console()
      # calc Pca Model
      pcaResult <- pcaMethods::pca(t(outRem), nPcs=2, ...)
      # id poss outliers
      possOutliers <- pcaResOut(pcaResult, outTol)
      pcaResults[[nIter]] <- possOutliers
      # remove any outliers
      if(any(possOutliers$possOut)){
      message(sum(possOutliers$possOut), " outliers identified PCA model ", nIter)
      flush.console()
      outRem <- outRem[, possOutliers$possOut == F]
      }
     } # end while loop
    # if necessary calculate last PCA model
    if(any(possOutliers$possOut)){
    message("Calculating last PCA model ", nIter + 1, "...")
    flush.console()
    # calc Pca Model
    pcaResult <- pcaMethods::pca(t(outRem), nPcs=2, ...)
    # id poss outliers
    possOutliers <- pcaResOut(pcaResult, outTol)
    possOutliers$possOut[possOutliers$possOut == T] <- F 
    pcaResults[[nIter + 1]] <- possOutliers
    }
   } else {
     message("No outliers identified PCA model ", nIter)
     flush.console()
   }
  # remove any empty list elements
  pcaResults <- pcaResults[sapply(pcaResults, length) == 3]
  # extract names non outlying samples
  nameNonOut <- pcaResults[[length(pcaResults)]]$possOut
  nameNonOut <- names(nameNonOut)[nameNonOut == F]
  #outlier removed index
  outRemIndx <- c(seq(1, ncol(peakTable), 1)[-obsIndx], obsIndx[match(nameNonOut, obsNames)])
  # return pca Results and outliers removed peakTable
  return(list(outRem=peakTable[, outRemIndx], pcaResults=pcaResults))
}