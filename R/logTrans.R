#' log transform peak table
#' 
#' @description log transform zero filled peak table to correct for multiplicative noise.
#' 
#' @param peakTable either a data.frame, full file path as a character string to a  .csv file of a peak table in the form observation (samples) in columns and
#' variables (Mass spectral signals) in rows. If argument is not supplied a GUI file selection window will open and a .csv file can be selected.
#' @param obsNames character vector of observation (i.e. sample/ QC/ Blank) names to identify appropriate observation (sample) columns.
#' @param base the base with respect to which logarithms are computed \code{\link{log}}. Defaults to e=exp(1)
#' 
#' @return a data frame identical to peakTable with observations log transformed.
#' @export
logTrans <- function(peakTable=NULL, obsNames=NULL, base=exp(1)){
  
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
    stop("peakTable observations contain NAs or zeros unable to log transform, use the ?zeroFill function")
  }
  message("log transforming to the base ", round(base, digits=3), "...\n")
  flush.console()
  # convert chars to numeric
  obsTable <- apply(obsTable, 2, as.numeric)
  obsTable <- apply(obsTable, 1:2, log, base)
  # replace log transformed with obsTable
  peakTable[, obsIndx] <- obsTable
  # return zero filled table
  return(peakTable)
}