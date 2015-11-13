#' zero fill peak table
#' 
#' @description identify and replace zero/ missing values (NA) in peak table with half
#' the mimimum non-zero observed value.
#' 
#' @param peakTable either a data.frame, full file path as a character string to a  .csv file of a peak table in the form observation (samples) in columns and
#' variables (Mass spectral signals) in rows. If argument is not supplied a GUI file selection window will open and a .csv file can be selected.
#' @param obsNames character vector of observation (i.e. sample/ QC/ Blank) names to identify appropriate observation (sample) columns.
#' @param value numeric value to fill zero/ missing values (NA). By default
#' half the mimimum non-zero observed peak intensity is used.
#' 
#' @return a data frame identical to peakTable with zero/ NAs filled with half
#' the mimimum non-zero observed value.
#' @export
zeroFill <- function(peakTable=NULL, obsNames=NULL, value=NULL){
  # error handling
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
  # replace any NAs with zero
  obsTable <- peakTable[, obsIndx]
  obsTable[is.na(obsTable)] <- 0
  # replace any zero with half the minimum not zero value
  message("zero filling with ", 
          ifelse(is.null(value), 'half the minimum non-zero value', value), '\n')
  flush.console()
  if(is.null(value)){
  minNotzero <- sort(unlist(obsTable))
  obsTable[obsTable <= 2] <- (minNotzero[minNotzero >= 2][1])/2
  } else {
  obsTable[obsTable == 0] <- value 
  }
  # replace zero filled with obsTable
  peakTable[, obsIndx] <- obsTable
  # return zero filled table
  return(peakTable)
}