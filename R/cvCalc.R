#' calculate the analytical coefficient of variation from a set of observations
#' 
#' @description calculate the coefficient of variation (CV\%) for all LC-MS features
#' of peak table using a set of  observations (e.g. pooled quality control samples).
#' The CV\% is calculated according to the following equation:
#' 
#' \eqn{CV\% = (\delta / \mu) * 100}
#' 
#' where \eqn{\delta} is the standard deviation and \eqn{\mu} is the mean
#' @param peakTable either a data.frame, full file path as a character string to a  .csv file of a peak table in the form observation (samples) in columns and
#' variables (Mass spectral signals) in rows. If argument is not supplied a GUI (tcltk) file selection window will open and a .csv file can be selected.
#' @param obsNames character vector (i.e. sample/ QC/ Blank) names to identify appropriate observation (sample) columns
#' with to calculate the coefficient of variance (CV\%).
#' @param thresh numeric the minimum CV\% to retain an LC-MS variable (default = 30 i.e. 30\%).
#' 
#' @return a data frame similar to peakTable with any rows (LC-MS features) removed
#' which are above the maximum CV\% threshold. In addition a new column is added
#' 'coeffVar' which contains the CV\% calculated for each variable.
#' @export
cvCalc <- function(peakTable=NULL, obsNames=NULL, thresh=30){
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
  # calc coeff. of var.
  message("calculating coefficient of variation...\n")
  flush.console()
 
  peakTable$coeffVar <- apply(peakTable[, obsIndx], 1, function(Var){
                             (sd(as.numeric(Var))/ mean(as.numeric(Var))) * 100}) 
  # subset and return error if none below thresh
  cvIndx <- which(peakTable$coeffVar < thresh)
  if(length(cvIndx) ==  nrow(peakTable)){
  stop('No features were found to be below the CV% threshold of ', thresh)  
  }
  message(length(cvIndx), ' (', round((length(cvIndx)/ nrow(peakTable)) * 100, digits=1), '%)',
          ' features were below the CV% threshold of ', thresh, '\n')
  flush.console()
  # subset
  peakTable <- peakTable[cvIndx, , drop=F]
  # return peakTable
  return(peakTable)
}