#' substract blanks (i.e. negative controls/ background) from peak table
#' 
#' @description calculate either the mean or median fold change between samples
#' and blanks for every LC-MS variable and retain only LC-MS features above
#' the blank fold change cut-off.
#' 
#' @param peakTable either a data.frame, full file path as a character string to a  .csv file of a peak table in the form observation (samples) in columns and
#' variables (Mass spectral signals) in rows. If argument is not supplied a GUI file selection window will open and a .csv file can be selected.
#' @param blankNames character vector of blank (i.e. negative control) names to identify
#' appropriate observation (blank) columns.
#' @param sampNames character vector of sample names to identify appropriate observation (sample) columns.
#' @param method character either 'median' or 'mean' fold change calculation.
#' @param thresh numeric sample:blank fold change cut-off. Any LC-MS
#' features below this threshold will be removed. (default = 2). 
#' 
#' @return a data frame similar to peakTable with any rows (LC-MS features) removed
#' which are below the minimum blank fold change threshold.
#' @export
blankSub <- function(peakTable=NULL, blankNames=NULL, sampNames=NULL, method="mean",
                     thresh=2){
  # error handling
  if(is.null(blankNames)){
    stop('argument blankNames is missing with no default')
  } else if(is.null(sampNames)){
    stop('argument sampNames is missing with no default')
  } else if(!(method == 'mean' | method == 'median')){
    stop("argument method must be one of either 'mean' or 'median'")
  }
  # error handling or read from csv function
  peakTable <- tableCheckRead(peakTable, stringsAsFactors=F)
  # match blankNames to peak table colnames
  blankIndx <- match(blankNames, colnames(peakTable))
  # if none matched then stop
  if(length(blankIndx) == 0){
    stop("none of the blankNames were matched to column names of the peakTable argument. Please check and try again.")
  }
  # match sampNames to peak table colnames
  sampIndx <- match(sampNames, colnames(peakTable))
  # if none matched then stop
  if(length(sampIndx) == 0){
    stop("none of the sampNames were matched to column names of the peakTable argument. Please check and try again.")
  }
  # calculate fold change and create new column
  message("performing blank substraction using...\n", length(sampIndx), " samples\n",
          length(blankIndx), " blanks\n")
  flush.console()
  
  if(method == 'mean'){
  peakTable[ ,'meanFCsampBlank'] <- apply(peakTable, 1, function(Var){
  mean(as.numeric(Var[sampIndx])) / mean(as.numeric(Var[blankIndx]))})
  # above thresh indx
  abThrIndx <- which(peakTable[ ,'meanFCsampBlank'] > thresh) 
  message(length(abThrIndx), ' (', round((length(abThrIndx)/nrow(peakTable)) * 100, digits=1), 
          '%) of the LC-MS features were above the mean sample:blank fold change threshold of ', 
          thresh, '\n')
  flush.console()
  # error handling
  if(length(abThrIndx) == 0){
    stop('No LC-MS features were found to be above the mean sample:blank fold change threshold of ', 
         thresh, '\n')
  }
  # subset table
  peakTable <- peakTable[abThrIndx, ]
  
  } else {
  peakTable[ ,'medianFCsampBlank'] <- apply(peakTable, 1, function(Var){
  median(as.numeric(Var[sampIndx])) / median(as.numeric(Var[blankIndx]))})
  # above thresh indx
  abThrIndx <- which(peakTable[ ,'medianFCsampBlank'] > thresh) 
  message(length(abThrIndx), ' (', round((length(abThrIndx)/nrow(peakTable)) * 100, digits=1), 
          '%) of the LC-MS features were above the median sample:blank fold change threshold of ', 
          thresh, '\n')
  flush.console()
  # error handling
  if(length(abThrIndx) == 0){
    stop('No LC-MS features were found to be above the median sample:blank fold change threshold of ', 
         thresh, '\n')
  }
  # subset table
  peakTable <- peakTable[abThrIndx, ]
  }
  # return peakTable
  return(peakTable)
}