#' pooled quality control based signal drift/ attenuation smoothing
#' 
#' @description attempts to correct for systematic drifts in LC-MS datasets using
#' intervally injected pooled quality control samples representative of all samples
#' of a experiment.
#' 
#' @param peakTable either a data.frame, full file path as a character string to a  .csv file of a peak table in the form observation (samples) in columns and
#' variables (Mass spectral signals) in rows. If argument is not supplied a GUI file selection window will open and a .csv file can be selected.
#' @param sampNames character vector of sample names to identify appropriate observation (sample) columns.
#' @param qcNames character vector of quality control (QC) names to identify appropriate observation (sample) columns.
#' @param smoothSpan smoother span, proportion of total data points used for smoothing. See \code{\link{lowess}}. default = 0.2
#' @param ... additional arguments to \code{\link{lowess}}. 
#' 
#' @details systematic drifts in mass spectral signal can be corrected using
#' intervally injected pooled quality control samples for each mass spectral signal
#' variable in the peak table. Sample intensities are first replaced
#' with NAs and cubic spline interpolation (\code{\link{na.spline}} from the zoo package)
#' is used to replace missing values between QC samples. The \code{\link{lowess}} function
#' is then used to model systematic drift/ attenuation within each mass spectral
#' signal variable. The degree of smoothing determined by the smoothSpan argument. 
#'  
#' 
#' @return a data frame identical to peakTable with signal drift/ attenuation adjusted.
#' @source \url{http://www.nature.com/nprot/journal/v6/n7/abs/nprot.2011.335.html} 
#' Procedures for large-scale metabolic profiling of serum and plasma using gas 
#' chromatography and liquid chromatography coupled to mass spectrometry
#' @seealso \code{\link{lowess}}, \code{\link{na.spline}}.
#' @export
lowessSmooth <- function(peakTable=NULL, sampNames=NULL, qcNames=NULL, 
                         smoothSpan=0.2, ...){
  
  if(is.null(sampNames)){
    stop('argument sampNames is missing with no default')
  } else if(is.null(qcNames)){
    stop('argument qcNames is missing with no default')
  } 
  # error handling or read from csv function
  peakTable <- tableCheckRead(peakTable, stringsAsFactors=F)
  # match sample names to peak table colnames
  sampIndx <- match(sampNames, colnames(peakTable))
  # if less than all matched then stop
  if(length(sampIndx) < length(sampNames)){
    stop(length(sampIndx), " of ", length(sampNames), 
         " observation names were matched in the peakTable column names, check the sampNames and peakTable column names")
  }
  # match QC names to peak table colnames
  qcIndx <- match(qcNames, colnames(peakTable))
  # if less than all matched then stop
  if(length(qcIndx) < length(qcNames)){
    stop(length(qcIndx), " of ", length(qcNames), 
         " observation names were matched in the peakTable column names, check the qcNames and peakTable column names")
  }
  
  # subset table
  obsTable <- peakTable[, sort(c(sampIndx, qcIndx))]
  # check if any zeros or NAs 
  if(any(is.na(obsTable) | obsTable == 0)){
    stop("peakTable observations contain NAs or zeros use the ?zeroFill function")
  }
  # match sample names to peak table colnames
  sampIndxObs <- match(sampNames, colnames(obsTable))
  # match QC names to peak table colnames
  qcIndxObs <- match(qcNames, colnames(obsTable))
  
  message("applying LOWESS signal drift/ attenuation smoothing...")
  flush.console()
  # convert to numeric 
  obsTable <- apply(obsTable, 2, as.numeric)
  
  ##replace sample columns intensities with NA
  obsTable[, sampIndxObs] <- NA
  spline <- t(apply(obsTable, 1, zoo::na.spline)) #cubic spline interpolation for missing values
  #apply Lowess smoothing on QCs
  Lowess <- t(apply(spline, 1, function(Var) lowess(Var, f=smoothSpan, ...)$y)) 
  ###identify minimum values of each QC variable to replace negative or less than half minimum   Lowess values 
  MinQC <- apply(obsTable, 1, function(x) min(x[is.na(x)==F]))
  NegLowessIndx <- which(apply(Lowess, 1, function(x) any(x < 0)))
  ###replace negative Lowess values with half minimum QC sample if at the end of run order sequence
  NegLowessEndStart <- rep(0, nrow(Lowess))
  for(j in NegLowessIndx)
  {
    ##if only one to 3 Lowess values at the end or start of the injection sequence are less than zero then fill with next nearest positive value
    NegLowessEndStart[j] <- (length(which(Lowess[j, ] < MinQC[j]/2))<4 & 
                               any(which(Lowess[j, ] < MinQC[j]/2) == ncol(Lowess)) | any(which(Lowess[j, ] < MinQC[j]/2) == 1)) * 1
    Lowess[j, which(Lowess[j, ] < MinQC[j]/2)] <- MinQC[j]/2
  }    
  ###reinsert normalised or non-normalised samples###
  median.curve.df <- as.matrix(apply(Lowess, 1, median))
  median.curve.df <- median.curve.df[, rep(1, ncol(Lowess))]
  curve.prop <- median.curve.df/ Lowess
  obsTable <- apply(peakTable[, sort(c(sampIndx, qcIndx))], 2, as.numeric)
  obsTable <- obsTable * curve.prop
  
  
  # replace with normalized obsTable
  peakTable[, sort(c(sampIndx, qcIndx))] <- obsTable
  # return zero filled table
  return(peakTable)
}