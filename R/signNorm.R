#' Normalize mass spectral data peak table
#' 
#' @description normalizes mass spectral signal either by the median fold change 
#' (probabilistic quotient) or total ion signal methods
#' 
#' @param peakTable either a data.frame, full file path as a character string to a  .csv file of a peak table in the form observation (samples) in columns and
#' variables (Mass spectral signals) in rows. If argument is not supplied a GUI file selection window will open and a .csv file can be selected.
#' @param obsNames character vector of observation (i.e. sample/ QC/ Blank) names to identify appropriate observation (sample) columns.
#' @param method either "medFC" for median fold change or "totIon" for total ion signal normalization.
#' also a custom vector of factors equal in length to the obsNames argument
#'  with which to normalize the data can also be supplied. default = "medFC".
#' 
#' @return a data frame identical to peakTable with samples normalized according to specificied method.
#' @source \url{http://www.ncbi.nlm.nih.gov/pubmed/21526840} 
#' Optimized preprocessing of ultra-performance liquid chromatography/mass spectrometry 
#' urinary metabolic profiles for improved information recovery.
#' @export
signNorm <- function(peakTable=NULL, obsNames=NULL, method="medFC"){
  
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
  # convert to numeric 
  obsTable <- apply(obsTable, 2, as.numeric)
  
  # if custom vector of factors
  if(length(method) > 1){
  # error handling
  if(length(obsNames) != length(method)){
    stop("The length of custom normalization factor must be the same as the obsNames argument")
  }
  message("custom normalization method...\n")
  flush.console()
  # sweep out the custom vector 
  obsTable <- sweep(obsTable, 2, method, "/")
  } else { # undertake default methods  
  message(ifelse(method == "medFC", "Median fold change", "Total ion"), 
          " normalization...\n")
  flush.console()

  if(method == "medFC"){
    
    normalize.medFC <- function(mat) {
      # Perform median fold change normalisation
      #           X - data set [Variables & Samples]
      medSam <- apply(mat, 1, median)
      medSam[which(medSam==0)] <- 0.0001
      mat <- apply(mat, 2, function(mat, medSam){
        medFDiSmpl <- mat/medSam
        vec <- mat/median(medFDiSmpl)
        return(vec)
      }, medSam)
      return (mat)
    }
    
    obsTable <- normalize.medFC(obsTable)
  } else if (method == "totIon"){
   totSig <- colSums(obsTable)
   # calc factor
   totSig <- totSig/ min(totSig)
   obsTable <- t(apply(obsTable, 1, function(Var){
     Var/ totSig
   }))
  } else {
    stop(method, " is not a recognized normalization method")
  }
  }
  # replace with normalized obsTable
  peakTable[, obsIndx] <- obsTable
  # return zero filled table
  return(peakTable)
}