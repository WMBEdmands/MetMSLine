#' preProc - combined pre-processing steps of the MetMSLine package.
#'
#' @description Wrapper function performs all multiparametric preprocessing steps
#' for large-scale high-resolution LC-MS metabolomic datasets. Combines the
#' \code{\link{zeroFill}},  \code{\link{signNorm}}, \code{\link{loessSmooth}},
#' \code{\link{blankSub}}, \code{\link{cvCalc}} and \code{\link{logTrans}}
#' functions.
#'
#' @details the wrapper function performs the following steps of partially optional
#' data pre-processing steps, not supplying certain function arguments affects which
#' preprocessing steps will be conducted:
#'
#' 1. zero filling (see: \code{\link{zeroFill}} for further details).
#'
#' 2. signal normalization (optional, see: \code{\link{signalNorm}} for further details).
#'
#' 3. blank substraction (see: \code{\link{blankSub}} for further details).
#'
#' 4. pooled QC-based loess smoothing (optional, see: \code{\link{loessSmooth}} for further details).
#'
#' 5. Coefficient of variation calculation and filtration
#' (optional, see: \code{\link{cvCalc}} for further details).
#'
#' 6. log transformation (see: \code{\link{logTrans}} for further details).
#'
#' @param peakTable either a data.frame, full file path as a character string to a  .csv file of a peak table in the form observation (samples) in columns and
#' variables (Mass spectral signals) in rows. If argument is not supplied a GUI file selection window will open and a .csv file can be selected.
#' @param obsNames character vector of observation (i.e. sample/ QC/ Blank) names to identify appropriate observation (sample) columns.
#' @param sampNames character vector of sample names to identify appropriate observation (sample) columns. If
#' either the sampNames or qcNames arguments are not supplied then the pooled QC
#'  \code{\link{loessSmooth}} function is not performed.
#' @param qcNames character vector of quality control (QC) names to identify appropriate observation (quality control) columns. If the qcNames argument is not supplied neither the \code{\link{loessSmooth}}
#' or \code{\link{cvCalc}} functions will be carried out.
#' @param blankNames character vector of blank (i.e. negative control) names to identify
#' appropriate observation (blank) columns. If this and the sampNames argument are supplied the preProc function
#' will conduct a background subtraction (samples:blanks). This is carried
#' out by calculating the fold change (either mean or median: see blankFCmethod) between the samples
#' (see sampNames argument) as the numerator and the blanks as the denominator for each LC-MS feature.
#' Any LC-MS features lower than the fold change threshold (blankFCthresh) will be removed.
#' The fold change threshold can be set as required (default = 2, see blankFCthresh).
#' @param zeroFillValue numeric value to fill zero/ missing values (NA). By default
#' half the mimimum non-zero observed peak intensity is used for \code{\link{zeroFill}} function.
#' @param cvThresh numeric the minimum CV\% to retain an LC-MS variable (default = NULL).
#' If this argument and qcNames are not supplied the CV\% will not be calculated
#' using the quality control samples and no LC-MS features below the CV\%
#' threshold will not be removed. Argument for the \code{\link{cvCalc}} function.
#'  @param outputDir optional directory path to save output images before and after QC
#' smoothing. A subdirectory will be created in which to save the png images. Argument for the \code{\link{loessSmooth}} function.
#' @param smoothSpan numeric (values between 0-1) fixed smoothing span. If supplied a this fixed smoothing
#' span parameter will override the cross validated feature-by-feature smoothing
#' span optimization. Argument for the \code{\link{loessSmooth}} function.
#' @param folds numeric (default=7, i.e. 7-fold cross validation) n-fold cross validation.
#' Argument for the \code{\link{loessSmooth}} function.
#' @param baseLogT the base with respect to which logarithms are computed \code{\link{log}}. Defaults to e=exp(1). Argument for the \code{\link{logTrans}} function.
#' @param normMethod either "medFC" for median fold change or "totIon" for total ion signal normalization.
#' also a custom vector of factors equal in length to the obsNames argument with
#' which to normalize the data can also be supplied. default = "medFC". Argument for the \code{\link{bnormMethod}} function.
#' @param blankFCmethod character either 'median' or 'mean' fold change calculation. Argument for the \code{\link{blankSub}} function.
#' @param blankFCthresh numeric sample:blank fold change cut-off. Any LC-MS
#' features below this threshold will be removed. Argument for the \code{\link{blankSub}} function.
#' @return a data frame identical to peakTable argument (with potentially fewer
#' rows/ LC-MS features if blank substraction or CV\% calculation/ threshold filtration
#' have been performed, see: \code{\link{blankSub}}, \code{\link{cvCalc}}).
#' @examples
#' # read in table of synthetic MS1 metabolomic profiling data
#' peakTable <- system.file("extdata", "synthetic_MS_data.csv", package = "MetMSLine")
#' peakTable <- read.csv(peakTable, header=T, stringsAsFactors=F)
#' # all observation names
#' obsNames <- colnames(peakTable)[grep("QC_|sample_|blank_", colnames(peakTable))]
#' # quality control names
#' qcNames <- colnames(peakTable)[grep("QC_", colnames(peakTable))]
#' # remove all but last column conditioning QC
#' qcNames <- qcNames[-c(9:1)]
#' # sample names only those bounded by qcs
#' sampNames <- colnames(peakTable)[grep("sample_", colnames(peakTable))]
#' sampNames <- sampNames[-c(length(sampNames):{length(sampNames) - 3})]
#' # blank (negative control) names
#' blankNames <- colnames(peakTable)[grep("blank_", colnames(peakTable))]
#' # detect number of cores using parallel package
#' nCores <- parallel::detectCores()
#' # conduct LC-MS data preprocessing
#' preProc_peakTable <- preProc(peakTable, obsNames, sampNames, qcNames, blankNames,
#'                              cvThresh=30, nCores=nCores)
#'
#' @export
preProc <- function(peakTable=NULL, obsNames=NULL, sampNames=NULL, qcNames=NULL,
                    blankNames=NULL, zeroFillValue=NULL, normMethod=NULL,
                    cvThresh=NULL, nCores=NULL, outputDir=NULL, smoothSpan=NULL,
                    folds=7, baseLogT=exp(1), blankFCmethod='mean',
                    blankFCthresh=2){
  #error handling
  if(is.null(obsNames)){
    stop('argument obsNames is missing with no default')
  }
  # error handling or read from csv function
  peakTable <- tableCheckRead(peakTable, stringsAsFactors=FALSE)
  # match obsNames to peak table colnames
  obsIndx <- match(obsNames, colnames(peakTable))
  # if less than all matched then stop
  if(length(obsIndx) < length(obsNames)){
    stop(length(obsIndx), " of ", length(obsNames),
         " observation names were matched in the peakTable column names, check the obsNames and peakTable column names")
  }

  # 1. zero fill
  peakTable <- zeroFill(peakTable, obsNames, value=zeroFillValue)
  # 2. if necessary normalize
  if(!is.null(normMethod)){
  peakTable <- signNorm(peakTable, sampNames, method=normMethod)
  } else {
  message('no normalization will be conducted by normMethod function...\n')
  flush.console()
  }
  # 3. if necessary peform blank subtraction
  if(!is.null(sampNames) & !is.null(blankNames)){
    peakTable <- blankSub(peakTable, blankNames, sampNames, blankFCmethod,
                          blankFCthresh)
  } else {
    message('no blank (negative control) based background subtraction will be performed...\n')
    flush.console()
    message('two vectors of sample names (sampNames) and blanks/ negative controls (blankNames) must be supplied...\n')
    flush.console()
  }
  # 4. if necessary peform loess smoothing
  if(!is.null(sampNames) & !is.null(qcNames)){
  peakTable <- loessSmooth(peakTable, sampNames, qcNames, nCores=nCores,
                           outputDir=outputDir, folds=folds,
                           smoothSpan=smoothSpan)
  # remove any negative variables
  if(!is.null(peakTable$negVals)){
    message("Automatically removing ", sum(peakTable$negVals),
            " LC-MS features containing negative values following loess smoothing...:\n\n",
            paste0('EIC ', peakTable[peakTable$negVals, 1][1:ifelse(sum(peakTable$negVals) > 10, 10, sum(peakTable$negVals))], collapse = '\n'),
            '...(limited to first 10 peak groups)\n')
    flush.console()
    peakTable <- peakTable[peakTable$negVals == FALSE, ]
    # remove negVals column
    peakTable$negVals <- NULL
   }

  } else {
  message('no QC based signal attenuation correction will be performed...\n')
  flush.console()
  message('two vectors of sample names (sampNames) and quality controls (qcNames) must be supplied...\n')
  flush.console()
  }

  # 5. Calc. coeff. of Var. if necessary
  if(!is.null(cvThresh) & !is.null(qcNames)){
    peakTable <- cvCalc(peakTable, qcNames, cvThresh)
  } else {
    message('No coefficient of variation (CV%) based LC-MS subtraction will be performed...\n')
    flush.console()
    message('both The cvThresh argument and a vector of quality control names (qcNames) must be supplied...\n')
    flush.console()
  }
  # 6. Log transform
  peakTable <- logTrans(peakTable, obsNames, base=baseLogT)
  # return final peak table
  return(peakTable)
} # end function
