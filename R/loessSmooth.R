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
#' @param nCores number of cores to use for parallel processing with the foreach
#' and snow packages.  
#' @param outputDir optional directory path to save output images before and after QC
#' smoothing. A subdirectory will be created in which to save the png images.
#' @param smoothSpan fixed smoothing span. If supplied a this fixed smoothing
#' span parameter will override the cross validated feature-by-feature smoothing
#' span optimization.
#' @param folds numeric (default=7, i.e. 7-fold cross validation).
#' @details systematic drifts in mass spectral signal can be corrected using
#' intervally injected pooled quality control samples for each mass spectral signal
#' variable in the peak table. The optimum span parameter for each LC-MS feature 
#' is identified using 7-fold cross-validation 
#' (see: \code{\link{crossval}} from the bootstrap package). A new column is
#' added to the peakTable "smoothSpanLoessFit" containing the result of the 
#' optimal value of span for the loess function. 
#'  For each LC-MS variable \code{\link{loess}} models
#' of the quality control data are calculated using direct surface fitting (see: \code{\link{loess.control}}) and
#' the loess fit using 7-fold cross validation is carried out by a slight
#' modification of the \code{\link{loess.wrapper}} function from the \code{\link{bisoreg}}
#' package. 
#' The optimum loess model is then used to predict (see: \code{\link{predict.loess}}) 
#' intensity values for both the QC and sample injections according
#' to the degree of smoothing (span) in the final optimum loess model. 
#' The product of the original LC-MS variable intensity and the deviation of the 
#' loess predicted values from their median is then calculated to adjust the 
#' signal drift. 
#' This function attempts to
#' adjust systematic drift/ attenuation within each LC-MS feature. 
#' The degree of smoothing determined by the smoothSpan argument however results
#' may vary slightly using the K-fold cross validation method. In order to maintain
#' absolute reproducibility it is recommend to set the seed prior to using this 
#' function (see: \code{\link{set.seed}}). It is recommended that parallel processing
#' is considered to reduce the computational time cost using the argument nCores.
#'  
#' @return a data frame identical to peakTable with signal drift/ attenuation adjusted.
#' @source \url{http://www.nature.com/nprot/journal/v6/n7/abs/nprot.2011.335.html} 
#' Procedures for large-scale metabolic profiling of serum and plasma using gas 
#' chromatography and liquid chromatography coupled to mass spectrometry
#' @seealso \code{\link{loess}}, \code{\link{loess.wrapper}}, \code{\link{bisoreg}},
#'  \code{\link{crossval}}, \code{\link{predict.loess}}, \code{\link{set.seed}}.
#' @export
loessSmooth <- function(peakTable=NULL, sampNames=NULL, qcNames=NULL, 
                        nCores=NULL, outputDir=NULL, smoothSpan=NULL, folds=7){
  # error handling
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
  
  # check if all samples bounded by qcs
  if(max(qcIndxObs) < max(sampIndxObs) | min(qcIndxObs) > min(sampIndxObs)){
    stop('All samples must be bounded by quality control samples.')
  }
  # convert to numeric in case character or factor
  obsTable <- data.frame(apply(obsTable, 2, as.numeric), stringsAsFactors=F)
  # new Df observatin and QC indx for loess modelling
  qcsTable <- data.frame(t(obsTable[, qcIndxObs]), qcIndxObs, stringsAsFactors=F)
  
  if(folds > length(qcIndxObs)){
    message('number of QCs less than number of folds (cross-validation) reducing to ', length(qcIndxObs))
    flush.console()
    folds <- length(qcIndxObs)
  }
  
  message("Calculating LOESS fit (", folds, "-fold CV) using ", length(qcIndxObs), " QC samples and signal drift/attenuation smoothing...\n")
  flush.console()
  
  qcDummyIdx <- {1:ncol(obsTable) %in% qcIndxObs} + 1
  # show plots before and after smoothing
  if(!is.null(outputDir)){
    lSmDir <- paste0(outputDir, '/loessSmooth/')
    dir.create(lSmDir)
    plotFiles <- paste0(lSmDir, '/', 'EIC_', peakTable[, 1], '.png')
    message('Generating ', nrow(peakTable), ' output images in output directory:\n',
            lSmDir, '\n\n')
    flush.console()
  }
  # if nCores !null then start snow cluster
  if(!is.null(nCores)){
    if(!require(foreach)){
      stop('The foreach package (install.packages("foreach")) must be installed to perform parallel computation...')
    }
    message(paste0("Starting SNOW cluster with ", nCores, 
                   " local sockets...\n"))
    flush.console()
    cl <- parallel::makeCluster(nCores)
    doSNOW::registerDoSNOW(cl)
    message(folds, "-fold cross validaton loess fitting will be applied to ", 
            nrow(peakTable), " LC-MS features. Please wait...\n")
    flush.console()
   
    cat(paste0('Progress (', nrow(peakTable), ' LC-MS features):\n'))
    progSeq <- round({nrow(peakTable) * seq(0, 1, 0.05)}, 0)
    progSeq[1] <- 1
    progress <- function(n){if(n %in% progSeq){cat(paste0(round({n/nrow(peakTable)} * 100, 0), '%  '))}}
    opts <- list(progress=progress)
    
    
    obsTable <- foreach(Var=1:nrow(obsTable), .packages=c("bootstrap", "zoo"), 
                         .combine='rbind', .options.snow=opts) %dopar% {
      loessFitKfoldCV(qcIndxObs, qcsTable[, Var], 
                      obsTable[Var, ], smoothSpan=smoothSpan, folds = folds, 
                      plotFile=switch({!is.null(outputDir)} + 1, NULL, plotFiles[Var]),
                      colours=c("black", "red")[qcDummyIdx])
    }
    parallel::stopCluster(cl)
    # rename last smoothSpanLoessFit column
    colnames(obsTable)[ncol(obsTable)] <- 'smoothSpanLoessFit'
    # or single thread  
  } else {
  # new column for saving optimal span value by 7-fold cross validation
  obsTable$smoothSpanLoessFit <- 0
  # est. progress bar
  pb <- txtProgressBar(min=0, max=nrow(obsTable), style=3)
  # calc. loess models each LC-MS variable return predict LC-MS intensities
  for(Var in 1:nrow(obsTable)){
    # set progress bar
    setTxtProgressBar(pb, Var)
    obsTable[Var, ] <- loessFitKfoldCV(qcIndxObs, qcsTable[, Var], 
                                        obsTable[Var, -ncol(obsTable)],
                                        smoothSpan=smoothSpan,
                                        folds = folds, 
                                        plotFile=switch({!is.null(outputDir)} + 1, NULL, plotFiles[Var]),
                                        colours=c("black", "red")[qcDummyIdx])
  }
  }
  
  # warn user if any variables contain negative values
   negVals <- apply(obsTable[, 1:(ncol(obsTable) - 1)], 1, function(Var) any(Var < 0))
   lNegValsIndx <- length(which(negVals))
   if(lNegValsIndx > 0){
     message(paste0(lNegValsIndx, ' (', (lNegValsIndx/ length(negVals) * 100), '%) of the LC-MS variables contain one of more negative values following loess smoothing...\n\n', 
                                 'A column of logicals "negVals" will be added to the returned table indicating these...\n\n',
                                 'Please examine and remove these potentially problematic features before proceeding with further analysis...\n\n'))
     # create new column in peak table
       peakTable$negVals <- negVals  
       }
  # replace with adjusted values obsTable
  peakTable[, sort(c(sampIndx, qcIndx))] <- obsTable[, 1:(ncol(obsTable) - 1)]
  # add optimal span parameter column
  peakTable$smoothSpanLoessFit <- obsTable[, ncol(obsTable)]
  # return zero filled table
  return(peakTable)
}