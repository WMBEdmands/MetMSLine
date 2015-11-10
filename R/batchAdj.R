#' linear model adjustment for batch
#' 
#' @description linear modelling of a peak table adjusting for batch \code{\link{lm}}.
#' i.e. peakTable ~ batchEff. Residuals are used as normalization factors for
#' each variable variable.
#' 
#' @param peakTable either a data.frame, full file path as a character string to a  .csv file of a peak table in the form observation (samples) in columns and
#' variables (Mass spectral signals) in rows. If argument is not supplied a GUI file selection window will open and a .csv file can be selected.
#' @param obsNames character vector (i.e. sample/ QC/ Blank) names to identify appropriate observation (sample) columns.
#' @param batchEff a vector or matrix/ data.frame of co-variates implicated in batch membership.
#' If a matrix/data.frame is supplied the linear model will be constructed with two or more fixed effects.
#' (i.e. in the form peakTable ~ batchEff1 + batchEff2). Categorical variables (e.g. batch 1 and batch 2) must
#' be encoded as factors and continuous variables encoded as numeric for the batch adjustment
#' to work properly.
#' @param ... additional variables to \code{\link{lm}}. will be added to 
#' @return a list containing 
#' 1. a data frame identical to peakTable with variables adjusted for batch effect.
#' 2. a matrix of the linear model residuals.
#' @export
batchAdj <- function(peakTable=NULL, obsNames=NULL, batchEff=NULL){
  
  if(is.null(obsNames)){
    stop('argument obsNames is missing with no default')
  }
  if(is.null(batchEff)){
    stop('argument batchEff is missing with no default.\n At least one batch effect vector or a matrix/data.frame containing multiple interacting batch effects must be supplied')
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
  message("Adjusting for batch effects...")
  flush.console()
  # calc linear model
  obsTable <- apply(obsTable, 2, as.numeric)
  if(is.null(dim(batchEff))){
  batch.lm <- lm(t(obsTable) ~ batchEff)
  } else {
  tObsTable <- as.matrix(t(obsTable))
  form <- as.formula(paste0("tObsTable ~ ", paste0("batchEff$", colnames(batchEff), collapse=" + ")))
  batch.lm <- lm(form)  
  }
  # extract residuals and scale each variable to minimimum value
  batchRes <- batch.lm$residuals
#   batchRes.scaled <- t(apply(batchRes, 2, function(x){1 + (x + abs(min(x)))}))
  # batch adjusted
  obsTable.adj <- obsTable/ t(batchRes) #.scaled
  # replace peak table values with residuals
  peakTable[, obsIndx] <- t(batchRes)
  batchAdj <- peakTable
  peakTable[, obsIndx] <- obsTable.adj
  # return zero filled table
  return(list(peakTable=peakTable, batchRes=batchAdj))
}