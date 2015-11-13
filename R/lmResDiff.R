#' linear model residual differences by subtraction and continuous var (i.e.
#'days to diagnosis)
#' @param batchAdjusted a matrix or data.frame output of batch adjusted residuals obtained
#' from linear modelling of analytical batch covariate \code{\link{batchAdj}}. 
#' The matrix/ data.frame must take the form observations (samples) in columns and
#' mass spectral variable residuals in rows.
#' @param obsNames character vector of observation (i.e. sample/ QC/ Blank) names to identify appropriate observation (sample) columns.
#' @param resDiffId a character or logical vector containing the pair/ class identity of observations.
#' This identity will be used to substract the batchAdjusted residuals. (e.g. case-control status).
#' @param contVar a character/ factor/ numeric vector containing a continuous variable which will be modelled
#' \code{\link{lm}}
#' against the residual difference values. This vector will be subset and only the TRUE, maximum factor level
#' from resDiffId will be included.
#' @param minPintercept the minimum p value threshold for the intercept value of \code{\link{lm}}. 
#' @param minPintercept the minimum p value threshold for the intercept value of \code{\link{lm}}. 
#' @param minPcontVar he minimum p value threshold for the continuous variable contVar of \code{\link{lm}}
#' @param MTC multiple testing correction method see(\code{\link{p.adjust}}). default = "BH".
#' @return returns a list containing two elements:
#' 1. resultsTable = a data.frame containing the variables found to be below both
#' p value thresholds (minPintercept and minPcontVar).
#' 2. resDiff = a data.frame containing the substraction of the resDiffId groups.
#' @export
lmResDiff <- function(batchAdjusted=NULL, obsNames=NULL, resDiffId=NULL, 
                      contVar=NULL, minPintercept=0.005, minPcontVar=0.005, 
                      MTC="BH"){
  # error handling
  if(is.null(obsNames)){
    stop('argument obsNames is missing with no default')
  } else if(is.null(batchAdjusted)){
    stop("the argument batchAdjusted is missing with no default")
  } else if(is.null(resDiffId)){
    stop("the argument resDiffId is missing with no default")
  } else if(is.null(contVar)){
    stop("the argument contVar is missing with no default")
  } else if(length(resDiffId) != length(contVar)){
    stop("the vectors resDiffId and contVar are of differing lengths, please check and try again")
  } else if(length(resDiffId) != length(obsNames)){
    stop("the of the observation (sample) names (obsNames) do not match the length of the resDiffId argument, please check and try again")
  }
  
  # error handling or read from csv function
  batchAdjusted <- tableCheckRead(batchAdjusted, type="batch adjusted peak table", 
                                  stringsAsFactors=F)
  # match obsNames to peak table colnames
  obsIndx <- match(obsNames, colnames(batchAdjusted))
  # if less than all matched then stop
  if(length(obsIndx) < length(obsNames)){
    stop(length(obsIndx), " of ", length(obsNames), 
         " observation names were matched in the peakTable column names, check the obsNames and peakTable column names")
  }
  # subset table
  obsTable <- batchAdjusted[, obsIndx]
  # check if any zeros or NAs 
  if(any(is.na(obsTable))){
    stop("peakTable observations contain NAs use the ?zeroFill function")
  }
  
  # transpose obsTable
  obsTable <- t(obsTable) 
  # subtract residuals 1 and 2
  dfRes1 <- obsTable[resDiffId == T, ]
  dfRes2 <- obsTable[resDiffId == F, ]
  resDifferences <- dfRes1 - dfRes2
  # linear model res differences and continuous var
  message("calculating linear model and extracting coefficients...")
  flush.console()
  residContVarLm <- lm(resDifferences ~ contVar[resDiffId == T])
  
  # extract coefficients and R squared
  residContVarDf <- as.data.frame(do.call("rbind", summary(residContVarLm)))
  residContVarDf$R <- sqrt(as.numeric(residContVarDf$r.squared))
  residContVarDf <- as.data.frame(residContVarDf[, c("coefficients","r.squared")])
  residContVarDf <- as.data.frame(t(apply(residContVarDf, 1, unlist)))
  residContVarDf <- as.data.frame(residContVarDf[,c("coefficients1", "coefficients7", 
                                                    "coefficients2", "coefficients8", 
                                                    "r.squared")])
  colnames(residContVarDf) <- c("Intercept Coefficient", "Intercept p Value", 
                                "ContVar Coefficient", "ContVar p Value", "rSquared")
  # control for false discovery
  residContVarDf$AdjPValueContVar <- p.adjust(residContVarDf[,"ContVar p Value"], MTC)
  # merge analysis results and original batch adjusted data frame/ matrix
  residContVarDf <- cbind(residContVarDf, batchAdjusted)
  
  # restrict to lowest p-values intercept and TTD
  SFeatInterContVar <- residContVarDf[residContVarDf[, "Intercept p Value"] < minPintercept & 
                                        residContVarDf[, "ContVar p Value"] < minPcontVar, , drop=F]
  return(list(resultsTable=SFeatInterContVar, resDiff=resDifferences))
}