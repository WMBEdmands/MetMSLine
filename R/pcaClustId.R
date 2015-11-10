#' Identify clusters in a PCA plot from a list of co-variates using PAM clustering.
#'  
#' @description attempts to the identify clusters in a pca from a data frame
#' of covariates, using partitioning around the medoid clustering.
#' 
#' @param pcaResult a \code{\link{pcaRes}} class object
#' @param peakTable optional if pcaResult argument not supplied. Either a data.frame, full file path as a character string to a  .csv file of a peak table in the form observation (samples) in columns and
#' variables (Mass spectral signals) in rows. If argument is not supplied a GUI file selection window will open and a .csv file can be selected.
#' @param coVarTable either a data.frame, full file path as a character string to a  .csv file of a co-variates table in the form observation (sample) names in the first column and
#' co-variates from the 2nd column onward. If argument is not supplied a GUI file selection window will open and a .csv file can be selected. 
#' @param obsNames character vector of observation (i.e. sample/ QC/ Blank) names to identify appropriate observation (sample) columns.
#' @param ... additional arguments to \code{\link{pamk}}.
#' 
#' @details potential clusters in \code{\link{pcaRes}} scores are identified using partitioning
#' around the medoid clustering (\code{\link{pamk}}) from the fpc package with an
#' estimation of the number of clusters. Given a data.frame of co-variates/ sample information, 
#' the most likely explanatory co-variate AND/OR potential two-factor interactions will be established using linear modelling 
#' \code{\link{lm}}.
#' 
#' Co-variates containing missing values or only one unique value will not be considered.
#' The linear model consists of response ~ terms where response is the clusters
#' established by \code{\link{pamk}} and terms are the factor levels of the
#' co-variate table. The best explanatory co-variate for the PCA clustering 
#' is defined as the linear model with the highest coefficient of determination (R2).
#' 
#' @return a list containing three elements:
#' 
#' 1. a list pamkClust containing three elements returned from \code{\link{pamk}}
#' 
#' 2. a list lmCoVarClust containing the linear models obtained from \code{\link{lm}}.
#' 
#' 3. a named character vector rSquaredLm containing the coefficients of determination
#' from the linear models. Named with each covariate or two-factor interaction
#' considered.
#' 
#' @seealso \code{\link{lowess}}, \code{\link{na.spline}}, \code{\link{pcaOutId}}.
#' @export
pcaClustId <- function(pcaResult=NULL, peakTable=NULL, coVarTable=NULL, 
                       obsNames=NULL, ...){
  # if a pcaRes object is not supplied then peakTable
  if(is.null(pcaResult)){
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
  
  # calc first Pca model
  message("Calculating PCA model ", nIter, "...")
  flush.console()
  # calc Pca Model
  pcaResult <- pcaMethods::pca(t(obsTable), nPcs=2, ...)
  
  } else if(class(pcaResult) != "pcaRes"){
  # check pcaRes class object 
  stop('pcaResult arguments is not a "pcaRes" class object ?pcaRes')
  }
  # error handling or read from csv function
  coVarTable <- tableCheckRead(coVarTable, type="co-variates")
  # check scores rownames match coVariates table
  scoresIndx <- match(coVarTable[, 1], rownames(pcaResult@scores))
  # error handling
 if(all(is.na(scoresIndx))){
    stop("No co-variate table sample names in the first column match sample names in the PCA scores")
  }
  # id covartable index
  coVarIndx <- which(is.na(scoresIndx) == F)
  # remove any NAs from scoresIndx
  scoresIndx <- scoresIndx[complete.cases(scoresIndx)]
 
  if(length(scoresIndx) < nrow(coVarTable)){
    warning("Only ", length(scoresIndx), 
         " sample names in the first column of the coVarTable match sample names in the PCA scores")
  }
  # PAMK clustering using scores
  pamkClust <- try(fpc::pamk(pcaResult@scores[scoresIndx, ]), silent=T)
  # catch if there is an error with pamk (only one cluster identified)
  if(class(pamkClust) == "try-error"){
    message("No PCA score plot clusters were identified...")
    flush.console()
    return(list(pamkClust=NULL, lmCoVarClust=NULL, rSquaredLm=NULL))
  } else {
  nClust <- length(unique(pamkClust$pamobject$clustering))
  message(nClust, " PCA scores clusters identified by PAM")
  # check for missing values or only one factor level in co-variates table
  missingIndx <- apply(coVarTable, 2, function(coVar){ 
    any(is.na(coVar)) | any(coVar == "") | length(unique(coVar)) == 1})
  if(any(missingIndx)){
  cat("The following coVarTable columns contain missing values or only one value and will not be considered:\n",
    paste0(seq(1, sum(missingIndx), 1), ". ", colnames(coVarTable)[missingIndx], "\n"))
  if(all(missingIndx)){
    stop("All coVarTable columns contain missing values")
  }
  }
  # subset coVartable with missing Index co-variates table
  coVarTable <- coVarTable[coVarIndx, missingIndx == F, drop=F]
  # all two factor interactions
  factors <- colnames(coVarTable)[2:ncol(coVarTable)]
  factors <- paste0('coVarTable[, "', factors, '"]')
  factorsM <- outer(factors, factors, paste, sep=" + ")
  factorsM[upper.tri(factorsM, diag = T)] <- ""
  factorsM <- as.vector(factorsM)
  factors <- c(factors, factorsM[factorsM != ""])
  # convert all character to numeric factor levels
  coVarTable <- apply(coVarTable, 2 , function(x) as.numeric(as.factor(x)))
  # convert pamK to numeric
  pamClusts <- as.numeric(pamkClust$pamobject$clustering)
  # lm clustering covariates
  lmRes <- lapply(factors, function(coVar){ 
                  form <- as.formula(paste0("pamClusts ~ ", coVar))
                  return(lm(form))})
  # extract R2 and rank covars
  rSquaredLm <- format(sapply(lmRes, function(coVar) round(summary(coVar)$r.squared, digits=4)), 
                       scientific=F)
  names(rSquaredLm) <- gsub("coVarTable\\[, |\\]", "", factors)
  cat("largest coefficient of determination (r2) from linear modelling of all\n", 
      " co-variates and their two-factor interactions to PCA scores clusters: \n") 
  print(rSquaredLm[which.max(rSquaredLm)])
 
  # return pca Results and outliers removed peakTable
  return(list(pamkClust=pamkClust, lmCoVarClust=lmRes, rSquaredLm=rSquaredLm))
  }
}