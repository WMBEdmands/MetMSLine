#' loess with 7 - fold CV function
#' @param x predictor values
#' @param y response values
#' @param original LC-MS feature.
#' @param folds numeric (default=7, i.e. 7-fold cross validation).
loessFitKfoldCV <- function(x, y, feature, folds=7){
  # calc. loess model fit using 7 fold cross validation
tmpLoess <- suppressWarnings(loessWrapperMod(x, y, folds = folds))
# extract optimum span value
spanTmp <- tmpLoess$pars$span
# predict all qc and sample intensities using predict loess
tmpPred <- predict(tmpLoess, newdata=1:length(feature))
# multiply the orignal variable by the deviation of the loess predicted value 
# from its median
cbind(feature * (median(tmpPred)/ tmpPred), spanTmp)
}