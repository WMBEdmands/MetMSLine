#' loess with n - fold CV function
#' @param x predictor values
#' @param y response values
#' @param feature LC-MS feature.
#' @param smoothSpan fixed smoothing span. If supplied a this fixed smoothing
#' span parameter will override the cross validated feature-by-feature smoothing
#' span optimization.
#' @param folds numeric (default=7, i.e. 7-fold cross validation).
#' @param plotFile character full path to plot file png (optional). If supplied
#' a plot file will be generated at this location showing the before and after
#' loess smoothing result.
#' @param colours vector of colours for datapoints (i.e. red for QCs and black for
#' samples). (optional if plot file supplied then this argument is required). 
#' Argument is passed to \link{plot} col. 
#' @return vector smoothed feature and smoothing span parameter selected concatenated
#' to end.
loessFitKfoldCV <- function(x, y, feature, smoothSpan=NULL, folds=7, 
                            plotFile=NULL, colours=NULL){

# calc. loess model fit using n-fold cross validation
if(is.null(smoothSpan)){
tmpLoess <- suppressWarnings(loessWrapperMod(x, y, folds = folds))
} else {
  if(!is.numeric(smoothSpan)){
    stop('smoothSpan argument must be a numeric value between 0-1')
  }
  if(smoothSpan > 1 | smoothSpan < 0){
    stop('smoothSpan argument must be a numeric value between 0-1')
  }
tmpLoess <- suppressWarnings(loess(y ~ x, span = smoothSpan, surface='direct'))
}
# check if worked
# extract optimum span value
spanTmp <- tmpLoess$pars$span
# interpolate values of curve
curveTmp <- vector('numeric', length(feature))
# all na
curveTmp[] <- NA
# add fitted values
curveTmp[tmpLoess$x] <- tmpLoess$fitted
# cubic spline interpolation
curveTmp <- zoo::na.spline(curveTmp)
# multiply the orignal variable by the deviation of the loess predicted value 
# from its median
adjFeat <- feature * {median(curveTmp)/ curveTmp}
# if necc generate plot
if(!is.null(plotFile)){
  if(grepl('\\.png$', plotFile) == FALSE){
    stop('The plot file name must end in .png')
  }
  rangeTmp <- range(range(feature), range(adjFeat))
  png(plotFile, width=2400, height=1200, res=150)
  par(mfrow=c(1, 2))
  plot(1:length(feature), feature, main='raw', pch=21, col='black', bg=colours,
       xlab='sequence', ylab='peak area/height', ylim=rangeTmp, cex=1.5)
  # fit line
  points(1:length(feature), curveTmp, type='p', col='black', bg='blue', pch=24, 
         cex=0.6)
  plot(1:length(feature), adjFeat, 
       main=paste0('corr (smoothSpan = ', spanTmp, ')'),
       pch=21, col='black', bg=colours, xlab='sequence', 
       ylab='peak area/height', ylim=rangeTmp, cex=1.5)
  dev.off()  
}
return(cbind(adjFeat, spanTmp))
}