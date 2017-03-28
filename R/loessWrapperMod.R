#' adapted from \link{\code{bisoreg}} package
#' 
loessWrapperMod<- function (x, y, data, span.vals = seq(0.25, 1, by = 0.05), folds = 5){
  # sort from highest to lowest span values
  orderTmp <- order(span.vals, decreasing = TRUE)
  span.vals <- span.vals[orderTmp]
  mae <- numeric(length(span.vals))
  theta.fit <- function(x, y, span) loess(y ~ x, span = span, surface='direct')
  theta.predict <- function(fit, x0) predict(fit, newdata = x0)
  ii = 0
  for (span in span.vals){
    ii <- ii + 1
    y.cv <- bootstrap::crossval(x, y, theta.fit, theta.predict, span = span, 
                                ngroup = folds)$cv.fit
    fltr <- !is.na(y.cv)
    mae[ii] <- mean(abs(y[fltr] - y.cv[fltr]))
  }
  # remove duplicate values
  dupIdx <- duplicated(mae) == FALSE
  span.vals <- span.vals[dupIdx]
  mae <- mae[dupIdx]
  span <- span.vals[which.min(mae)]
  out <- loess(y ~ x, span = span, surface='direct')
  while(is.na(out$s) | out$s == Inf){
  # remove min
  span.vals <- span.vals[-which.min(mae)]
  mae <- mae[-which.min(mae)]
  if(length(span.vals) == 0){
    break
  }
  span <- span.vals[which.min(mae)]
  out <- loess(y ~ x, span = span, surface='direct')
  }
  return(out)
}