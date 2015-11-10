#' rtCorrClust
#'
#' @description performs retention time clustering thin intra-retention time cluster
#' correlation coefficient (spearman's Rho) clustering. 
#' 
#' @param peakTable either a data.frame, full file path as a character string to a  .csv file of a peak table in the form observation (samples) in columns and
#' variables (Mass spectral signals) in rows. If argument is not supplied a GUI file selection window will open and a .csv file can be selected.
#' @param obsNames character vector of observation (i.e. sample/ QC/ Blank) names to identify appropriate observation (sample) columns.
#' @param rtThresh retention time threshold for retention time clustering.
#' @param corrThresh correlation coefficient threshold (non-parametric Spearman's 
#' Rho) to group features within a retention time cluster.
#' @param minFeat minimum number of features with a Retention time/ correlation
#' cluster to consider it a group (default = 2, i.e. a cluster must contain at
#' least 2 features to be considered a group).
#' @param hclustMethod hierarchical clustering method to \code{\link{hclust.vector}} 
#' method of fastcluster package (default = "median").
#' @param distMeas distance measure for retention time clustering (default = "euclidean"). 
#' see \code{\link{hclust.vector}}.
#' @param corMethod correlation method to \code{\link{cor}} (default = "spearman).
#' Within retention clusters dissimilarity is computed as 1-correlation coefficient.
#' 
#' @details The \code{\link{cutree}} function is used to identify retention time clusters based on a cut height
#' equal to retention time (rtThresh) supplied.
#' Within retention time cluster correlation clustering is based on 1 - correlation
#' coefficient dissimilarity. The cut height of intra retention time cluster
#' correlation clusters is based on a minimum correlation coefficient threshold
#' (default = 0.6). Following retention time/ correlation cluster identification
#' any groups containing less than a mimimum feature number are removed 
#' (minFeat default = 2), a default value for the minFeat of 2 means that only
#' isolated single features with no relation to another feature are removed.
#' @return a list containing 2 elements: 
#' 
#'  1. the original peak table with the results of retention time 
#'  correlation clustering included in additional columns "rtgroup", "CorrClusts", 
#' "RtCorrClust" as a data frame and details of the most intense 
#'  feature for each group. In addition for each group details of the best sample 
#'  to consider for subsequent MSn fragmentation analysis based on the most 
#'  top 5 most intense samples is provided. This "best" sample for MS/MS is also
#'  selected whilst concurrently minimising the number of samples to potentially 
#'  reinject for targetted MSn fragmentation analysis.
#'  
#'  2. A data frame of the weighted means and details of the most intense 
#'  feature for each group. In addition for each group details of the best sample 
#'  to consider for subsequent MSn fragmentation analysis based on the most 
#'  top 5 most intense samples is provided. This "best" sample for MS/MS is also
#'  selected whilst concurrently minimising the number of samples to potentially 
#'  reinject for targetted MSn fragmentation analysis. 
#'   @seealso \code{\link{hclust.vector}}, \code{\link{cutree}}.
#' @export
rtCorrClust <- function(peakTable=NULL, obsNames=NULL, rtThresh=NULL, corrThresh=0.6, 
                        minFeat=2, hclustMethod="median", distMeas="euclidean",
                        corMethod="spearman"){
  
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
  # replace any NAs with zero
  obsTable <- peakTable[, obsIndx]
  
  #rownames(obsTable) <- peakTable[, 1]
  message("hierarchical clustering peak group retention times... \n")
  flush.console()
  
  hr <- fastcluster::hclust.vector(peakTable[, "rtmed"], metric=distMeas, 
                                   method=hclustMethod)
 if(is.null(rtThresh)){
   # if contains rtmin and rtmax columns
   if(all(c("rtmin", "rtmax") %in% colnames(peakTable))){
     rtThresh <- median((peakTable[, "rtmax"] - peakTable[, "rtmin"])/2)
   }
  }
  # cut tree according to rtThresh
  rtgroups <- cutree(hr, h=rtThresh)
   
  message("intra RT group correlation clustering ", length(unique(rtgroups)),
          " rt groups... \n")
  flush.console()
  
  RtCorrClusts <- unlist(by(obsTable, rtgroups, function(x){
    if(nrow(x) >= 2){
    cor.m <- suppressWarnings(cor(t(x), method=corMethod))
    cor.m[is.na(cor.m)] <- 0
    hr.tmp <- hclust(as.dist(1 - abs(cor.m)), 
                     method="average", members=NULL)
    corGroups.tmp <- cutree(hr.tmp, h=1 - corrThresh)
    freq.tmp <- table(corGroups.tmp)
    freq.tmp <- as.numeric(names(freq.tmp)[which(freq.tmp < minFeat)])
    if(length(freq.tmp) > 0){
      corGroups.tmp[corGroups.tmp %in% freq.tmp] <- 0
    }
    } else {
      corGroups.tmp <-as.numeric(rep(0, nrow(x)))
      names(corGroups.tmp) <- rownames(x)
    }
   return(corGroups.tmp)
  }))
  # rt corr clust array from by to data.frame
  RtCorrClusts <- data.frame(Feat=gsub("^[0-9]*\\.", "", names(RtCorrClusts)), 
                             rtgroup=gsub("\\.X.+|\\.[0-9]*$", "", names(RtCorrClusts)),
                             CorrClusts=RtCorrClusts, stringsAsFactors=F)
  # match row order
  RtCorrClusts <- RtCorrClusts[match(rownames(obsTable), RtCorrClusts$Feat), ]
  # bind with group Val
  RtCorrClusts <- data.frame(peakTable, RtCorrClusts, stringsAsFactors=F)
  # create index and unque Rt clust group names
  clustIndx <- which(RtCorrClusts$CorrClusts != 0)
  RtCorrClusts$RtCorrClust <- ""
  RtCorrClusts$RtCorrClust[clustIndx] <- paste0(RtCorrClusts$rtgroup[clustIndx],
                                                "_", RtCorrClusts$CorrClusts[clustIndx])
  # calculate top 5 intensities
   top5int <- as.data.frame(t(apply(RtCorrClusts[, obsIndx], 1, function(x){
   # calculate average intensity for max weight
   meanInt <- mean(as.numeric(x))
   # identify top 5 max Int samples 
   top5int <- colnames(RtCorrClusts)[obsIndx][head(order(x, decreasing=T), n=5)]
   return(c(meanInt, top5int))})))
   
   colnames(top5int) <- c("meanInt", paste0("mostIntSamp_", 1:5))

  # calculate best sample for MS/MS based on frequency of top 5
   # frequency of most intense samples
   freqMostInt <- sort(table(unlist(top5int[, 2:6])), decreasing=T)
   # add a column for the best sample to select for MS/MS
   top5int$bestSampMSMS <- apply(top5int[, 2:6], 1, 
                                 function(x){
                                 indx.tmp <- which(names(freqMostInt) %in% x)
                                 bestSamp.tmp <- names(freqMostInt)[indx.tmp[1]]
                                 return(bestSamp.tmp)})
  # add top5int data frame into RtCorrclusts
  RtCorrClusts <- cbind(RtCorrClusts[, setdiff(1:ncol(RtCorrClusts), obsIndx)],
                        top5int, RtCorrClusts[, obsIndx])
  obsIndx <- match(obsNames, colnames(RtCorrClusts)) 
  # identify rt correlation clusters
  # calculate weighted mean clusters
  message("Calculating weighted mean for ", 
          length(unique(RtCorrClusts$RtCorrClust[clustIndx])), 
          " pseudospectra accounting for ", length(clustIndx), " of ", 
          nrow(RtCorrClusts), " total features ")
  flush.console()
  # calculate weighted.means pseudospectra and return max intensity feature
  # and sample
  wMeanPspec <- by(RtCorrClusts[clustIndx, ], 
                   as.factor(RtCorrClusts$RtCorrClust[clustIndx]), function(x){
                     # calc weights each variable
                     wts.tmp <- rowSums(x[, obsIndx])
                     wt.mean.tmp <- apply(x[, obsIndx], 2, function(y){
                                          weighted.mean(y, wts.tmp)})
                     # concatenate max feature intensity and weighted means
                     maxFeat.indx <- which.max(wts.tmp)
                     # req columns from xcms diff
                     reqCols <- seq(1, ncol(RtCorrClusts), 1)
                     reqCols <- setdiff(reqCols, obsIndx)
                     # sum features in pseudospectra
                     nFeatPspec <- nrow(x)
                     names(nFeatPspec) <- "nFeatPspec"
                     wt.mean.tmp <- c(x[maxFeat.indx, reqCols], nFeatPspec, 
                                      wt.mean.tmp)
                     return(wt.mean.tmp)})
  # convert by array to data frame
  wMeanPspec <- data.frame(do.call(rbind, wMeanPspec), stringsAsFactors=F)
  # coerce list to character
  wMeanPspec <- data.frame(lapply(wMeanPspec, as.character), stringsAsFactors=F)
  # return a list containing the original xcms diff report with rt corr clust
  # info and the weighted mean data frame
  return(list(diffRepAdd=RtCorrClusts, wMeanPspec=wMeanPspec))
}