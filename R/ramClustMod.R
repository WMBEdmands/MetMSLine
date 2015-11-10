#' \code{\link{ramclustR}} - modified for MetMSLine
#' 
#' @param peakTable optional if pcaResult argument not supplied. Either a data.frame, 
#' full file path as a character string to a  .csv file of a peak table in the 
#' form observation (samples) in columns and variables (Mass spectral signals) 
#' in rows. If argument is not supplied a GUI file selection window will open 
#' and a .csv file can be selected.
#' @param obsNames character vector of observation (i.e. sample/ QC/ Blank) names 
#' to identify appropriate observation (sample) columns.
#' @param st Sigma t value controlling the rate at which feature similarity 
#' decays with respect to differences in retention times between feature pairs. 
#' When an xcms diff report is used as input, this parameter is by default set to 
#' half the peak width of all the peaks in the xcmsSet. Can be manually overridden 
#' to tune performance when necessary. Larger values are less selective, smaller 
#' values are more selective.
#' @param sr Sigma r value controlling the rate at which feature similarity 
#' decays with respect to correlational simiarity between a pair of features. 
#' numeric e.g. 0.5, Sigma r: the rate of similarity decay in response to 
#' correlational r . Can be manually overridden to tune performance. Larger 
#' values are less selective, smaller values are more selective.
#' @param maxt numeric e.g 20, the maximum time allowed for similarity between 
#' features to be greater than zero. For example, if two features are separated 
#' by more than 20 seconds, the similarity score is automatically set to zero, 
#' ensuring that they end up in different clusters. If you are importing data 
#' from another vendor, ensure that either you change the feature times to 
#' seconds (default in ramclustR) or you change this value to a value 
#' representing minutes, rather than seconds (e.g. 0.33). Used primarily to 
#' speed up compuation time by not calculating values above this threshold.
#' @param deepSplit boolean e.g. =FALSE, access to the deepSplit function in 
#' \code{\link{dynamicTreeCut}}. False generates fewer spectra, 
#' TRUE splits spectra more readily and thereby generates more spectra. 
#' For more detailed description, see documentation for the 
#' \code{\link{dynamicTreeCut}} package.
#' @param blocksize integer number of features processed in one block =2000, 
#' ramclustR makes attempts to minimize memory demands to maintain processing 
#' speed by avoiding memory constraints. It does so by processing the full 
#' similarity matrix in blocks of blocksize dimensions. This has two benefits: 
#' one is that a given PC is less likely to run out of memory during processing 
#' and the second is that we need not process blocks that contain not feature 
#' pairs with a delta rt of greater than maxt, thereby avoiding needless 
#' calculations and processing time.
#' @param minModuleSize the number of features required for a cluster to be 
#' reported. That is, if a feature is clustered in a group with less than 
#' minModuleSize features, it will not be exported. If you wish to analyze 
#' singletons, set minModuleSize=1.
#' @param linkage The linkage method used to perform heirarchical 
#' fastcluster-based heirarchical clustering. default="average".
#' 
#' @return a list containing 2 elements: 
#' 
#'  1. the original peak table with the results of retention time 
#'  correlation clustering included in additional column "RtCorrClust" as a 
#'  data frame.
#'  
#'  2. A data frame of the weighted means and details of the most intense 
#'  feature for each group. In addition for each group details of the best sample 
#'  to consider for subsequent MSn fragmentation analysis based on the most 
#'  top 5 most intense samples is provided. This "best" sample for MS/MS is also
#'  selected whilst concurrently minimising the number of samples to potentially 
#'  reinject for targetted MSn fragmentation analysis. 
#'  
#' @source \url{http://pubs.acs.org/doi/abs/10.1021/ac501530d} RAMClust: A Novel 
#' Feature Clustering Method Enables Spectral-Matching-Based Annotation for 
#' Metabolomics Data.
#' @export
ramClustMod <- function (peakTable=NULL, obsNames=NULL,
          st = 3, sr = 0.5, maxt = 20,  mult = 5, deepSplit = FALSE, blocksize = 2000, 
          hmax = 0.3, minModuleSize = 2, linkage = "average"){

# package dependencies
# require(xcms, quietly = TRUE)
# require(ff, quietly = TRUE)
# require(fastcluster, quietly = TRUE)
# require(dynamicTreeCut, quietly = TRUE)

if(is.null(obsNames)){
  stop('argument obsNames is missing with no default')
} 
# error handling or read from csv function
peakTable <- MetMSLine:::tableCheckRead(peakTable, stringsAsFactors=F)
# match obsNames to peak table colnames
obsIndx <- match(obsNames, colnames(peakTable))
# if less than all matched then stop
if(length(obsIndx) < length(obsNames)){
  stop(length(obsIndx), " of ", length(obsNames), 
       " observation names were matched in the peakTable column names, check the obsNames and peakTable column names")
}
# replace any NAs with zero
obsTable <- t(peakTable[, obsIndx])

# if contains rtmin and rtmax columns
  if(all(c("rtmin", "rtmax") %in% colnames(peakTable))){
     st <- median((peakTable[, "rtmax"] - peakTable[, "rtmin"])/2)
   }

times <- peakTable[, "rtmed"]
mzs <-  peakTable[, "mzmed"]
xcmsOrd <- order(times)
obsTable <- obsTable[, order(times)]
mzs <- mzs[order(times)]
times <- times[order(times)]
featnames <- paste(mzs, "_", times, sep = "")
dimnames(obsTable)[[2]] <- featnames
n <- ncol(obsTable)
vlength <- (n * (n - 1))/2
nblocks <- floor(n/blocksize)
ffmat <- ff::ff(vmode = "double", dim = c(n, n), init = 0)
gc()
eval1 <- expand.grid(0:nblocks, 0:nblocks)
names(eval1) <- c("j", "k")
eval1 <- eval1[which(eval1[, "j"] <= eval1[, "k"]), ]
bl <- nrow(eval1)
a <- Sys.time()
cat("\n", paste("calculating ramclustR similarity: nblocks = ", 
                bl))
cat("\n", "finished:")
RCsim <- function(bl) {
  cat(bl, " ")
  j <- eval1[bl, "j"]
  k <- eval1[bl, "k"]
  startc <- min((1 + (j * blocksize)), n)
  if ((j + 1) * blocksize > n) {
    stopc <- n
  }
  else {
    stopc <- (j + 1) * blocksize
  }
  startr <- min((1 + (k * blocksize)), n)
  if ((k + 1) * blocksize > n) {
    stopr <- n
  }
  else {
    stopr <- (k + 1) * blocksize
  }
  if (startc <= startr) {
    mint <- min(abs(outer(range(times[startr:stopr]), 
                          range(times[startc:stopc]), FUN = "-")))
    if (mint <= maxt) {
      temp1 <- round(exp(-(((abs(outer(times[startr:stopr], 
                                       times[startc:stopc], FUN = "-"))))^2)/(2 * 
                                                                                (st^2))), digits = 20)
      temp2 <- round(exp(-((1 - (pmax(cor(obsTable[, startr:stopr], 
                                          obsTable[, startc:stopc]))))^2)/(2 * (sr^2))), 
                     digits = 20)
      temp <- 1 - (temp1 * temp2)
      temp[which(is.nan(temp))] <- 1
      temp[which(is.na(temp))] <- 1
      temp[which(is.infinite(temp))] <- 1
      ffmat[startr:stopr, startc:stopc] <- temp
      rm(temp1)
      rm(temp2)
      rm(temp)
      gc()
    }
    if (mint > maxt) {
      ffmat[startr:stopr, startc:stopc] <- 1
    }
  }
  gc()
}
system.time(sapply(1:bl, RCsim))
b <- Sys.time()
cat("\n", "\n")
cat(paste("RAMClust feature similarity matrix calculated and stored:", 
          round(difftime(b, a, units = "mins"), digits = 1), "minutes"))
gc()
blocksize <- mult * round(blocksize^2/n)
nblocks <- floor(n/blocksize)
remaind <- n - (nblocks * blocksize)
RC <- vector(mode = "integer", length = vlength)
if(nblocks != 0){
pb <- txtProgressBar(min=0, max=nblocks, style=3)
}
for (k in 0:(nblocks)) {
  if(nblocks != 0){
  Sys.sleep(0.01)
  setTxtProgressBar(pb, k)
  }
  startc <- 1 + (k * blocksize)
  if ((k + 1) * blocksize > n) {
    stopc <- n
  }
  else {
    stopc <- (k + 1) * blocksize
  }
  temp <- ffmat[startc:nrow(ffmat), startc:stopc]
  temp <- temp[which(row(temp) - col(temp) > 0)]
  if (exists("startv") == FALSE) 
    startv <- 1
  stopv <- startv + length(temp) - 1
  RC[startv:stopv] <- temp
  gc()
  startv <- stopv + 1
  rm(temp)
  gc()
}
rm(startv)
gc()
message("Converting RAMClustR similarity matrix to a distance object...this
        is a potentially computationally intensive step please wait...")
flush.console()
RC <- structure(RC, Size = (n), Diag = FALSE, Upper = FALSE, 
                method = "RAMClustR", Labels = featnames, class = "dist")
gc()
c <- Sys.time()
cat("\n", "\n")
cat(paste("RAMClust distances converted to distance object:", 
          round(difftime(c, b, units = "mins"), digits = 1), "minutes"))
ff::delete.ff(ffmat)
rm(ffmat)
gc()
system.time(RC <- fastcluster::hclust(RC, method = linkage))
gc()
d <- Sys.time()
cat("\n", "\n")
cat(paste("fastcluster based clustering complete:", round(difftime(d, 
                                                                   c, units = "mins"), digits = 1), "minutes"))

message("cutting dendrogram dynamically...")
flush.console()
if (minModuleSize == 1) {
  clus <- dynamicTreeCut::cutreeDynamicTree(RC, maxTreeHeight = hmax, deepSplit = deepSplit, 
                                            minModuleSize = 2)
  sing <- which(clus == 0)
  clus[sing] <- max(clus) + 1:length(sing)
}
if (minModuleSize > 1) {
  clus <- dynamicTreeCut::cutreeDynamicTree(RC, maxTreeHeight = hmax, 
                                            deepSplit = deepSplit, 
                                            minModuleSize = minModuleSize)
}
gc()
rm(RC)

peakTable$RtCorrClusts <- clus[order(xcmsOrd)]

zeroIndx <- which(peakTable$RtCorrClusts != 0)
message("Calculating weighted mean for ", 
        length(unique(peakTable$RtCorrClusts[zeroIndx])), 
        " pseudospectra accounting for ", length(peakTable$RtCorrClusts[zeroIndx]), 
        " of ", nrow(peakTable), " total features ")
flush.console()
# calculate weighted.means pseudospectra and return max intensity feature
# and sample
wMeanPspec <- by(peakTable[zeroIndx, ], 
                 as.factor(peakTable$RtCorrClust[zeroIndx]), function(x){
                   # calc weights each variable
                   wts.tmp <- rowSums(x[, obsIndx])
                   wt.mean.tmp <- apply(x[, obsIndx], 2, function(y){
                     weighted.mean(y, wts.tmp)})
                   # concatenate max feature intensity and weighted means
                   maxFeat.indx <- which.max(wts.tmp)
                   # required columns from xcms diff
                   reqCols <- seq(1, ncol(peakTable), 1)
                   reqCols <- colnames(peakTable)[setdiff(reqCols, obsIndx)]
                   # calculate average intensity for max weight
                   meanInt <- mean(as.numeric(x[maxFeat.indx, obsIndx]))
                   names(meanInt) <- "meanInt"
                   # identify top 5 max Int samples 
                   top5int <- head(colnames(sort(x[maxFeat.indx, obsIndx], 
                                                 decreasing=T)), n=5)
                   names(top5int) <- paste0("mostIntSamp_", 1:5)
                   # sum features in pseudospectra
                   nFeatPspec <- nrow(x)
                   names(nFeatPspec) <- "nFeatPspec"
                   wt.mean.tmp <- c(x[maxFeat.indx, reqCols], meanInt, 
                                    nFeatPspec, top5int, wt.mean.tmp)
                   return(wt.mean.tmp)})
# convert by array to data frame
wMeanPspec <- data.frame(do.call(rbind, wMeanPspec), stringsAsFactors=F)
# frequency of most intense sample
mostIntIndx <- grep("mostIntSamp_", colnames(wMeanPspec))
freqMostInt <- sort(table(unlist(wMeanPspec[, mostIntIndx])), decreasing=T)
# add a column for the best sample to select for MS/MS
bestSampMSMS <- apply(wMeanPspec[, mostIntIndx], 1, function(x){
  mostIntSamp.tmp <- unlist(x)
  indx.tmp <- which(names(freqMostInt) %in% 
                      mostIntSamp.tmp)
  bestSamp.tmp <- names(freqMostInt)[indx.tmp[1]]
  return(bestSamp.tmp)})
# rearrange column order and add in best sample for MS/MS
sampIdLogi <- colnames(wMeanPspec) %in% obsNames 
RtCorrClust <- sort(unique(peakTable$RtCorrClust))
RtCorrClust <- RtCorrClust[RtCorrClust != 0]
wMeanPspec <- data.frame(RtCorrClust=RtCorrClust, 
                         wMeanPspec[, sampIdLogi == F], bestSampMSMS, 
                         wMeanPspec[, sampIdLogi == T], stringsAsFactors=F)
# coerce list to character
wMeanPspec <- data.frame(lapply(wMeanPspec, as.character), stringsAsFactors=F)
# return a list containing the original xcms diff report with rt corr clust
# info and the weighted mean data frame
return(list(diffRepAdd=peakTable, wMeanPspec=wMeanPspec))
# RC$featclus <- clus
# RC$frt <- times
# RC$fmz <- mzs
# RC$xcmsOrd <- xcmsOrd
# msint <- rep(0, length(RC$fmz))
# message("Calculating weighted means for each pseudospectral cluster...")
# flush.console()
# for (i in 1:ncol(obsTable)) {
#   msint[i] <- weighted.mean(obsTable[, i], obsTable[, i])
# }
# RC$msint <- msint
# # if (ExpDes[[2]]["MSlevs", 1] == 2) {
# #   msmsint <- rep(0, length(RC$fmz))
# #   for (i in 1:ncol(obsTable)) {
# #     msmsint[i] <- weighted.mean(data2[, i], data2[, i])
# #   }
# #   RC$msmsint <- msmsint
# # }
# clrt <- aggregate(RC$frt, by = list(RC$featclus), FUN = "mean")
# RC$clrt <- clrt[which(clrt[, 1] != 0), 2]
# clrtsd <- aggregate(RC$frt, by = list(RC$featclus), FUN = "sd")
# RC$clrtsd <- clrtsd[which(clrtsd[, 1] != 0), 2]
# RC$nfeat <- as.vector(table(RC$featclus)[2:max(RC$featclus)])
# RC$nsing <- length(which(RC$featclus == 0))
# e <- Sys.time()
# cat("\n", "\n")
# cat(paste("dynamicTreeCut based pruning complete:", round(difftime(e, 
#                                                                    d, units = "mins"), digits = 1), "minutes"))
# f <- Sys.time()
# cat("\n", "\n")
# cat(paste("RAMClust has condensed", n, "features into", max(clus), 
#           "spectra in", round(difftime(f, a, units = "mins"), digits = 1), 
#           "minutes", "\n"))
# RC$cmpd <- paste("C", 1:length(RC$clrt), sep = "")
# RC$ann <- RC$cmpd
# RC$annconf <- rep("", length(RC$clrt))
# RC$annnotes <- rep("", length(RC$clrt))
# RC$MSdata <- obsTable
# 
#   cat("\n", "\n", "... collapsing features into spectra")
#   wts <- colSums(obsTable[])
#   RC$SpecAbund <- matrix(nrow = nrow(obsTable), ncol = max(clus))
#   for (ro in 1:nrow(RC$SpecAbund)) {
#     for (co in 1:ncol(RC$SpecAbund)) {
#       RC$SpecAbund[ro, co] <- weighted.mean(obsTable[ro, 
#                                                   which(RC$featclus == co)], wts[which(RC$featclus == 
#                                                                                          co)])
#     }
#   }
#   dimnames(RC$SpecAbund)[[1]] <- dimnames(RC$MSdata)[[1]]
#   dimnames(RC$SpecAbund)[[2]] <- paste("C", 1:ncol(RC$SpecAbund), 
#                                        sep = "")
#  
#   g <- Sys.time()
#   cat("\n", "\n")
#   cat(paste("RAMClustR has collapsed feature quantities\n             into spectral quantities:", 
#             round(difftime(g, f, units = "mins"), digits = 1), 
#             "minutes", "\n"))
# 
# rm(obsTable)
# 
# #   if (length(dimnames(RC$SpecAbund)[[1]]) > length(unique(dimnames(RC$SpecAbund)[[1]]))) {
# #     RC$SpecAbundAve <- aggregate(RC$SpecAbund[, 1:ncol(RC$SpecAbund)], 
# #                                  by = list(dimnames(RC$SpecAbund)[[1]]), FUN = "mean", 
# #                                  simplify = TRUE)
# #     dimnames(RC$SpecAbundAve)[[1]] <- RC$SpecAbundAve[, 1]
# #     RC$SpecAbundAve <- as.matrix(RC$SpecAbundAve[, 2:ncol(RC$SpecAbundAve)])
# #     dimnames(RC$SpecAbundAve)[[2]] <- dimnames(RC$SpecAbund)[[2]]
# #     gc()
# #   }
# # 
# # gc()
# 
# return(RC)
}