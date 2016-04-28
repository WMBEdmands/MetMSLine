#' Generate cytoscape map files from peakTable
#'
#' @param peakTable either a data.frame, full file path as a character string to a  .csv file of a peak table in the form observation (samples) in columns and
#' variables (Mass spectral signals) in rows. If argument is not supplied a GUI file selection window will open and a .csv file can be selected.
#' @param obsNames character vector of observation (i.e. sample/ QC/ Blank) names to identify appropriate observation (sample) columns.
#' @param outDir character full path of cytoscape file output directory. The function
#' will create an additional subdirectory within this output called "cytoScapeFiles" into
#' which files will be saved.
#' @param EICnos character name of column in peakTable of unique feature numbers. 
#' This unique number (e.g. extracted ion chromatogram number) will be used to 
#' merge the peak table with the node attributes file.
#' @param fileNameId character additional identificatory character string to 
#' add to the cytoscape output file names. default = "". e.g. "blankSubtracted"
#' or "filtered".
#' @param corrThresh correlation coefficient threshold to group features within
#'  a retention time cluster.
#' @param corrMethod character correlation method see \code{\link{cor}} for details. default "spearman".
#' @param delta numeric maximum p-value (following multiple testing correction) above 
#' which the null hypothesis (no correlation) is rejected.
#' @param MTC character Multiple Testing Correction default is "BH", see \code{\link{p.adjust.methods}} for
#' details of options. ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none").
#' Any p-values after multiple testing correction above the value of delta will have their
#' corresponding correlation coefficents replaced with zero.
#' 
#' @return writes 3 text files: network (.sif), node and edge attribute files as 
#' plain text (.txt) files in the output directory (see outDir).
#' 
#' @export
cytoScapeMap <- function(peakTable=NULL, obsNames=NULL, outDir=NULL, EICnos=NULL,
                         fileNameId="", corrThresh=0.6, corrMethod="spearman", 
                         delta=0.05, MTC="BH"){
  
  # error handling
  if(is.null(obsNames)){
    stop('argument obsNames is missing with no default')
  } else if(is.null(outDir)){
    stop('argument outDir is missing with no default')
  } else if(is.null(EICnos)){
    stop('argument EICnos is missing with no default')
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
    stop("peakTable observations contain NAs or zeros unable to log transform, use the ?zeroFill function")
  }
  # calc correlation matrix
  message("Calculating correlation matrix for ", nrow(obsTable), " features")
  flush.console()
  cor.prob <- function(X, dfr = nrow(X) - 2, Method=corrMethod) {
    R <- cor(X, method=Method)
    r2 <- R^2
    Fstat <- r2 * dfr / (1 - r2)
    P <- 1 - pf(Fstat, 1, dfr)
    list(R, P)
  }
  # returns list R and P values
  corProbRes <- cor.prob(t(obsTable), Method=corrMethod)
  cor.m <- corProbRes[[1]]
  cor.m[apply(corProbRes[[2]], 2, p.adjust, 
              n=sum(lower.tri(cor.m)), method=MTC) > delta] <- 0
  colnames(cor.m)  <- peakTable[, EICnos]
  row.names(cor.m) <- peakTable[, EICnos]
  # replace upper tri with zero
  cor.m[upper.tri(cor.m, diag=T)] <- 0
  # ID features above below corrThresh
  sif <- apply(cor.m, 2, function(x){
    pos.tmp <- which(x > corrThresh)
    neg.tmp <- which(x < -corrThresh)
    if(length(pos.tmp) > 0){
    names(pos.tmp) <- x[pos.tmp]
    }
    if(length(neg.tmp) > 0){
    names(neg.tmp) <- x[neg.tmp]
    }
    return(list(pos=pos.tmp, neg=neg.tmp))
  })
  # melt list result
  sif.df <- reshape2::melt(sif)
  # add correlation value
  sif.df[, 4] <- as.numeric(gsub('.+pos\\.|.+neg\\.', '', names(unlist(sif))))
  # create sif file names
  sif.df[, 5] <- names(sif)[sif.df[, 1]]
  # n nodes and edges
  nodeAttr <- data.frame(name=as.numeric(unique(c(sif.df[, 3], sif.df[, 5]))))
  message(nrow(nodeAttr), 
          " nodes with ", nrow(sif.df), " edges identified at a corrThresh of ", 
          corrThresh)
  flush.console()
  # add no corr features
  noCorrFeat <- setdiff(peakTable[, EICnos], unique(c(sif.df[, 3], sif.df[, 5])))
  nodeAttr <- rbind(nodeAttr, data.frame(name=as.numeric(noCorrFeat)))
  nodeAttr <- merge(nodeAttr, peakTable, by.x="name", by.y=EICnos)
  nodeAttr <- nodeAttr[, -unique(unlist(lapply(obsNames, grep, colnames(nodeAttr))))]
  idCols <- data.frame(name=nodeAttr$name, matrix("", ncol=2, nrow=nrow(nodeAttr)))
  colnames(idCols)[2:3] <- c("Id", "Comments") 
  nodeAttr$name <- NULL
  nodeAttr <- cbind(idCols, nodeAttr)
  # write sif file
  cytoScape_dir <- paste0(outDir, '/cytoScapeFiles')
  suppressWarnings(dir.create(cytoScape_dir))
  res.dir <- paste0(cytoScape_dir, 
                    "/cytoScape_", corrThresh, "_n", length(obsNames), 
                    "_samples", "_", fileNameId)
  message("writing results output sif, node and edge attribute files to mzXML file directory :\n", 
          dirname(res.dir))
  flush.console()
  # write sif
  # add in features with no correlation above thresh
  noCorrFeatM <- matrix(0, nrow=length(noCorrFeat), ncol=5)
  colnames(noCorrFeatM) <- colnames(sif.df)
  noCorrFeatM[, 3] <- noCorrFeat
  noCorrFeatM[, 5] <- noCorrFeat
  sif.df <- rbind(sif.df, noCorrFeatM)
  write.table(paste0(sif.df[, 3], " \t tm \t ", sif.df[, 5], "\n"),
              paste0(res.dir, ".sif"), col.names=F, quote=F, row.names=F)
  # create node names
  sif.df[, 6] <- paste0(sif.df[, 3], " (tm) ", sif.df[, 5])
  
  edgeAttr <- data.frame(name=sif.df[, 6], correl_direction=sif.df[, 2],
                         corr=sif.df[, 4], corrRounded=abs(round(sif.df[, 4], digits=1)))
  # write edge attributes
  write.table(edgeAttr,
              paste0(res.dir,"_edgeAttr.txt"), sep="\t", row.names=F, quote=F)
  # write node attribute table
  write.table(nodeAttr,
              paste0(res.dir,"_nodeAttr.txt"), sep="\t", row.names=F, quote=F)
  # return correlation matrix
  #return(cor.m)
} 