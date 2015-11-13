## ---- eval=F-------------------------------------------------------------
#  install.packages(c( 'ff', 'dynamicTreeCut', 'data.table', 'reshape2', 'foreach'))

## ---- include=F----------------------------------------------------------
library(MetMSLine)

## ---- collapse=TRUE------------------------------------------------------
# read in table of synthetic MS1 metabolomic profiling data
peakTable <- system.file("extdata", "synthetic_MS_data.csv", package = "MetMSLine")
peakTable <- read.csv(peakTable, header=T, stringsAsFactors=F)


# load synthetic co-variates table in comma delimited csv file
coVariates <- system.file("extdata", "synthetic_coVar_data.csv", package = "MetMSLine")
coVariates <- read.csv(coVariates, header=T)

## ------------------------------------------------------------------------
colnames(peakTable)

## ------------------------------------------------------------------------
# observation names (i.e. sample names) by regular expression (?grep)
# all observation names
obsNames <- colnames(peakTable)[grep("QC_|sample_|blank_", colnames(peakTable))]
# print all the observation names
obsNames

## ------------------------------------------------------------------------
qcNames <- colnames(peakTable)[grep("QC_", colnames(peakTable))] 
# print the QC names
qcNames

## ------------------------------------------------------------------------
sampNames <- colnames(peakTable)[grep("sample_", colnames(peakTable))]
# print the sample names
sampNames

## ------------------------------------------------------------------------
blankNames <- colnames(peakTable)[grep("blank_", colnames(peakTable))]
# print the blank names
blankNames

## ------------------------------------------------------------------------
# detect number of cores using parallel package
nCores <- parallel::detectCores()
# conduct LC-MS data preprocessing
preProc_peakTable <- preProc(peakTable, obsNames, sampNames, qcNames, blankNames,
                             cvThresh=30, nCores=nCores)

## ------------------------------------------------------------------------
tail(colnames(preProc_peakTable), n=3)

