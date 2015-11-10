#' Generate target list (.csv) file from subset peak table
#' 
#' @param peakTable either a data.frame, full file path as a character string to a  .csv file of a peak table in the form observation (samples) in columns and
#' variables (Mass spectral signals) in rows. If argument is not supplied a GUI file selection window will open and a .csv file can be selected.
#' @param filename character full path file name to save the target list file to.
#' @param deltaMz numeric delta mass accuracy (ppm) for precursor selection.
#' @param deltaRt numeric delta retention time in minutes
#' @param isoWidth character isolation width for quadrupole. either 'Narrow (~1.3 m/z)' 
#' or 'Medium (~4 m/z)'.
#' @export
targMS2csv <- function(peakTable=NULL, filename=NULL, manufacturer="Agilent", deltaMz=20, deltaRt=1,
                       isoWidth = c("Narrow (~1.3 m/z)", "Medium (~4 m/z)")){
  # error handling or read from csv function
  if(is.null(filename)){
    stop('argument filename is missing with no default')
  }
  peakTable <- tableCheckRead(peakTable, stringsAsFactors=F)
  # check if table contains mzmed and rtmed column names
  if(any((c('mzmed', 'rtmed') %in% colnames(peakTable)) == F)){
    stop('peakTable must contain both the column names "mzmed" and "rtmed", denoting 
         the m/z ratio and retention time respectively')
  }
  if(manufacturer == 'Agilent'){
  #column names
  cnames <- c("On", "Prec. m/z", "Delta m/z (ppm)", "Z", "Prec. Type", 
              "Ret. time (min)", "Delta ret. time (min)", "Iso. width", 
              "Collision energy") 
  targTable <- matrix("", nrow = nrow(peakTable) + 2, ncol = length(cnames))
  targTable[1, 1] <- "AutoPreferredExcludeMSMSTable"
  targTable[2, ] <- cnames
  r <- 3:(nrow(peakTable) + 2)
  targTable[r, 1] <- "TRUE"
  targTable[r, 2] <- round(peakTable[, "mzmed"], 4)
  targTable[r, 3] <- deltaMz
  targTable[r, 4] <- 1
  targTable[r, 5] <- "Preferred" 
  targTable[r, 6] <-  round(peakTable[, "rtmed"]/60, 3)
  targTable[r, 7] <- deltaRt
  targTable[r, 8] <-   isoWidth
  idx <- order(as.numeric(targTable[r, 6]))
  targTable[r, ] <- targTable[r[idx], ]
  write.table(targTable, file = filename, append = FALSE, quote = FALSE, 
              sep = ",", eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
              col.names = FALSE)
  }
} # end function