#' peak/ co-variate table error handling and/ or read from character string
#' @param tmpTable  either a data.frame, full file path as a character string to a text file (.csv/.txt/.tsv) of a table 
#' If argument is not supplied a GUI file selection window will open and a text file (.csv/.txt/.tsv) can be selected.
#' @param type description of table type as character string for gui windows and error messages. 
#' @param ... additional arguments to \code{\link{fread}}.
tableCheckRead <- function(tmpTable=NULL, type="peak-picker", ...){
  
  if(is.null(tmpTable)){
  message("file chooser window opened...")
  flush.console()
  tmpTable <- tcltk::tclvalue(tcltk::tkgetOpenFile(filetypes = "{{comma delimited text file} {.csv}} {{All files} *}",
                                                    title=paste0("Select your ", type, " table")))
  # read in table if necessary
  tmpTable <- data.table::fread(tmpTable, header=T, ...)   
} else if(is.character(tmpTable)){
  # read in table from character string if necessary
  tmpTable <- data.table::fread(tmpTable, header=T, ...)  
} 
if(!is.data.frame(tmpTable)){
  stop(paste0(type, " argument is not a data.frame"))
}
return(tmpTable)
}