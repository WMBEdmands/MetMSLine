#' manhattan plot function
#'
#' @param x character string of a single column name representing the x axis (e.g. order or mass-to charge ratio).
#' @param y character string of a single column name (if data argument supplied) containing the statistical test p-values the y axis.
#' @param data data.frame containing the x and y columns
#' @param cols numeric vector of integers, hexadecimal colours or names of colours for data points.
#' @param plotOutputFile character string full path to write output image file (see\link{\code{ggsave}}).
#' @param minPval numeric minimum p-value to plot significance line and change point shapes.
#' @param labelFontSize numeric label and title font size (supplied to \link{\code{theme_bw}})
#'  
#' @param ... additional arguments to \link{\code{ggsave}}.
#' 
#' @export
manhattanPlot <- function(x=NULL, y=NULL, data=NULL, cols=NULL,
                          plotOutputFile=NULL, minPval=0.05, labelFontSize=15,
                          ...){
  # error handling
  if(is.null(x) | is.null(y)){
    stop('an argument must be supplied for both x and y...\n')  
  } 
  if(is.null(data)){
    stop('the argument data is missing with no default and must be a data.frame...\n') 
  }

  if(!require(ggplot2)){
    message('attempting to install ggplot2...\n')
    flush.console()
    install.packages('ggplot2')
    if(!require(ggplot2)){
      stop('The package ggplot2 must be installed to utilize this function...\n')
    }
  } 
  
  colnames(data)[colnames(data) == x] <- 'x'
  colnames(data)[colnames(data) == y] <- 'y'
  shapes <- (data[, 'y'] <= minPval) + 1
  
  data[, 'y'] <- -log10(data[, 'y'])
  
  if(is.null(cols)){
  cols <- shapes
  cols <- ifelse(cols == 2, 'red', 'black')
  }
  shapes <- ifelse(shapes == 2, 1.5, 1)
  g_plotTmp <- ggplot(aes(x=x, y=y), data=data) +
    geom_point(colour=as.factor(cols), fill=as.factor(cols), size=shapes) +
    xlab(x) +
    ylab("-Log10 (q-value)") +
    ggtitle(paste0(y, '\n(n=', sum(shapes == 1.5), ' <= ', minPval, ')')) + 
    guides(fill=F) +
    guides(colour=F) +
    guides(size=F) +
    guides(shape=F) +
    geom_hline(yintercept=-log10(minPval), colour="red", size=1)+
    theme_bw(labelFontSize)+
    theme(plot.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border=element_rect(colour="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
   
  if(!is.null(plotOutputFile)){
    ggsave(plotOutputFile, plot=g_plotTmp, ...)
  }
  return(g_plotTmp)
} # end function