#' automatic covariate (y-variable) type based univariate statistical test selection.
#' 
#' @description function determines whether a co-variate is two sample, continuous or up to
#' n distinct classes (three or more sample) and selects the statistical method to apply accordingly.
#' 
#'      
#' @param peakTable either a data.frame, full file path as a character string to a  .csv file of a peak table in the form observation (samples) in columns and
#' variables (Mass spectral signals) in rows. If argument is not supplied a GUI file selection window will open and a .csv file can be selected.
#' 
#' @param obsNames character vector of observation (i.e. sample/ QC/ Blank) names to identify appropriate observation (sample) columns.
#' @param  coVariate numeric, character or factor vector of covariates (y-variables).
#'       
#' @param continuous default = TRUE, if TRUE continuous variables are treated as
#'   such and correlation analysis is selected (spearman or pearson is selected
#'   based on whether the data is log-transformed or not see Logged argument).
#'   if FALSE continuous variables are split into a two sample variable.
#'       
#' @param minSampPerClass default =3, minimum number of samples for a unique class 
#'   else multiple class samples will be considered as continuous variables.
#'  (N.B. at least 5 for Mann-Whitney-U)
#'       
#' @param Logged default = TRUE, has the data already been log transformed, if TRUE
#'   parametric methods will be selected (i.e. t-test (two sample), ANOVA 
#'   (three or more sample), pearson correlation) if FALSE then non-parametric methods 
#'   will be selected (i.e. t-test (two sample), ANOVA (three or more sample), pearson 
#'   correlation). NB. if residuals from \code{\link{batchAdj}} contain negative
#'   values then the negatives will be shifted to the right with a constant value
#'    with the smallest lowest value equal to one prior to fold change calculation.
#' @param base numeric this will be used to exponentiate log transformed data
#' for mean/ median fold change calculation. 
#'     
#' @param MTC Multiple Testing Correction default is "none", see \code{\link{p.adjust.methods}} for
#'   details of options. ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none") 
#'
#' @return a named list (according to stat method selected) containing the 5 followning
#'     named elements :
#'     
#'     varType = type of variable identified.
#'     
#'     method = statistical method applied (either parametric or non-parametric 
#'     according to whether data log transformed).
#'     
#'     MTC = multiple testing correction method applied
#'     
#'     coVariate = the co-variate supplied, if continuous is False the high/low two sample
#'     variable created for any continuous variables will appear here.
#'      
#'     result = depending on type of variable supplied the result element will contain
#'     the following:
#'     
#'       1. "Continuous" 4 column data frame of:
#'        
#'                      i. adjusted/ multiple testing corrected p-/q-values 
#'                         (dependent on MTC method) 
#'                         
#'                      ii. raw p.value (from cor.test function) 
#'                      
#'                      iii. estimated correlation coefficient (from cor.test function) 
#'                      
#'                      iv. mean or median fold change value (dependent on whether
#'                          parametric or non-parametric method used). In the case
#'                          of continuous variables a two sample high-/low- variable 
#'                          is created (using the cut function) 
#'                          
#'       2. "two sample" 3 column data frame of:
#'        
#'                     i. adjusted/ multiple testing corrected p-/q-values 
#'                        (dependent on MTC method) 
#'                        
#'                     ii. raw p.value (from either t.test or wilcox.test function) 
#'                     
#'                     iii. mean or median fold change value (dependent on whether
#'                         parametric or non-parametric method used). 
#'                         
#'       3. "three or more sample" variable column (dependent on number of comparisons) 
#'            data frame of: 
#'                        
#'                     i. adjusted/ multiple testing corrected p-/q-values 
#'                     (dependent on MTC method).
#'                      
#'                     ii. raw p.value (from either ANOVA or wilcox.test function).
#'                      
#'                     iii - number pairwise class comparisons. 
#'                           mean or median fold change columns. 
#'                           All pairwise class comparisons of fold change are carried out. 
#'                           (dependent on whether parametric or non-parametric method used).  
#'                     
#'   raw p.values/ multiple testing corrected p-/q-values (dependent on method)
#'   and fold changes (if applicable).
#'   
#'   A class must contain a minimum number of samples to be considered.
#'   if continuous equals true continuous variables 
#' @export
coVarTypeStat <- function(peakTable=NULL, obsNames=NULL, coVariate=NULL, continuous=T, 
                          nMaxClasses=3, minSampPerClass=3, Logged=T, base=2,
                          MTC="none"){
  options(stringsAsFactors = FALSE)
  #error handling
  if(length(obsNames) != length(coVariate)){
    stop("the obsNames and coVariate vector lengths are different")
  }
  # error handling or read from csv function
  peakTable <- tableCheckRead(peakTable, stringsAsFactors=F)
  
  if(is.null(coVariate)){
    stop("argument coVariate is missing with no default")
  } 
  
  # convert to factor
  coVariate <- as.factor(coVariate)
  # n levels factor
  valueFreq <- table(coVariate)
  # test if only one factor level
  if(length(valueFreq) == 1){
  stop("argument coVariate contains only 1 value and therefore 1 class")
  }
    
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
    # convert chars to numeric
    obsTable <- apply(obsTable, 2, as.numeric)  
  
    # select coVariate type
    varType <- ifelse(length(valueFreq) == 2, "two sample", "three or more sample")
    varType <- ifelse(varType == "three or more sample" & min(valueFreq) <= minSampPerClass, 
                      "Continuous", varType)
    varType <- ifelse(varType == "three or more sample" & length(valueFreq) > nMaxClasses,
                      "Continuous", varType)
    
    # continuous variables first if continuous is False then change to two sample
    if(varType == "Continuous"){
      # if Continuous False then create two sample variable from continuous
      varType <- ifelse(continuous == F, "two sample", varType)
      # high low two sample using cut
      binom.cont <- as.numeric(cut(as.numeric(coVariate), 2))
      # if still Continuous  then do correlation
      if(varType == "Continuous"){
        coVariate <- as.numeric(coVariate)
        # if parametric pearson else if non-parametric spearman
        method.tmp <- ifelse(Logged == T, "pearson", "spearman")
        # apply correlation coeff and calc FC
        testCoVar <- as.data.frame(t(apply(obsTable, 1, function(Var){
          # cor.test then extract pvalue and estimate corr coef 
          cor.tmp <- cor.test(Var, coVariate, method = method.tmp)
          cor.tmp <- as.numeric(cor.tmp[c("p.value", "estimate")])
          # scale residuals for example if any values are negative for fold
          # change calculation by adding a constant
          if(any(Var <0)){
          Var <- Var + (abs(min(Var)) + 1)
          }  
          # fold change
          if(Logged == T){
            # if parametric mean fold change
            
            FC <- unique(ave(base ^ Var, binom.cont))
            FC <- FC[1]/ FC[2]
            } else {
            # if non-parametric median fold change
            FC <- unique(ave(Var, coVariate, FUN=median))
            FC <- FC[1]/ FC[2]
          }
          return(c(cor.tmp, FC))
        })))
        # add column names
        colnames(testCoVar) <- c("p.value", "corrCoeff", paste0("FoldChange ", 
                                                                paste0(unique(binom.cont), 
                                                                collapse=":")))
        testCoVar <- data.frame(Adj.p.value = p.adjust(testCoVar$p.value, 
                                                       method = MTC) , testCoVar)
        # return list of method and result
        testCoVar <- list(varType = varType, 
                          method = method.tmp,
                          MTC = MTC,
                          coVariate = coVariate,
                          result = testCoVar)
      } else {
        coVariate <- binom.cont
        valueFreq <- table(coVariate)
      }
    }   
    
    # decide univariate test
    if(varType == "two sample"){
      # test 
      testCoVar <- as.data.frame(t(apply(obsTable, 1, function(Var){
        if(Logged == T){
          # parametric
          p.value.tmp <- try(t.test(Var ~ coVariate)$p.value, silent=T)
          # catch error
          if(class(p.value.tmp) == 'try-error'){
          p.value.tmp <- 1  
          }
          # scale residuals for example if any values are negative for fold
          # change calculation by adding a constant
          if(any(Var < 0)){
            Var <- Var + (abs(min(Var)) + 1)
          }  
          # if parametric exponentiate mean fold change
          FC <- unique(ave(base ^ Var, coVariate))
          FC <- FC[1]/ FC[2]
        } else {
          # non-parametric
          P.value.tmp <- try(wilcox.test(Var ~ coVariate)[["p.value"]], silent=T)
          # catch error
          if(class(p.value.tmp) == 'try-error'){
            p.value.tmp <- 1  
          }
          # scale residuals for example if any values are negative for fold
          # change calculation by adding a constant
          if(any(Var <0)){
            Var <- Var + (abs(min(Var)) + 1)
          }
          # if non-parametric median fold change
          FC <- unique(ave(Var, coVariate, FUN=median))
          FC <- FC[1]/ FC[2]
        }
        
        return(c(p.value.tmp, FC))
      })))
      # add colnames
      colnames(testCoVar) <- c("p.value", paste0("FoldChange ", 
                                                 paste0(unique(coVariate), 
                                                        collapse=":")))
      testCoVar <- data.frame(Adj.p.value = p.adjust(testCoVar$p.value, 
                                                     method = MTC) , testCoVar)
      # return list of method and result
      testCoVar <- list(varType = varType, 
                        method = ifelse(Logged == T, "t.test", "wilcox.test"),
                        MTC = MTC,
                        coVariate = coVariate,
                        result = testCoVar)
      # if three or more sample ANOVA if parametric Mann-Whitney U if non-parametric
      # fold change is all comparisons
    } else if (varType == "three or more sample"){
      
      # indx for fold change all classes by expand grid
      indx <- expand.grid(1:length(valueFreq), 1:length(valueFreq))
      indx <- indx[indx[, 2] >= indx[, 1], , drop = F] # remove duplicates
      indx <- indx[indx[, 2] - indx[, 1]!=0, , drop = F] 
      
      # test 
      testCoVar <- as.data.frame(t(apply(obsTable, 1, function(Var){
        if(Logged == T){
          # parametric
          p.value.tmp <- try(aov(Var ~ coVariate), silent=T)
        # catch error
        if(class(p.value.tmp) == 'try-error'){
          p.value.tmp <- 1  
        } else {
          # extract p.value summary aov
          p.value.tmp <- summary(p.value.tmp)[[1]][1, 5]
        }
        } else {
          # non-parametric Mann-Whitney U
          P.value.tmp <- try(wilcox.test(Var, coVariate)[["p.value"]], silent=T)
        # catch error
        if(class(p.value.tmp) == 'try-error'){
          p.value.tmp <- 1  
        }
        }
        # fold change each class comparison  
        FC <- sapply(1:nrow(indx), function(classComp){
          name1.tmp <- names(valueFreq)[indx[classComp, 1]]
          name2.tmp <- names(valueFreq)[indx[classComp, 2]]
          if(Logged == T){
            # scale residuals for example if any values are negative for fold
            # change calculation by adding a constant
            if(any(Var < 0)){
              Var <- Var + (abs(min(Var)) + 1)
            }
            # if parametric mean fold change
            FC.tmp <- mean(base ^ Var[which(coVariate == name1.tmp)])/
              mean(base ^ Var[which(coVariate == name2.tmp)])
          } else {
            # scale residuals for example if any values are negative for fold
            # change calculation by adding a constant
            if(any(Var <0)){
              Var <- Var + (abs(min(Var)) + 1)
            }
            # if non-parametric median fold change
            FC.tmp <- median(Var[which(coVariate == name1.tmp)])/
              median(Var[which(coVariate == name2.tmp)])
          }
          return(FC.tmp)
        })
        return(c(p.value.tmp, FC))
      })))
      # Fold change comparison names
      FC.names <- paste0("FoldChange_Class_", names(valueFreq)[indx[, 1]], 
                         "_vs_", names(valueFreq)[indx[, 2]])
      # add colnames
      colnames(testCoVar) <- c("p.value", FC.names)
      # mult testing correction
      testCoVar <- data.frame(Adj.p.value = p.adjust(testCoVar$p.value, 
                                                     method = MTC) , testCoVar)
      # return list of method and result
      testCoVar <- list(varType = varType, 
                        method = paste0(length(valueFreq), "_way_",
                                        ifelse(Logged == T, "ANOVA", 
                                               "Mann_Whitney_U")),
                        MTC = MTC,
                        coVariate = coVariate,
                        result = testCoVar)
    }
    return(testCoVar)
} # end function