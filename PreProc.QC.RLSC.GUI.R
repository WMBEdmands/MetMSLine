
require(fgui)
require(tcltk2)


###Initial message and  checklist###
tkmessageBox(message = "You are commencing MetMSLine data processing.
MetMSLine is capable of operating with any peak-picking software output (.txt, .csv or .tsv format), however the following checklist must be STRICTLY adhered before commencing:
               
1. The 1st peak-picker output table column must be the EIC number or unique peak/feature identifier.

2. The 2nd peak-picker output table column must be the peak apex or median/average m/z value for the variable.
               
3. The 3rd peak-picker output table column must be the peak apex or median/average Retention time value for the variable.

4. Retention times must be in seconds and not minutes.

5. All of the sample intensity/observation columns follow after the MS variable information columns.
               
6. Most importantly the sample columns MUST be in acquisition order with the column conditioning QCs and intervally injected QCs in the correct positions." )

ReturnVal <- tkmessageBox(title = "study parent directory",
                          message = "1. First select your study parent directory", icon = "info", type = "ok")
study.dir<-tk_choose.dir(default = "", caption = "Select directory")

ReturnVal <- tkmessageBox(title = "Peak-picker output directory",
                          message = "2. Next select your peak-picker output directory", icon = "info", type = "ok")

Peak.picker.output.dir<-tk_choose.dir(default = "", caption = "Select directory")


PreProc.QC.RLSC<-function(Peak.picker.output.file="e.g. XCMS_output.tsv, Mzmine_output.txt etc.",first.QC.name="e.g. QC_1_CCQC",CCQC=10,QCInterval=5,MFC.norm=TRUE,smoother.span=0.2,RSD=30,alpha.gLog=1,scatter.plots=TRUE) {
  ###load package dependencies
  require (zoo) 
  
  first.QC.name<-gsub("-",".",first.QC.name)
  
  setwd(Peak.picker.output.dir)
  
  ###save parameters used for data generation###
  Parameters<-data.frame(Peak.picker.output=Peak.picker.output.file,CCQC,QC.Interval=QCInterval,MFC.norm,smoother.span=smoother.span,QC.RSD.cutoff=RSD)
 
  message("Reading peak.picker.output.file")#,quote=F)
  flush.console()
  #print("Reading peak.picker.output.file",quote=F)
  ##automatically identify file extension##
  if(substr(Peak.picker.output.file,nchar(Peak.picker.output.file)-2,nchar(Peak.picker.output.file))=="tsv"){
  X<-read.table(Peak.picker.output.file,sep="\t",header=T)
  } else if (substr(Peak.picker.output.file,nchar(Peak.picker.output.file)-2,nchar(Peak.picker.output.file))=="txt"){
    X<-read.table(Peak.picker.output.file,sep="\t",header=T)
  } else if (substr(Peak.picker.output.file,nchar(Peak.picker.output.file)-2,nchar(Peak.picker.output.file))=="csv"){
    X<-read.csv(Peak.picker.output.file,header=T)
  }
  message("Peak.picker.output.file read into R")#,quote=F)
  flush.console()
  ###create PreProc.QC.RLSC results subdirectory to keep everything tidy!####
  
  dirname<-paste(study.dir,"/PreProc.QC.RLSC.results/",sep="")
  dir.create(dirname)
  setwd(dirname)
  
  date<-Sys.time()
  date<-gsub("-",".",date)
  write.csv(Parameters,paste("Parameters",substr(date,1,10),".csv",sep=" "),row.names=FALSE)
  
last.CCQC<-(grep(first.QC.name,colnames(X))+CCQC)
#peakcolumn_no<-14+SGroups
peakcolumnsIndex<-as.numeric(1:(grep(first.QC.name,colnames(X))-1))##index of peak information columns
columnvector<-as.numeric(last.CCQC:length(X)) #all injections numeric vector
RAW_QCIndices<- seq(last.CCQC,length(X),QCInterval)


####zero filling####
  
  ##Remove rows containing missing values
  X<-na.omit(X)
  peakcolumns<-X[,peakcolumnsIndex]###subset peak variable information columns

  X.2<-X #create X.2 matrix
  X.2<-X.2[,columnvector]#remove peak columns and start dataframe from last column conditioning QC
  X.2[is.na(X.2)]<-0 ##replace n/a with zero

message("zero filling...")#,quote=F)
flush.console()

  Yzerofilled<-replace(X.2,X.2==0,XminNotzero<-min(apply(X.2, 2, function(x) min(x[x>0])))/2) #replace all zeros with half lowest prior to log transform
  X.2<-cbind(peakcolumns,Yzerofilled)#rebind columns

message("...done")#,quote=F)
flush.console()

  peakcolumn_no<-grep(first.QC.name,colnames(X)) #begin at last column conditioning QC
  QCIndices<- seq(peakcolumn_no,length(X.2),QCInterval)#index of QC samples
  samplevector<-as.numeric(peakcolumn_no:length(X.2)) #all injections numeric vector in Y
  samples<-samplevector[-match(QCIndices,samplevector)] # sample index
 
  samples.df<-X.2[,samples]
 
  #####MEDIAN FOLD CHANGE NORMALISATION######
  if(MFC.norm==TRUE){

  normalize.medFC <- function(mat) {
    # Perform median fold change normalisation
    #           X - data set [Variables & Samples]
    medSam <- apply(mat, 1, median)
    medSam[which(medSam==0)] <- 0.0001
    mat <- apply(mat, 2, function(mat, medSam){
      medFDiSmpl <- mat/medSam
      vec<-mat/median(medFDiSmpl)
      return(vec)
    }, medSam)
    return (mat)
  }
  
  samples.df<-normalize.medFC(samples.df)
  }


message("smoothing signal...")#,quote=F)
flush.console()

X.2[,samples]<-NA #replace sample columns with missing values
QCs<-X.2[,QCIndices]


#Lowess smoothing

Lowess<-apply(QCs,1,lowess,f=smoother.span) #apply Lowess smoothing on QCs, arguments can be added here
Lowess<- data.frame(matrix(unlist(Lowess), nrow=length(X$name), byrow=T)) #coerce list result to dataframe
Lowess<-Lowess[-c(1:length(QCs[1,]))] # remove X coordinates
X.2[,QCIndices]<-Lowess # replace QC values with LOESS smoothed

Z<-X.2[,samplevector] #create interpolation matrix

spline<-apply(Z,1,na.spline) #cubic spline interpolation for missing values
spline<-t(spline)

curve<-X.2[,samplevector]<-spline
  
###reinsert normalised or non-normalised samples###
if(CCQC!=0)
{
X[,(samples+(CCQC-1))]<-samples.df
} else {
  X[,(samples+(CCQC))]<-samples.df  
}
median.curve.df<-as.matrix(apply(curve,1,median))

median.curve.df<-median.curve.df[,rep(1,ncol(curve))]
  
curve.prop<-median.curve.df/curve

corrected<-X[,columnvector]*curve.prop

message("...done")
flush.console()

  if (sum(corrected[,1])==nrow(corrected)){
  stop("smoothing parameter f is too small 
       (increase the proportion of datapoints used for signal drift correction)")
  
}

  ###identify QC variables with zero standard deviation####
  QC_SD_zero<-apply(QCs,1,sd)
  X.2<-X.2[QC_SD_zero!=0,]
  Yzerofilled<-Yzerofilled[QC_SD_zero!=0,]
  QCs<-QCs[QC_SD_zero!=0,]
  X<-X[QC_SD_zero!=0,]
  peakcolumns<-peakcolumns[QC_SD_zero!=0,]
  curve<-curve[QC_SD_zero!=0,]
  corrected<-corrected[QC_SD_zero!=0,]
  
  ######
  ##Negative variable identification following correction
  
  negatives_mat<-data.matrix((corrected<0)*1)##identify variables with negative values following QC.LSC correction
  NegMat_rowsums<-rowSums(negatives_mat) ##take row sums 
  
  negative_variables<-corrected[NegMat_rowsums>=1,] ##subset negative variables
  peak_neg_variables<-peakcolumns[NegMat_rowsums>=1,]
  positive_variables<-corrected[NegMat_rowsums==0,] ##subset positive variables for gLog transform
  peak_pos_variables<-peakcolumns[NegMat_rowsums==0,]
  
  Raw_pos_variables<-Yzerofilled[NegMat_rowsums==0,]
  Raw_neg_variables<-Yzerofilled[NegMat_rowsums>=1,]
  
  Curve_pos_variables<-curve[NegMat_rowsums==0,]
  Curve_neg_variables<-curve[NegMat_rowsums>=1,]
                
  
#Reproducibility calculation
#Raw data

RAW_QCIndices<- seq(last.CCQC,length(X),QCInterval)
RAW_QCs<-X[,RAW_QCIndices]
RAW_QCs<-RAW_QCs[NegMat_rowsums==0,]

  Var_raw<-apply(RAW_QCs,1,var)#variance
  
  stddev_raw<-sqrt(Var_raw)#standard deviation
  
  average_raw<-apply(RAW_QCs,1,mean) #average
  
  RSD_raw<-(stddev_raw/average_raw)*100
  
  RSD_raw_below<-ifelse (RSD<=RSD_raw & RSD_raw>=0,0,1)
  
  Sum_RAW_Reprodfeatures<-sum(RSD_raw_below)#sum reproducible features before correction
    
  RSD_TUS_QCs_RAW<-(sqrt(var(apply(RAW_QCs,2,sum)))/mean(apply(RAW_QCs,2,sum)))*100

  RSD_TUS_QCs_RAW_below_threshold<-(sqrt(var(apply(RAW_QCs[RSD_raw_below==1,],2,sum)))/mean(apply(RAW_QCs[RSD_raw_below==1,],2,sum)))*100
  
# Reproducibility following correction

  QCsCorrIndices<-seq(1,length(samplevector),QCInterval) # index of QCs from corrected matrix
  QCscorr<-corrected[,QCsCorrIndices] # create matrix of corrected QCs
  QCscorr<-QCscorr[NegMat_rowsums==0,]
  
  Var_corr<-apply(QCscorr,1,var)#variance

  stddev_corr<-sqrt(Var_corr)#standard deviation

  average_corr<-apply(QCscorr,1,mean) #average

  RSD_corr<-((stddev_corr/average_corr)*100)

  RSD_corr_below<-ifelse (RSD<=RSD_corr&RSD_corr>=0,0,1)

  Sum_Corr_Reprodfeatures<-sum(RSD_corr_below)#sum reproducible features before correction
  
  RSD_TUS_QCscorr<-(sqrt(var(apply(QCscorr,2,sum)))/mean(apply(QCscorr,2,sum)))*100
  
  RSD_TUS_QCscorr_below_threshold<-(sqrt(var(apply(QCscorr[RSD_corr_below==1,],2,sum)))/mean(apply(QCscorr[RSD_corr_below==1,],2,sum)))*100
  
  ###Generalized Log transform###
message("Generalized Log transformation...")#,quote=F)
flush.console()

  log.corrected<-apply(positive_variables,2,function(x){log2((x+sqrt(x)^2+alpha.gLog^2)/2)}) #log transform on all samples and QCs

message("...done")
flush.console()

  ###fold change difference %CV following smoothing create scatterplots from greatly changed signals####
  
  foldCVIndex<-((RSD_raw/RSD_corr)>=2)
  
  sumfoldCVIndex<-sum(((RSD_raw/RSD_corr)>=2)*1)
    
  RawfoldCV<-as.data.frame(Raw_pos_variables[foldCVIndex,]) ##matrix >2fold reduction in CV% following curve
  
  QCcorrfoldCV<-as.data.frame(positive_variables[foldCVIndex,]) ##matrix >2fold reduction in CV% following curve
  
  QCcorrLogfoldCV<-as.data.frame(log.corrected[foldCVIndex,])
  
  peakfoldCV<-as.data.frame(peak_pos_variables[foldCVIndex,])
  
  CurvefoldCV<-as.data.frame(Curve_pos_variables[foldCVIndex,])
  
  RSD_rawfoldCV<-RSD_raw[foldCVIndex]
  
  RSD_corrfoldCV<-RSD_corr[foldCVIndex]
  
  tests_raw<-cbind(stddev_raw,average_raw,RSD_raw,RSD_raw_below)
  
  tests_corr<-cbind(stddev_corr,average_corr,RSD_corr,RSD_corr_below)
  
  Raw_data_RSD<-X[NegMat_rowsums==0,]
  

  Corrected<-cbind(peak_pos_variables,tests_raw,tests_corr,positive_variables)

  Corrected.LogT<-cbind(peak_pos_variables,tests_raw,tests_corr,log.corrected)
  
  Sum_reproducible_features<-cbind(Sum_RAW_Reprodfeatures,RSD_TUS_QCs_RAW,RSD_TUS_QCs_RAW_below_threshold,Sum_Corr_Reprodfeatures,RSD_TUS_QCscorr,RSD_TUS_QCscorr_below_threshold,sumfoldCVIndex)
 
  ###.csv file creation###
  
 # Curve<-cbind(peak_pos_variables,Curve_pos_variables)
#  write.csv (Curve,file="Curve.csv")
  write.csv(Corrected,"Corrected.csv",row.names=FALSE)
  write.csv (Corrected.LogT,file="Corrected.LogT.csv",row.names=FALSE)
  write.csv (Sum_reproducible_features,file="Sum_reproducible_features.csv",row.names=FALSE)
  #write.csv(Raw_data_RSD,"Raw_data_RSD.csv",row.names=FALSE)
  
if (scatter.plots==TRUE){  
  message("SAVING SCATTERPLOTS...")#,quote=F)
  flush.console()
  
  if(sumfoldCVIndex>1){
    ####Plot univariate scatterplot smoothing####
  QCdummyMindex<-seq(1,length(RawfoldCV),QCInterval) # dummy matrix of QC injection position for PCA modelling
  QCdummyM<-rep(1,ncol(RawfoldCV))
  QCdummyM[QCdummyMindex]<-2
  Plot_names<-as.character(paste("M",round(peakfoldCV[,2],digits=4),"T",round(peakfoldCV[,3],digits=1),sep=""))
  
  dirname<-paste(dirname,"QC.LSC.Scatterplots")
  dir.create(dirname)
  setwd(dirname)
  
  pb<-txtProgressBar(min=0,max=nrow(RawfoldCV),style=3)#,width=300)#title="Scatter plot progress bar"
  
  
  for (i in 1:nrow(RawfoldCV)) { 

    Sys.sleep(0.1)
   setTxtProgressBar(pb,i)
    # setTkProgressBar(pb,i,label=paste( round(i/nrow(RawfoldCV)*100, 0),"% done"))
    flush.console()
    
    Titlename<-Plot_names[i]
    RSD_raw_plot<-as.character(round(RSD_rawfoldCV[i],digits=2))
    RSD_corr_plot<-as.character(round(RSD_corrfoldCV[i],digits=2))
       AQ_order<-1:length(RawfoldCV[1,])
  
    ##Raw data## showing LOWESS curve
    
    png(paste(Titlename,".png"),width=1200,height=1200,res=275)
    plot(AQ_order,RawfoldCV[i,], main=paste("QC.LSC",Titlename),sub=paste("CV%",RSD_raw_plot),xlab="Acquisition_order", ylab="MS_signal_intensity",pch=19,col=c("black","red")[QCdummyM])
   
    points(AQ_order,CurvefoldCV[i,],col="blue",cex=0.25)
   
    legend("topleft",pch=c(19,19,17),legend=c("Samples","QCs","LSC.curve"),col=c("black","red","blue"))
  
    dev.off()
    
    ####CORRECTED FEATURE FIGURES###
   
    png(paste(Titlename,"CORR",".png"),width=1200,height=1200,res=275)
   
    plot(AQ_order,QCcorrfoldCV[i,], main=paste("QC.LSC",Titlename,"Corrected"),sub=paste("CV%",RSD_corr_plot,"(smooth.span=",as.character(round(smoother.span,digits=2)),")"), xlab="Acquisition_order", ylab="Corrected_signal",pch=19,col=c("black","red")[QCdummyM])
        
    dev.off()    
  }
  close(pb)
  message("...done")
  flush.console()
 }
}
  if(sumfoldCVIndex<2){print("smoothing parameter f too large
                              (decrease the proportion of datapoints used for signal drift correction)")}

message("Pre-processing finished please examine results")#,quote=F)
flush.console()

}

guiv(PreProc.QC.RLSC,
     argText=list(MFC.norm=c("Median Fold change normalisation ? "),
                  Peak.picker.output.file=c("Peak-picking software output table (.csv, .txt or .tsv) ? "),
                  first.QC.name=c("What is the precise name of your first column conditioning QC ? "),
                  CCQC=c("Number of Column Conditioning QC samples at the beginning of acquisition ?"),
                  QCInterval=c("Quality control injection interval ? "),
                  smoother.span=c("Smoother span for LOWESS signal attenuation smoothing ? "),
                  RSD=c("Relative Standard Deviation cut-off for pooled QC signal filtration ? "),
                  alpha.gLog=c("alpha value for generalized Log transformation ? "),
                  scatter.plots=c("Output Scatter plots showing effect of LOWESS smoothing ? ")),
    argOption=list(MFC.norm=c("FALSE","TRUE"),scatter.plots=c("TRUE","FALSE")),
    argSlider=list(CCQC=c(0,30,1),QCInterval=c(1,40,1),smoother.span=c(0.01,1,0.01),RSD=c(0,50,1),
                   alpha.gLog=c(0.1,3,0.1)),helps=NULL)
   
###END###