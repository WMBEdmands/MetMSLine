PreProc.QC.RLSC<-function(X="XCMS_output.tsv",CCQC=10,SGroups=3,QCInt=10,MFC.norm=TRUE,f=1/5,RSD=30,a=1,scatter.plots=TRUE,wd="D:\\R_data_processing\\STUDY NAME\\",XCMS_dir="D:\\R_data_processing\\STUDY NAME\\XCMS\\") {
# Performs QC based LOWESS curve signal correction. X-XCMS diffreport, CCQC - number of column conditioning QCs, 
# QCInt - QC injection interval (ie. every 4th sample), f = the smoother span (proportion of points in the plot which influence the smooth at each value), larger values = more smoothness.
# a<-1 alpha generalized log transform
# returns matrix of original, curve values and corrected signal intensity data #  
  
  setwd(XCMS_dir)
  
  ###save parameters used for data generation###
  Parameters<-data.frame(XCMS.output=X,CCQC,SGroups,QC.Interval=QCInt,MFC.norm,smoother.span=f,QC.RSD.cutoff=RSD)
  
  ##automatically identify file extension##
  if(substr(X,nchar(X)-2,nchar(X))=="tsv"){
  X<-read.table(X,sep="\t",header=T)
  }
  if(substr(X,nchar(X)-2,nchar(X))=="txt"){
    X<-read.table(X,sep="\t",header=T)
  }
  if(substr(X,nchar(X)-2,nchar(X))=="csv"){
    X<-read.csv(X,header=T)
  }
  
  ###create PreProc.QC.RLSC results subdirectory to keep everything tidy!####
  
  dirname<-paste(wd,"PreProc.QC.RLSC.results\\",sep="")
  dir.create(dirname)
  setwd(dirname)
  
  date<-Sys.time()
  date<-gsub("-",".",date)
  write.csv(Parameters,paste("Parameters",substr(date,1,10),".csv",sep=" "),row.names=FALSE)
    
require (zoo)
last.CCQC<-CCQC+14+SGroups
XCMScolumn_vector<-14+SGroups
XCMScolumnsIndex<-as.numeric(1:XCMScolumn_vector)##index of XCMS information columns
XCMScolumns<-X[,XCMScolumnsIndex]###subset XCMS variable information columns
columnvector<-as.numeric(last.CCQC:length(X)) #all injections numeric vector
  RAW_QCIndices<- seq(last.CCQC,length(X),QCInt)


####zero filling####
  
  ##Remove rows containing missing values
  X<-na.omit(X)
  
  X.2<-X #create X.2 matrix
  X.2<-X.2[,columnvector]#remove XCMS columns and start dataframe from last column conditioning QC
  X.2[is.na(X.2)]<-0 ##replace n/a with zero
  Yzerofilled<-replace(X.2,X.2==0,XminNotzero<-min(apply(X.2, 2, function(x) min(x[x>0])))/2) #replace all zeros with half lowest prior to log transform
  X.2<-cbind(XCMScolumns,Yzerofilled)#rebind columns
  

  XCMScolumn_vector<-XCMScolumn_vector+1 #begin at last column conditioning QC
  QCIndices<- seq(XCMScolumn_vector,length(X.2),QCInt)#index of QC samples
  samplevector<-as.numeric(XCMScolumn_vector:length(X.2)) #all injections numeric vector in Y
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
  
X.2[,samples]<-NA #replace sample columns with missing values
QCs<-X.2[,QCIndices]
  

#Lowess smoothing

Lowess<-apply(QCs,1,lowess,f=f) #apply Lowess smoothing on QCs, arguments can be added here
Lowess<- data.frame(matrix(unlist(Lowess), nrow=length(X$name), byrow=T)) #coerce list result to dataframe
Lowess<-Lowess[-c(1:length(QCs[1,]))] # remove X coordinates
X.2[,QCIndices]<-Lowess # replace QC values with LOESS smoothed

Z<-X.2[,samplevector] #create interpolation matrix

spline<-apply(Z,1,na.spline) #cubic spline interpolation for missing values
spline<-t(spline)

curve<-X.2[,samplevector]<-spline
  
###reinsert normalised or non-normalised samples###
X[,(samples+(CCQC-1))]<-samples.df
corrected<-X[,columnvector]/X.2[,samplevector]

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
  XCMScolumns<-XCMScolumns[QC_SD_zero!=0,]
  curve<-curve[QC_SD_zero!=0,]
  corrected<-corrected[QC_SD_zero!=0,]
  
  ######
  ##Negative variable identification following correction
  
  negatives_mat<-data.matrix((corrected<0)*1)##identify variables with negative values following QC.LSC correction
  NegMat_rowsums<-apply(negatives_mat,1,sum) ##take row sums 
  
  negative_variables<-corrected[NegMat_rowsums>=1,] ##subset negative variables
  XCMS_neg_variables<-XCMScolumns[NegMat_rowsums>=1,]
  positive_variables<-corrected[NegMat_rowsums==0,] ##subset positive variables for gLog transform
  XCMS_pos_variables<-XCMScolumns[NegMat_rowsums==0,]
  
  Raw_pos_variables<-Yzerofilled[NegMat_rowsums==0,]
  Raw_neg_variables<-Yzerofilled[NegMat_rowsums>=1,]
  
  Curve_pos_variables<-curve[NegMat_rowsums==0,]
  Curve_neg_variables<-curve[NegMat_rowsums>=1,]
                
  
#Reproducibility calculation
#Raw data

RAW_QCIndices<- seq(last.CCQC,length(X),QCInt)
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

  QCsCorrIndices<-seq(1,length(samplevector),QCInt) # index of QCs from corrected matrix
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
  
  log.corrected<-apply(positive_variables,2,function(x){log2((x+sqrt(x)^2+a^2)/2)}) #log transform on all samples and QCs
  
  ###fold change difference %CV following smoothing create scatterplots from greatly changed signals####
  
  foldCVIndex<-((RSD_raw/RSD_corr)>=2)
  
  sumfoldCVIndex<-sum(((RSD_raw/RSD_corr)>=2)*1)
    
  RawfoldCV<-as.data.frame(Raw_pos_variables[foldCVIndex,]) ##matrix >2fold reduction in CV% following curve
  
  QCcorrfoldCV<-as.data.frame(positive_variables[foldCVIndex,]) ##matrix >2fold reduction in CV% following curve
  
  QCcorrLogfoldCV<-as.data.frame(log.corrected[foldCVIndex,])
  
  XCMSfoldCV<-as.data.frame(XCMS_pos_variables[foldCVIndex,])
  
  CurvefoldCV<-as.data.frame(Curve_pos_variables[foldCVIndex,])
  
  RSD_rawfoldCV<-RSD_raw[foldCVIndex]
  
  RSD_corrfoldCV<-RSD_corr[foldCVIndex]
  
  tests_raw<-cbind(stddev_raw,average_raw,RSD_raw,RSD_raw_below)
  
  tests_corr<-cbind(stddev_corr,average_corr,RSD_corr,RSD_corr_below)
  
  Raw_data_RSD<-X[NegMat_rowsums==0,]
  
  QC_Corrected<-cbind(XCMS_pos_variables,tests_raw,tests_corr,log.corrected)
  
  Sum_reproducible_features<-cbind(Sum_RAW_Reprodfeatures,RSD_TUS_QCs_RAW,RSD_TUS_QCs_RAW_below_threshold,Sum_Corr_Reprodfeatures,RSD_TUS_QCscorr,RSD_TUS_QCscorr_below_threshold,sumfoldCVIndex)
 
  ###.csv file creation###
  
  Curve<-cbind(XCMS_pos_variables,Curve_pos_variables)
  write.csv (Curve,file="Curve.csv")
  write.csv (QC_Corrected,file="QC_Corrected.csv",row.names=FALSE)
  write.csv (Sum_reproducible_features,file="Sum_reproducible_features.csv",row.names=FALSE)
  write.csv(Raw_data_RSD,"Raw_data_RSD.csv",row.names=FALSE)
  
if (scatter.plots==TRUE){  
  if(sumfoldCVIndex>1){
    ####Plot univariate scatterplot smoothing####
  QCdummyMindex<-seq(1,length(RawfoldCV),QCInt) # dummy matrix of QC injection position for PCA modelling
  QCdummyM<-rep(1,ncol(RawfoldCV))
  QCdummyM[QCdummyMindex]<-2
  Plot_names<-as.character(XCMSfoldCV[,2])
  
  dirname<-paste(dirname,"QC.LSC.Scatterplots")
  dir.create(dirname)
  setwd(dirname)
  
  for (i in 1:nrow(RawfoldCV)) { 
    
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
   
    plot(AQ_order,QCcorrfoldCV[i,], main=paste("QC.LSC",Titlename,"Corrected"),sub=paste("CV%",RSD_corr_plot,"(smooth.span=",as.character(round(f,digits=2)),")"), xlab="Acquisition_order", ylab="Corrected_signal",pch=19,col=c("black","red")[QCdummyM])
        
    dev.off()    
  }
 }
}
  if(sumfoldCVIndex<2){print("smoothing parameter f too large
                              (decrease the proportion of datapoints used for signal drift correction)")}
  
}

###END###