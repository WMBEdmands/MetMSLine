
require(fgui)
require(tcltk2)

ReturnVal <- tkmessageBox(title = "Auto PCA output directory",
                          message = "Select the Auto PCA output directory", icon = "info", type = "ok")
wd<-tk_choose.dir(default = "", caption = "Select directory")


Auto.MV.Regress<-function( X="PCA.outliers.removed.csv", Yvar="Y.outliers.removed.csv", non.zero=2, Corr.thresh=0.3, pvalue=0.05, p.adjust.methods="none",  heatmap=TRUE, hclust.method="complete",dist.m.method="euclidean",Clust.ppm=10, Clust.RT.tol=2, mode="negative", box.prop=0.2,   Yunits="Y_units",  HMDBtol=0.005){
  
###load package dependencies
  require(gplots)
  require(ggplot2)
  require(sfsmisc)
  require(RColorBrewer)
  
setwd(wd)
###change arguments to lower case so not mistaken for function calls
if(hclust.method=="WARD.D")
{
  hclust.method<-"ward.D"
} else if (hclust.method=="WARD.D2")
{
  hclust.method<-"ward.D2"
} else {
hclust.method<-tolower(hclust.method)
}
dist.m.method<-tolower(dist.m.method)
p.adjust.methods<-tolower(p.adjust.methods)

message("Reading Auto.PCA output file...PLEASE WAIT")#,quote=F)
flush.console()

Samples<-read.csv(X,header=T)

message("...Done")#,quote=F)
flush.console()

###identify and store as objects EIC/ident, m/z and RT column names
EIC.column.name<-colnames(Samples)[1]
mzmed.column.name<-colnames(Samples)[2]
RTmed.column.name<-colnames(Samples)[3]

Y<-as.data.frame(read.csv(Yvar,header=T,row.names=1))

Ycolnames<-as.character(colnames(Y))

Yrownames<-as.character(row.names(Y))
####create Auto.MV.Regress subdirectory to keep everything tidy!####

wd<-paste(substr(wd,1,nchar(wd)-17),"/Auto.MV.Regress.results","/",sep="")

dir.create(wd)

setwd(wd)

###save all parameters used in a dated .csv file for future reference###
Parameters<-data.frame(X,Yvar,p.adjust.methods,box.plot.proportions=box.prop,
                       hclust.method,Pearson.corr.thresh=Corr.thresh,pvalue.cutoff=pvalue,
                       Cluster.RT.tol=Clust.RT.tol,Mass.Accuracy.ClusterID=Clust.ppm,Min.pos.values.Y=round(((nrow(Y)/100)*non.zero),digits=0))

date<-Sys.time()
date<-gsub("-",".",date)
write.csv(Parameters,paste("Parameters",substr(date,1,10),".csv",sep=" "),row.names=FALSE)


XCMScolumnsIndex<-c(1:which(colnames(Samples)=="RSD_corr_below"))
  
XCMScolumns<-Samples[,XCMScolumnsIndex]

Samples<-Samples[,-XCMScolumnsIndex]

###check if Y variable names are the same as sample column names if not stop the process###
if(length(which(colnames(Samples) %in% row.names(Y)==T)) != ncol(Samples))
{
  tkmessageBox(message ="The Y variable sample names in the 1st column do not match the 
PCA_outliers_removed.csv sample names! 

Please make sure the names match in the .csv files, change if necessary and start again.",title="ERROR")
  
  stop()
}
#SampleIndex<-colnames(t(Y)) %in% colnames(Samples) 

#Y<-data.frame(Y[SampleIndex,])

Yindices<-as.data.frame(Y>0) # all non zero Y-variables for regression

Yzeroindex<-as.data.frame(Y==0)

Yindices_sum<-as.data.frame(apply(Yindices,2,sum)) # sum of all non-zero Y-variables 

Yabove_10<-which(Yindices_sum>=((nrow(Y)/100)*non.zero)) # index all Y-variables with less than minimum proportion of non-zeros

Yindices<-as.data.frame(Yindices[,Yabove_10]) # remove all Y-variables from Yindices

Y<-as.data.frame(Y[,Yabove_10])

colnames(Y)<-Ycolnames[Yabove_10]

row.names(Y)<-Yrownames

##########################################################################################

results <- data.frame() # empty data.frame for storage of results above threshold
above_threshold_df<-data.frame()
box.plot.scores<-as.numeric()
Biomarker_numbers.df<-data.frame()


if (ncol(Y)>0){
  
  for (k in 1:ncol(Y)){
   
  foldername<-colnames(Y[k]) # new folder name
  
  dirname<-paste(wd,foldername,sep="") #directory name
  
  dir.create(dirname) # create new folder for individual Y-variable data processing
  
  setwd(dirname) # set working directory to new folder
  
  Ydata<-Y[,k] # Y variable

  Yindex<-Yindices[,k] # Non zero index

  
  Xsamplesincluded<-Samples[,Yindex] # X samples to include in regression
  
  
  Yincluded<-Ydata[Yindex] # Y variables to include in regression  
  
  Pcor<-cor(t(Xsamplesincluded),Yincluded,method=c("pearson")) #X-Y correlation
  
  if (any(is.na(Pcor)==TRUE)==FALSE){
  # Pcor probability function f test and t stat
 
  dfr<-ncol(Xsamplesincluded)-2 #degrees of freedom
  
  r2<-Pcor^2 #coefficient of determination
  
  Fstat<-r2*dfr/(1-r2) # F-stat
  
  Pcor_prob<-1-pf(Fstat,1,dfr) #p-value calculation
  
  Pcor.prob.adjusted<-p.adjust(Pcor_prob, method=p.adjust.methods, n=length(Pcor_prob)) ##multiple testing correction
  
  significant<-ifelse(Pcor.prob.adjusted<pvalue,Pcor,0) ##above significance threshold
  
  Pcor_results<-cbind(Pcor,Pcor_prob,Pcor.prob.adjusted,significant,ncol(Xsamplesincluded),foldername) 
  
  colnames(Pcor_results)<-c("Pcor","Pvalue",p.adjust.methods,"significant","Nsamples_included","Y_variable")
  
  Pcor_results<-cbind(XCMScolumns,Pcor_results)
  
  above_threshold<-which(significant>=Corr.thresh)
  
  above_threshold_dummy<-(significant>=Corr.thresh)*1
  
  Sum_above_threshold<-sum(significant>=Corr.thresh)
  
  Above_threshold_results<-Pcor_results[above_threshold,]
  
  Above_threshold_plots<-Xsamplesincluded[above_threshold,] ######subset for plotting linear regressions
    
  #########################################################################################################################################################
  ###subset the upper and lower classes of samples for box and whisker plot creation####
  
  if(sum(Yindices[,k])>=round((nrow(Y)*(box.prop)),digits=0)){
  high<-as.data.frame(Y[order(Y[,k],decreasing=TRUE)[1:round((nrow(Y)*(box.prop)),digits=0)],k])
  row.names(high)<-Yrownames[order(Y[,k],decreasing=TRUE)[1:round((nrow(Y)*(box.prop)),digits=0)]]
  colnames(high)<-"Box"
  low<-as.data.frame(Y[order(Y[,k],decreasing=FALSE)[1:round((nrow(Y)*(box.prop)),digits=0)],k])
  row.names(low)<-Yrownames[order(Y[,k],decreasing=FALSE)[1:round((nrow(Y)*(box.prop)),digits=0)]]
  colnames(low)<-"Box"
  } else {
    high<-as.data.frame(Y[Yindices[,k],k])
    row.names(high)<-Yrownames[Yindices[,k]]
    colnames(high)<-"Box"
    low<-as.data.frame(Y[order(Y[,k],decreasing=FALSE)[1:sum(Yindices[,k])],k])
    row.names(low)<-Yrownames[order(Y[,k],decreasing=FALSE)[1:sum(Yindices[,k])]]
    colnames(low)<-"Box"  
  }
  
  high[,1]<-TRUE
  low[,1]<-FALSE
  
  YBox<-rbind(high,low)
  YBox.row.names<-row.names(YBox)
  YBox<-YBox[order(YBox.row.names),]
  
  YBox<-ifelse(YBox==FALSE,"Y.low.zero","Y.high") ###change logical to category names
  
  #########################################################################################################################################################
    
  Above_threshold_Box_plots<-Samples[above_threshold,]
  
  #Quint.df<-t(Samples[,colnames(Samples) %in% food.table$X_DM])
  Above_threshold_Box_plots<-as.data.frame(Above_threshold_Box_plots[,colnames(Above_threshold_Box_plots) %in% YBox.row.names])
  
  Above_threshold_Box_plots<-Above_threshold_Box_plots[,order(colnames(Above_threshold_Box_plots))]
  
  Plot_names<-as.character(paste("M",round(Above_threshold_results[,2],digits=4),"T",round(Above_threshold_results[,3],digits=1),sep=""))
  

  #######Plotting above threshold correlation scatterplots#######
    
  if( Sum_above_threshold > 0) {
    
    message(paste("SAVING SCATTER AND BOX AND WHISKER PLOTS...",foldername,sep=""))#,quote=F)
    flush.console()
    
    pb<-txtProgressBar(min=0,max=nrow(Above_threshold_plots),style=3)#,width=300)#title="Scatter plot progress bar"
    
    
    for (i in 1:nrow(Above_threshold_plots)){ 
      
      ###progress bar for plots
      Sys.sleep(0.1)
      setTxtProgressBar(pb,i)
      flush.console()
     
      Titlename<-Plot_names[i]
      
     png(paste(Titlename,".",foldername,".png",sep=""),width=1200,height=1200,res=275)
            
     plot(log10(Yincluded),as.numeric(Above_threshold_plots[i,]), main=paste(foldername,Titlename),sub="(Log)",xlab=paste(foldername,Yunits), ylab="XObs_gLog_QC.LSC",xaxt="n",pch=19,cex=0.6)
         
      
     axis.labels<-round(lseq(min(Yincluded),max(Yincluded),10),digits=0)
     Log10Y<-log10(Yincluded)
     axis.points<-seq(min(Log10Y),max(Log10Y),length.out=10)
     axis(1, at=axis.points, labels = axis.labels)
      
      graphics.off()
      
      ####create box plots####
      boxplot.df<-data.frame(Box=YBox,XObs_gLog_QC.LSC=as.numeric(Above_threshold_Box_plots[i,]))
      
      ####T TEST HERE#####
      zero.boxplot<-boxplot.df[boxplot.df$Box=="Y.low.zero",2]
      nonzero.boxplot<-boxplot.df[boxplot.df$Box!="Y.low.zero",2]
      b1<-boxplot.stats(zero.boxplot)
      b2<-boxplot.stats(nonzero.boxplot)
      Above.zero.samples<-paste(round(sum.nonzero.above.zero<-sum(nonzero.boxplot>=b1$stats[[5]]),digits=0)," ","(",round(100*(sum.nonzero.above.zero/length(nonzero.boxplot)),digits=0),"%",")"," ",">1.5*IQR",sep="")
      Zero.outliers<-paste(zero.out<-length(which(zero.boxplot>b1$stats[[5]]))," ","(",round(100*(zero.out/length(zero.boxplot)),digits=0),"%",")"," ",">1.5*IQR",sep="") ###number of zero samples above 1.5*IQR threshold
      
      ###Box plot score calculation = percentage of non zero samples above zero group penalised by the proportion of samples above zero group whisker in zero group####
      box.plot.score<-((sum.nonzero.above.zero/length(nonzero.boxplot))-((zero.out/length(zero.boxplot))))
      box.plot.scores<-c(box.plot.scores,box.plot.score)
       
      g<-ggplot(boxplot.df,aes(x=Box,y=XObs_gLog_QC.LSC))+
        geom_boxplot(outlier.colour="red",outlier.size=0)+
        stat_summary(fun.y=mean, geom="point",shape=18, size=12,col="red",fill="red")+
        geom_jitter(position=position_jitter(w=0.15,h=0.15))+
        geom_hline(yintercept=b1$stats[[5]],colour="red",size=1)+
        annotate("text",label=Above.zero.samples,y=max(nonzero.boxplot)*1.1,x=1,size=3.5)+
        annotate("text",label=Zero.outliers,y=max(nonzero.boxplot)*1.1,x=2,size=3.5)+
        labs(title=Titlename)+
        theme_bw(20)
      
      ggsave(g,filename=paste(Titlename,".",foldername,"_BOXPLOT",".png",sep=""),width = 10, height = 10)
    }
    message("...Done")#,quote=F)
    flush.console() 
    
   }  
  
  ###Features above threshold aggregation#####
  
  Resultscolumns<-rbind(t(Pcor.prob.adjusted),t(Pcor),t(above_threshold_dummy))
  
  rownames(Resultscolumns)<-c(paste(foldername,"_pvalue_",pvalue,"_","MultTest_",p.adjust.methods,sep=""),paste(foldername,"_corr.coeff_nSamples_",ncol(Xsamplesincluded),sep=""),paste(foldername,"_Above_threshold",sep=""))
  
  Biomarker.number<-data.frame(cbind(ncol(Xsamplesincluded),foldername,Sum_above_threshold))
  
  Biomarker_numbers.df<-rbind(Biomarker_numbers.df,Biomarker.number)

  above_threshold_df<-rbind(above_threshold_df, t(above_threshold_dummy))
  
  results<-rbind(results,Resultscolumns)
  } 
}

results<-as.data.frame(t(results))
above_threshold_df<-as.data.frame(t(above_threshold_df))


Results_rowsums<-apply(above_threshold_df,1,sum) 

DummyMsignif<-as.data.frame(above_threshold_df[Results_rowsums>0,])

boxplot.score.df<-as.data.frame(matrix(0,ncol=ncol(above_threshold_df),nrow=nrow(DummyMsignif)))

boxplot.score.index<-DummyMsignif==1

###replace dummy matrix with box plot scores for each significant feature###

boxplot.score.df[boxplot.score.index]<-box.plot.scores

colnames(boxplot.score.df)<-paste(Biomarker_numbers.df[,2],rep(".boxplot.score",length.out=nrow(Biomarker_numbers.df)),sep="")

colnames(DummyMsignif)<-as.character(1:ncol(DummyMsignif))

significantmarker_data<-XCMScolumns[Results_rowsums>0,] #subset feature details above threshold

 #####Significant Feature HMDB weblink######

if(nrow(significantmarker_data)>2){
  
  significantmarker_mzMED<-significantmarker_data[,mzmed.column.name]
  HMDB.url<-data.frame()
  for (j in 1:length(significantmarker_mzMED)) {
    
    
    HMDB.url.link<-as.data.frame(paste("http://www.hmdb.ca/spectra/ms/search?utf8=%E2%9C%93&query_masses=",significantmarker_mzMED[j],"&tolerance=",HMDBtol,"&mode=",mode,"&commit=Search",sep=""))
    
    HMDB.url<-rbind(HMDB.url,HMDB.url.link)
    
  } 
  colnames(HMDB.url)<-"HMDB.url"
  significantmarker_data<-cbind(significantmarker_data,HMDB.url)  
}
colnames(significantmarker_data)[1]<-EIC.column.name

correlation.results<-results[Results_rowsums>0,]##sample data minus outliers for column binding following hmdb link generation

Ycorrelated.heatmap<-as.data.frame(correlation.results[,(as.logical(rep(c(0,1,0),length(Y))))])

sample.data<-Samples[Results_rowsums>0,]
##########################################################################################
##########################################################################################

significantmarker_data<-cbind(significantmarker_data,correlation.results,boxplot.score.df,sample.data)

setwd(wd)

####Heatmap######

if (heatmap==TRUE) {
  
  message("Hierarchical clustering and cluster ion identification...PLEASE WAIT")#,quote=F)
  flush.console()
  
  RowLabels<-as.character(significantmarker_data[,"name"])
 
  hmcols2<-rev(colorRampPalette(brewer.pal(10,"RdBu"))(256)) ##heatmap colour palette
  
  ####Create X-Y correlation heatmap#####
  if(ncol(Ycorrelated.heatmap)>1){

 #dist.m.method<- function(x) {
  #dist.co.x <- 1 - abs(x)
  #return(as.dist(dist.co.x))
#}

  pdf(paste("HM.Y_X",".pdf"))
  
  heatY<-heatmap.2(t(Ycorrelated.heatmap), dend="column",Colv=TRUE,Rowv=FALSE, symm=FALSE,scale="none",labRow=c(colnames(Ycorrelated.heatmap)),labCol=RowLabels,margins=c(7,7),col=hmcols2,trace="none",key=TRUE,keysize=1.5,cexCol=0.25,cexRow=0.4,hclustfun=function(x) hclust(x,method=hclust.method),distfun=function(x) dist(x,method=dist.m.method))#)
   
  dev.off()
  
} else if (ncol(Ycorrelated.heatmap)==1){
  message("WARNING: only one Y-variable for heatmap creation")
  flush.console()
}
  
  ###create X-X correlation matrix#####
  
  SignCor<-cor(t(sample.data),method=c("pearson"))
    
  pdf(paste("HM.X_X",".pdf"))
 
  ###create X-X correlation matrix#####
  
  heat<-heatmap.2(as.matrix(SignCor), dend="both",Colv=TRUE, symm=TRUE,scale="none",labRow=c(RowLabels),labCol=c(RowLabels),margins=c(2,2),col=hmcols2,trace="none",key=TRUE,cexCol=0.25,cexRow=0.25,hclustfun=function(x) hclust(x,method=hclust.method),distfun=function(x) dist(x,method=dist.m.method))#function(x) hclust(x,method=hclust.method)
  
  dev.off()
  
  
  ###Hierarchical clustering order####
  ReorderHierCl<-cbind(as.matrix(seq(1,ncol(SignCor),1)),as.matrix(heat$rowInd))
  ReorderHierCl<-ReorderHierCl[order(ReorderHierCl[,2]),]
  significantmarker_data<-cbind(ReorderHierCl,significantmarker_data)
  colnames(significantmarker_data)[1]<-"Hierclust.order"
  significantmarker_data<-significantmarker_data[order(significantmarker_data[,"Hierclust.order"]),]
  
  #####identify retention time clusters#####
  rtmed.cluster<-significantmarker_data[,RTmed.column.name]
  rtmed.cluster.shift<-c(rtmed.cluster[2],rtmed.cluster[-(length(rtmed.cluster))])
  rtmed.difference<-rtmed.cluster-rtmed.cluster.shift
  rtmed.cluster.mat<-data.frame(cbind(rtmed.cluster,rtmed.cluster.shift,rtmed.difference))
  
  ####match clusters based on similarity in retention time#########
  rtmed.cluster.mat$cluster.seq<-ifelse(x<-rtmed.cluster.mat$rtmed.difference<Clust.RT.tol & rtmed.cluster.mat$rtmed.difference>-Clust.RT.tol, cumsum(c(head(x, 1), tail(x, -1) - head(x, -1) == 1)), 0)
  rtmed.cluster.mat$cluster.seq<-c(ifelse((head(rtmed.cluster.mat$cluster.seq, -1) + tail(rtmed.cluster.mat$cluster.seq, -1) == tail(rtmed.cluster.mat$cluster.seq, -1)),tail(rtmed.cluster.mat$cluster.seq,-1),head(rtmed.cluster.mat$cluster,-1)),(tail(rtmed.cluster.mat$cluster.seq, 1)))
  significantmarker_data<-data.frame(cbind(rtmed.cluster.mat$cluster.seq,significantmarker_data))
  colnames(significantmarker_data)[1]<-"mz_clusters"
  
  ###order by mz then cluster####
  significantmarker_data<-significantmarker_data[order(-significantmarker_data[,"mz_clusters"],significantmarker_data[,mzmed.column.name]),]
 
  ###identify mass difference ####   
  significantmarker_mass.diff<-data.frame(c(0,tail(significantmarker_data$mzmed, -1) - head(significantmarker_data$mzmed, -1)))
  significantmarker_mass.diff[significantmarker_data[,"mz_clusters"]==0,]<-0 ###if no cluster replace with zero
  significantmarker_data<-cbind((seq(1,nrow(significantmarker_data),1)),significantmarker_data)
  colnames(significantmarker_data)[1]<-"cluster_ion_calc.order"
  
  ####identify isotopes from mass differences######
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff,2,function(x){ifelse(x<(0.984015583+0.003) & x>(0.984015583-0.003),1,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(0.99703+0.003) & x>(0.99703-0.003),2,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(1.00336+0.003) & x>(1.00336-0.003),3,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(1.007825032+0.003) & x>(1.007825032-0.003),4,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(1.979264556+0.003) & x>(1.979264556-0.003),5,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(1.9958+0.003) & x>(1.9958-0.003),6,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(2.00425+0.003) & x>(2.00425-0.003),7,x)}))
  

  ###remove isotopes before recalculation#####
  
  Isotope_mass.diff.index<-significantmarker_mass.diff.name %in% c(1,2,3,4,5,6,7)
  Isotopes<-significantmarker_data[Isotope_mass.diff.index==TRUE,]
  Isotope_mass.diff.name<-data.frame(significantmarker_mass.diff[Isotope_mass.diff.index==TRUE,])
  Isotope_name<-data.frame(significantmarker_mass.diff.name[Isotope_mass.diff.index==TRUE,]) 
  significantmarker_data<-significantmarker_data[Isotope_mass.diff.index==FALSE,]
  significantmarker_mass.diff<-data.frame(c(0,tail(significantmarker_data$mzmed, -1) - head(significantmarker_data$mzmed, -1)))
  
  ####replace hierarchical clusters with no retention time connection with a zero####
  significantmarker_mass.diff[significantmarker_data[,"mz_clusters"]==0,]<-0 ###if no cluster replace with zero
  
  ######recalculate potential isotopes#####
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff,2,function(x){ifelse(x<(0.984015583+0.003) & x>(0.984015583-0.003),1,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(0.99703+0.003) & x>(0.99703-0.003),2,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(1.00336+0.003) & x>(1.00336-0.003),3,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(1.007825032+0.003) & x>(1.007825032-0.003),4,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(1.979264556+0.003) & x>(1.979264556-0.003),5,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(1.9958+0.003) & x>(1.9958-0.003),6,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(2.00425+0.003) & x>(2.00425-0.003),7,x)}))
  
  ###remove isotopes again before recalculation of mass differences#####
  
  Isotope_mass.diff.index<-significantmarker_mass.diff.name %in% c(1,2,3,4,5,6,7)
  Isotopes<-rbind(Isotopes,significantmarker_data[Isotope_mass.diff.index==TRUE,])
  Isotope_mass.diff.name<-rbind(Isotope_mass.diff.name,data.frame(significantmarker_mass.diff[Isotope_mass.diff.index==TRUE,]))
  Isotope_name<-rbind(Isotope_name,data.frame(significantmarker_mass.diff.name[Isotope_mass.diff.index==TRUE,]))
  significantmarker_data<-significantmarker_data[Isotope_mass.diff.index==FALSE,]
  significantmarker_mass.diff<-data.frame(c(0,tail(significantmarker_data$mzmed, -1) - head(significantmarker_data$mzmed, -1)))
  
  
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff,2,function(x){ifelse(x<(a<-2.015650064)+((Clust.ppm/1000000)*a) & x>(b<-2.015650064)-((Clust.ppm/1000000)*b),8,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-2.999665647)+((Clust.ppm/1000000)*a) & x>(b<-2.999665647)-((Clust.ppm/1000000)*b),9,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-4.031300128)+((Clust.ppm/1000000)*a) & x>(b<-4.031300128)-((Clust.ppm/1000000)*b),10,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-12.0363855)+((Clust.ppm/1000000)*a) & x>(b<-12.0363855)-((Clust.ppm/1000000)*b),11,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-13.97926456)+((Clust.ppm/1000000)*a) & x>(b<-13.97926456)-((Clust.ppm/1000000)*b),12,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-14.01565006)+((Clust.ppm/1000000)*a) & x>(b<-14.01565006)-((Clust.ppm/1000000)*b),13,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-15.0234751)+((Clust.ppm/1000000)*a) & x>(b<-15.0234751)-((Clust.ppm/1000000)*b),14,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-15.01089904)+((Clust.ppm/1000000)*a) & x>(b<-15.01089904)-((Clust.ppm/1000000)*b),15,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-15.99491462)+((Clust.ppm/1000000)*a) & x>(b<-15.99491462)-((Clust.ppm/1000000)*b),16,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-17.0265491)+((Clust.ppm/1000000)*a) & x>(b<-17.0265491)-((Clust.ppm/1000000)*b),17,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-18.01056468)+((Clust.ppm/1000000)*a) & x>(b<-18.01056468)-((Clust.ppm/1000000)*b),18,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-21.981945)+((Clust.ppm/1000000)*a) & x>(b<-21.981945)-((Clust.ppm/1000000)*b),19,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-27.994915)+((Clust.ppm/1000000)*a) & x>(b<-27.994915)-((Clust.ppm/1000000)*b),20,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-28.006148)+((Clust.ppm/1000000)*a) & x>(b<-28.006148)-((Clust.ppm/1000000)*b),21,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-28.0313001)+((Clust.ppm/1000000)*a) & x>(b<-28.0313001)-((Clust.ppm/1000000)*b),22,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-28.9901636)+((Clust.ppm/1000000)*a) & x>(b<-28.9901636)-((Clust.ppm/1000000)*b),23,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-29.0027397)+((Clust.ppm/1000000)*a) & x>(b<-29.0027397)-((Clust.ppm/1000000)*b),24,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-29.0391252)+((Clust.ppm/1000000)*a) & x>(b<-29.0391252)-((Clust.ppm/1000000)*b),25,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-29.97417918)+((Clust.ppm/1000000)*a) & x>(b<-29.97417918)-((Clust.ppm/1000000)*b),26,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-30.01056468)+((Clust.ppm/1000000)*a) & x>(b<-30.01056468)-((Clust.ppm/1000000)*b),27,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-31.018498)+((Clust.ppm/1000000)*a) & x>(b<-31.018498)-((Clust.ppm/1000000)*b),28,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-31.98982924)+((Clust.ppm/1000000)*a) & x>(b<-31.98982924)-((Clust.ppm/1000000)*b),29,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-32.02621475)+((Clust.ppm/1000000)*a) & x>(b<-32.02621475)-((Clust.ppm/1000000)*b),30,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-34.0054793)+((Clust.ppm/1000000)*a) & x>(b<-34.0054793)-((Clust.ppm/1000000)*b),31,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-34.0530982)+((Clust.ppm/1000000)*a) & x>(b<-34.0530982)-((Clust.ppm/1000000)*b),32,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-37.955882)+((Clust.ppm/1000000)*a) & x>(b<-37.955882)-((Clust.ppm/1000000)*b),33,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-41.02669)+((Clust.ppm/1000000)*a) & x>(b<-41.02669)-((Clust.ppm/1000000)*b),34,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-42.010565)+((Clust.ppm/1000000)*a) & x>(b<-42.010565)-((Clust.ppm/1000000)*b),35,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-43.0058137)+((Clust.ppm/1000000)*a) & x>(b<-43.0058137)-((Clust.ppm/1000000)*b),36,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-43.0183897)+((Clust.ppm/1000000)*a) & x>(b<-43.0183897)-((Clust.ppm/1000000)*b),37,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-43.0547752)+((Clust.ppm/1000000)*a) & x>(b<-43.0547752)-((Clust.ppm/1000000)*b),38,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-43.9898292)+((Clust.ppm/1000000)*a) & x>(b<-43.9898292)-((Clust.ppm/1000000)*b),39,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-44.9977)+((Clust.ppm/1000000)*a) & x>(b<-44.9977)-((Clust.ppm/1000000)*b),40,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-46.0419)+((Clust.ppm/1000000)*a) & x>(b<-46.0419)-((Clust.ppm/1000000)*b),41,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-47.98474386)+((Clust.ppm/1000000)*a) & x>(b<-47.98474386)-((Clust.ppm/1000000)*b),42,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-56.0626003)+((Clust.ppm/1000000)*a) & x>(b<-56.0626003)-((Clust.ppm/1000000)*b),43,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-57.02146372)+((Clust.ppm/1000000)*a) & x>(b<-57.02146372)-((Clust.ppm/1000000)*b),44,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-57.0704253)+((Clust.ppm/1000000)*a) & x>(b<-57.0704253)-((Clust.ppm/1000000)*b),45,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-58.0530982)+((Clust.ppm/1000000)*a) & x>(b<-58.0530982)-((Clust.ppm/1000000)*b),46,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-59.0133043)+((Clust.ppm/1000000)*a) & x>(b<-59.0133043)-((Clust.ppm/1000000)*b),47,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-59.0371138)+((Clust.ppm/1000000)*a) & x>(b<-59.0371138)-((Clust.ppm/1000000)*b),48,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-60.0211294)+((Clust.ppm/1000000)*a) & x>(b<-60.0211294)-((Clust.ppm/1000000)*b),49,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-63.961904)+((Clust.ppm/1000000)*a) & x>(b<-63.961904)-((Clust.ppm/1000000)*b),50,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-75.032029)+((Clust.ppm/1000000)*a) & x>(b<-75.032029)-((Clust.ppm/1000000)*b),51,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-78.95850549)+((Clust.ppm/1000000)*a) & x>(b<-78.95850549)-((Clust.ppm/1000000)*b),52,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-79.956819)+((Clust.ppm/1000000)*a) & x>(b<-79.956819)-((Clust.ppm/1000000)*b),53,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-82.0530982)+((Clust.ppm/1000000)*a) & x>(b<-82.0530982)-((Clust.ppm/1000000)*b),54,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-86.000395)+((Clust.ppm/1000000)*a) & x>(b<-86.000395)-((Clust.ppm/1000000)*b),55,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-89.047679)+((Clust.ppm/1000000)*a) & x>(b<-89.047679)-((Clust.ppm/1000000)*b),56,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-97.96737954)+((Clust.ppm/1000000)*a) & x>(b<-97.96737954)-((Clust.ppm/1000000)*b),57,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-103.0091848)+((Clust.ppm/1000000)*a) & x>(b<-103.0091848)-((Clust.ppm/1000000)*b),58,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-107.0040994)+((Clust.ppm/1000000)*a) & x>(b<-107.0040994)-((Clust.ppm/1000000)*b),59,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-119.0041)+((Clust.ppm/1000000)*a) & x>(b<-119.0041)-((Clust.ppm/1000000)*b),60,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-121.019753)+((Clust.ppm/1000000)*a) & x>(b<-121.019753)-((Clust.ppm/1000000)*b),61,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-129.042594)+((Clust.ppm/1000000)*a) & x>(b<-129.042594)-((Clust.ppm/1000000)*b),62,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-146.036779)+((Clust.ppm/1000000)*a) & x>(b<-146.036779)-((Clust.ppm/1000000)*b),63,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-146.069142)+((Clust.ppm/1000000)*a) & x>(b<-146.069142)-((Clust.ppm/1000000)*b),64,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-152.010959)+((Clust.ppm/1000000)*a) & x>(b<-152.010959)-((Clust.ppm/1000000)*b),65,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-161.014668)+((Clust.ppm/1000000)*a) & x>(b<-161.014668)-((Clust.ppm/1000000)*b),66,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-162.052825)+((Clust.ppm/1000000)*a) & x>(b<-162.052825)-((Clust.ppm/1000000)*b),67,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-163.030318)+((Clust.ppm/1000000)*a) & x>(b<-163.030318)-((Clust.ppm/1000000)*b),68,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-176.032088)+((Clust.ppm/1000000)*a) & x>(b<-176.032088)-((Clust.ppm/1000000)*b),69,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-176.047344)+((Clust.ppm/1000000)*a) & x>(b<-176.047344)-((Clust.ppm/1000000)*b),70,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-178.041213)+((Clust.ppm/1000000)*a) & x>(b<-178.041213)-((Clust.ppm/1000000)*b),71,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-180.06339)+((Clust.ppm/1000000)*a) & x>(b<-180.06339)-((Clust.ppm/1000000)*b),72,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-192.027)+((Clust.ppm/1000000)*a) & x>(b<-192.027)-((Clust.ppm/1000000)*b),73,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-194.042655)+((Clust.ppm/1000000)*a) & x>(b<-194.042655)-((Clust.ppm/1000000)*b),74,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-198.014035)+((Clust.ppm/1000000)*a) & x>(b<-198.014035)-((Clust.ppm/1000000)*b),75,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-203.079373)+((Clust.ppm/1000000)*a) & x>(b<-203.079373)-((Clust.ppm/1000000)*b),76,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-215.994605)+((Clust.ppm/1000000)*a) & x>(b<-215.994605)-((Clust.ppm/1000000)*b),77,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-221.089937)+((Clust.ppm/1000000)*a) & x>(b<-221.089937)-((Clust.ppm/1000000)*b),78,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-239.993994)+((Clust.ppm/1000000)*a) & x>(b<-239.993994)-((Clust.ppm/1000000)*b),79,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-255.988909)+((Clust.ppm/1000000)*a) & x>(b<-255.988909)-((Clust.ppm/1000000)*b),80,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-258.004564)+((Clust.ppm/1000000)*a) & x>(b<-258.004564)-((Clust.ppm/1000000)*b),81,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-273.096087)+((Clust.ppm/1000000)*a) & x>(b<-273.096087)-((Clust.ppm/1000000)*b),82,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-273.999479)+((Clust.ppm/1000000)*a) & x>(b<-273.999479)-((Clust.ppm/1000000)*b),83,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-275.111737)+((Clust.ppm/1000000)*a) & x>(b<-275.111737)-((Clust.ppm/1000000)*b),84,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-305.068161)+((Clust.ppm/1000000)*a) & x>(b<-305.068161)-((Clust.ppm/1000000)*b),85,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-307.083811)+((Clust.ppm/1000000)*a) & x>(b<-307.083811)-((Clust.ppm/1000000)*b),86,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-352.06418)+((Clust.ppm/1000000)*a) & x>(b<-352.06418)-((Clust.ppm/1000000)*b),87,x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x<(a<-388.08531)+((Clust.ppm/1000000)*a) & x>(b<-388.08531)-((Clust.ppm/1000000)*b),88,x)}))
    
  #####identify unlabelled cluster differences######
  mass.diff.seq<-seq(1,88,1)
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x %in% mass.diff.seq,x,NA)}))
  
  ####Name Isotopes######
  if(nrow(Isotope_name)>0)
  {
  Isotope_name<-(apply(Isotope_name,2,function(x){ifelse(x==1,"[+O-NH3]",x)}))
  Isotope_name<-(apply(Isotope_name,2,function(x){ifelse(x==2,"[N15 isotope]",x)}))
  Isotope_name<-(apply(Isotope_name,2,function(x){ifelse(x==3,"[C13 isotope]",x)}))
  Isotope_name<-(apply(Isotope_name,2,function(x){ifelse(x==4,"[+H]",x)}))
  Isotope_name<-(apply(Isotope_name,2,function(x){ifelse(x==5,"[+O -CH2]",x)}))
  Isotope_name<-(apply(Isotope_name,2,function(x){ifelse(x==6,"[S34 isotope]",x)}))
  Isotope_name<-(apply(Isotope_name,2,function(x){ifelse(x==7,"[O18 isotope]",x)}))
  }
  ###Name fragments##### 
  
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==8,"[+H2]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==9,"[+OH-N]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==10,"[+H4]/[-H4]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==11,"[+O-C2H4]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==12,"[+O -H2]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==13,"[+CH2]/[-CH2]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==14,"[-CH3]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==15,"[+O2 -NH3]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==16,"[Hydroxylation/+O]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==17,"[M+NH4]/[-NH3]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==18,"[+H2O]/[-H2O]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==19,"[M+Na]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==20,"[-CO]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==21,"[-N2]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==22,"[-C2H4] ",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==23,"[+H-NO] ",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==24,"[-CHO]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==25,"[-C2H5] ",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==26,"[+O2-H2]/[+H2-O2]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==27,"[+OCH2]/[-CH2O] ",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==28,"[M+Methanol CH4O]/[-OCH3]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==29,"[2 x Hydroxylation]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==30,"[M+MeOH+H]/-[MeOH]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==31,"[+2 OH]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==32,"[M+NH3.NH4]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==33,"[M+K]/-[K]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==34,"[M+MeCN+H]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==35,"[Acetylation shift]/-[AcKetene]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==36,"[-HCNO]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==37,"[-CH3CO]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==38,"[-C3H7]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==39,"[-CO2]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==40,"[M+HCOOH]/[-HCOOH/+H-NO2]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==41,"[-EtOH]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==42,"[+O3]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==43,"[-C4H8]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==44,"[Glycyl conjugation shift]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==45,"[-C4H9]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==46,"[M+MeCN+NH4]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==47,"[-acetate]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==48,"[-Acetamide]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==49,"[-CH3COOH]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==50,"[-SO2]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==51,"[-Gly]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==52,"[phosphate - PO3 conjugation]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==53,"[SO3 conjugation shift]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==54,"[M+2ACN+H]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==55,"[-Malonyl]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==56,"[-Cys conj-ala]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==57,"[H2SO4- Sulfate conjugation/phosphate conjugation]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==58,"[Cysteinyl conjugation]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==59,"[Taurine conjugation]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==60,"[S-Cysteine conjugation shift]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==61,"[-Cys]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==62,"[-GSH AnhydroGlu]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==63,"[-Coumaroyl loss]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==64,"[-Ala-Gly/ GSH Glu loss]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==65,"[-Galloyl loss]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==66,"[N-AcetylCysteine conjugation shift]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==67,"[Glucose shift]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==68,"[-N-AcCys]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==69,"[AnhydroGlucuronide conjugation shift]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==70,"[-Feruloyl loss]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==71,"[-Cys-Gly loss]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==72,"[-Gluc loss]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==73,"[-Hydroxylation+Gluc]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==74,"[-Gluc]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==75,"[-AnhydroGluc + Na]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==76,"[-AnhydroGlucNAc]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==77,"[-Gluc + Na]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==78,"[-GlcNAc loss]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==79,"[-AnhydroGluc + SO2]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==80,"[-AnhydroGluc + SO3]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==81,"[-Gluc + SO2]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==82,"[-GSH GluAlaGly-2H]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==83,"[-Gluc + SO3]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==84,"[-GSH-GluAlaGly]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==85,"[S-Glutathione conjugation shift]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==86,"[-GSH]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==87,"[-DiAnhydroGlu]",x)}))
  significantmarker_mass.diff.name<-(apply(significantmarker_mass.diff.name,2,function(x){ifelse(x==88,"[-DiGlu]",x)}))
  
  ###rebind isotopes with potential in-source fragments#####
  significantmarker_data<-data.frame(cbind((seq(1,nrow(significantmarker_data),1)),significantmarker_mass.diff.name,significantmarker_mass.diff,significantmarker_data))
  Isotopes<-data.frame(cbind((rep(NA,nrow(Isotopes))),Isotope_name,Isotope_mass.diff.name,Isotopes))
  colnames(Isotopes)<-colnames(significantmarker_data)
  significantmarker_data<-as.data.frame(rbind(significantmarker_data,Isotopes))
  
  colnames(significantmarker_data)[1]<-"Isotope_rem.cluster_calc.order"
  colnames(significantmarker_data)[2]<-"Potential_cluster_ions"
  colnames(significantmarker_data)[3]<-"cluster_ions_mz.diff"
  
message("...Done")#,quote=F)
flush.console()

}
significantmarker_data<-significantmarker_data[order(significantmarker_data[,EIC.column.name]),]
###reorder according to mzmed rtmed and ID columns
mz.rt.cols<-which(colnames(significantmarker_data) %in% c(EIC.column.name,mzmed.column.name,RTmed.column.name))
data.cols<-c(1:ncol(significantmarker_data))[-mz.rt.cols]
significantmarker_data<-significantmarker_data[,c(1,mz.rt.cols,data.cols[c(2:length(data.cols))])]

#significantmarker_data<-subset(significantmarker_data,select=-c(X2))
Ycor<-Y

####split significant features into their respective subfolders and save also add column of scatterplot locations#####

for (k in 1:ncol(Ycor)){
  
  foldername<-colnames(Ycor[k]) # new folder name
 
  dirname<-paste(wd,foldername,sep="") #directory name
  
  setwd(dirname) # set working directory to new folder
  
  Yresultsindex<-as.logical(significantmarker_data[,paste(foldername,"_Above_threshold",sep="")])
  
  Yresults.signif<-significantmarker_data[Yresultsindex,]
  
  Above_thresh_csvlabel<-paste(foldername,"_results_above_threshold",".csv",sep="")
  
  if( nrow(Yresults.signif)>1) {
    
    Plot.url.names<-as.character(Yresults.signif[,"name"])
    
    Scatter.url<-data.frame()
    
    for (k in 1:nrow(Yresults.signif)){ 
      
      Plot.url.Feat.name<-Plot.url.names[k]
      
      Scatter.plot.url<-as.data.frame(paste(wd,foldername,"/",Plot.url.Feat.name,foldername,".png",sep=""))
      
      Scatter.url<-rbind(Scatter.url,Scatter.plot.url)
    }
    
    colnames(Scatter.url)<-"Scatter.plot.file.location"
    Yresults.signif<-cbind(Yresults.signif,Scatter.url)
    
    write.csv(Yresults.signif,Above_thresh_csvlabel,row.names=FALSE)
  } else if (nrow(Yresults.signif)<=1){
    
    Ynocor<-Ycor[,foldername]
  
  setwd(wd) 
######save new Y matrix as .csv if Y variables have been removed#####
    
  if (length(ncol(Ynocor)!=ncol(Ycor))!=0) {
    write.csv(Ynocor,"Y_uncorrelated_removed.csv")
  }  
}
}

setwd(wd)

write.csv(significantmarker_data,"Features_Above_threshold.csv",row.names=FALSE)

} else if (ncol(Y)==0) { 
  print("Error: insufficient Y variable values supplied, decrease non-zero parameter or try others")
  
}
}

guiv(Auto.MV.Regress,
     argText=list(X="Auto.PCA output file name ? ",
                  Yvar="Y variable (.csv) file name following outlier removal (e.g. Y.outliers.removed.csv) ",                
                  non.zero="Percentage of non-zero values to retain a Y-variable (else not considered) ? ",
                  Corr.thresh="Pearson Correlation coefficient threshold ? ",
                  pvalue="maximum p-value threshold ? ",
                  p.adjust.methods="multiple testing correction method ? ",
                  heatmap="Inter-feature HCA to identify highly correlated features/cluster ions ? ",
                  hclust.method="Hierarchical Clustering method ? ",
                  dist.m.method="Distance (dissimilarity) calculation method ? ",
                  Clust.RT.tol="Cluster ion retention time tolerance (in seconds) ", 
                  Clust.ppm="Cluster ion mass accuracy tolerance (in ppm) ",
                  mode="Which polarity was the data acquired in ? ",
                  box.prop="Quantile for high/low box plot generation (default is quintiles i.e. 0.2)",
                  Yunits="Y units to be displayed on the plots ? ",
                  HMDBtol="delta mass accuracy tolerance for HMDB hyperlink generation ? "),
     argSlider=list(non.zero=c(0,100,0.5),Corr.thresh=c(0,1,0.01),pvalue=c(1,0,0.0001),Clust.RT.tol=c(0,10,0.5),Clust.ppm=c(0,40,0.5),HMDB.tol=c(0,5,0.0001)),
     argOption=list(p.adjust.methods=c("NONE","BH","FDR","BONFERRONI","HOLM","HOCHBERG"),heatmap=c("TRUE","FALSE"),hclust.method=c("WARD.D","WARD.D2","SINGLE","COMPLETE","AVERAGE","MCQUITTY"),dist.m.method=c("EUCLIDEAN", "MAXIMUM", "MANHATTAN", "CANBERRA", "BINARY", "MINKOWSKI"),mode=c("positive","negative")),helps=NULL)

###END###
