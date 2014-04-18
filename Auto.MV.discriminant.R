Auto.MV.discriminant<-function(X="PCA.outliers.removed.MS.MS.spectra.csv",X.raw="Blank_15_005_pvalue_filtered.csv",Yvar="All_Y_variables_Dummy_Whole_Quint.csv",p.adjust.method="BH",pvalue=0.05,box.prop=1/5,Fold.Ch=2,
wd="Z:\\Raw_data\\EPIC_POS_Ur_n481_04_02_14\\",heatmap=TRUE, hclust.method="complete", dist.m.method="euclidean", Yunits="Y_units", HMDBtol=0.005, Clust.RT.tol=2, mode="positive", Clust.ppm=10){
  
  require(pROC)
  require(ggplot2)
  
  
XCMS.dir<-paste(wd,"XCMS\\",sep="")

setwd(XCMS.dir)

##automatically identify file extension##
###load data for fold change calculation###

if(substr(X.raw,nchar(X.raw)-2,nchar(X.raw))=="tsv"){
  X.raw<-read.table(X.raw,sep="\t",header=T)
} else if (substr(X.raw,nchar(X.raw)-2,nchar(X))=="txt"){
  X.raw<-read.table(X.raw,sep="\t",header=T)
} else if (substr(X.raw,nchar(X.raw)-2,nchar(X.raw))=="csv"){
  X.raw<-read.csv(X.raw,header=T)
}



wd<-paste(wd,"Auto.PCA.results\\",sep="")

setwd(wd)

Samples<-read.csv(X,header=TRUE)

Y<-read.csv(Yvar,header=TRUE,row.names=1)

###match feature names raw data and PCA outliers removed###

X.raw<-X.raw[X.raw$name %in% Samples$name,]

###create subdirectory###

wd<-paste(substr(wd,1,nchar(wd)-17),"Auto.MV.Discriminant.results","\\",sep="")

dir.create(wd)

setwd(wd)

XCMScolumnsIndex<-c(1:which(colnames(Samples)=="RSD_corr_below"))

XCMScolumns<-Samples[,XCMScolumnsIndex]

Samples<-Samples[,-XCMScolumnsIndex]

###match sample names raw data and PCA outliers removed###

X.raw<-X.raw[,colnames(X.raw) %in% colnames(Samples)]

###automatically decide if Y variables supplied are continuous or binary for### 
###discriminant analysis###


if (max(is.na(Y[,1]))==1){

  ###save all parameters used in a dated .csv file for future reference###
  Parameters<-data.frame(X,Yvar,p.adjust.methods,
                         hclust.method,dist.m.method,pvalue.cutoff=pvalue,min.Fold.change=Fold.Ch,
                         Cluster.RT.tol=Clust.RT.tol,Mass.Accuracy.ClusterID=Clust.ppm)
  
  date<-Sys.time()
  date<-gsub("-",".",date)
  write.csv(Parameters,paste("Parameters",substr(date,1,10),".csv",sep=" "),row.names=FALSE)
  
  
  Y.X.combined<-cbind(Y,t(Samples))
  
  Y.X.raw<-cbind(Y,t(X.raw))

  t.test.results<-data.frame()
  
  ROC.box.results<-data.frame()
  
  for (k in 1:ncol(Y)){
    
    Yvar.name<-colnames(Y[k])
    
    dirname<-paste(wd,Yvar.name,sep="") #directory name
    
    dir.create(dirname) # create new folder for individual Y-variable data processing
    
    setwd(dirname) # set working directory to new folder
    
    t.test.df<-Y.X.combined[which(Y.X.combined[,k]>=0),]    
    
    ###fold change calculation###
    
    fold.change<-apply(Y.X.raw[which(Y[,Yvar.name]==1),-c(1:ncol(Y))],2,mean)/
      apply(Y.X.raw[which(Y[,Yvar.name]==0),-c(1:ncol(Y))],2,mean)
    
    ####apply t test and create data frame####
    
    t.test.result<-data.frame(t(sapply(t.test.df[-(1:ncol(Y))], function(x) unlist(t.test(x~t.test.df[,Yvar.name])[c("p.value","statistic")]))))
    
    ###p.value multiple testing adjustment
    Adjust.p.value<-p.adjust(a<-t.test.result[,1],method=p.adjust.method, n=nrow(t.test.result)) 
        
    ###number below p value threshold and fold change###
    
    t.test.indx<-(Adjust.p.value<=pvalue & fold.change > Fold.Ch)*1
    
    t.test.result<-data.frame(t.test.result,Adjust.p.value,fold.change,t.test.indx)
    
    colnames(t.test.result)<-c(paste(Yvar.name,"p.value",sep="."),paste(Yvar.name,"t.stat",sep="."),
                               paste("p.adjusted.",Yvar.name,sep=""),paste(Yvar.name,"fold.Change",sep="."),
                               paste(Yvar.name,"Signif.thresh",sep="."))
     
    zero.boxplot<-Y.X.combined[which(Y.X.combined[,k]==0),-c(1:ncol(Y))]
    nonzero.boxplot<-Y.X.combined[which(Y.X.combined[,k]==1),-c(1:ncol(Y))]
    colnames(zero.boxplot)<-XCMScolumns$name
    colnames(nonzero.boxplot)<-XCMScolumns$name
    
    Samples.comb<-data.frame(t(zero.boxplot),t(nonzero.boxplot))
    
    ###BOX PLOT, MANHATTAN AND AREA UNDER THE CURVE###
    if (sum(t.test.indx)>1){
      zero.boxplot<-zero.boxplot[,t.test.indx==1]
      nonzero.boxplot<-nonzero.boxplot[,t.test.indx==1]
      Plot_names<-paste(XCMScolumns$XCMS_EIC,XCMScolumns$name,sep="_")[t.test.indx==1]
      box.plot.comb<-rbind(zero.boxplot,nonzero.boxplot)
      YBox<-c(rep("Y.Class.0",length.out=nrow(zero.boxplot)),rep("Y.Class.1",length.out=nrow(nonzero.boxplot)))
      
      ###empty vectors for box plot score and ROC results###
      box.plot.scores<-numeric()
      ROC.results<-data.frame()
      
      for (j in 1:ncol(zero.boxplot)){
        Titlename<-Plot_names[j]
        boxplot.df<-data.frame(Box=YBox,XObs_gLog_QC.LSC=as.numeric(box.plot.comb[,j]))  
        b1<-boxplot.stats(zero.boxplot[,j])
        b2<-boxplot.stats(nonzero.boxplot[,j])
        Above.zero.samples<-paste(round(sum.nonzero.above.zero<-sum(nonzero.boxplot[,j]>=b1$stats[[5]]),digits=0)," ","(",round(100*(sum.nonzero.above.zero/nrow(nonzero.boxplot)),digits=0),"%",")"," ",">1.5*IQR",sep="")
        Zero.outliers<-paste(zero.out<-length(which(zero.boxplot[,j]>b1$stats[[5]]))," ","(",round(100*(zero.out/nrow(zero.boxplot)),digits=0),"%",")"," ",">1.5*IQR",sep="") ###number of zero samples above 1.5*IQR threshold
        
        ###Box plot score calculation = percentage of non zero samples above zero group penalised by the proportion of samples above zero group whisker in zero group####
        box.plot.score<-((sum.nonzero.above.zero/length(nonzero.boxplot))-((zero.out/length(zero.boxplot))))
        box.plot.scores<-c(box.plot.scores,box.plot.score)
        
        ###ROC curve###
        png(paste(Titlename,Yvar.name,"ROC_curve",".png",sep="_"))#,width=1200,height=1200,res=275)
        ROC <- plot.roc(YBox,as.numeric(box.plot.comb[,j]),Box.plot.comb, main=paste(Titlename,sep=""), percent=TRUE, ci=TRUE, print.auc=TRUE) 
        dev.off()
        
        ci<-ROC$ci
        ROC.results<-rbind(ROC.results,ci)
        
        ####BOX and WHISKER plot#####
        g<-ggplot(boxplot.df,aes(x=Box,y=XObs_gLog_QC.LSC))+
          geom_boxplot(outlier.colour="red",outlier.size=0)+
          stat_summary(fun.y=mean, geom="point",shape=18, size=12,col="red",fill="red")+
          geom_jitter(position=position_jitter(w=0.15,h=0.15))+
          geom_hline(yintercept=b1$stats[[5]],colour="red",size=1)+
          annotate("text",label=Above.zero.samples,y=max(nonzero.boxplot[,j])*1.1,x=2,size=3.5)+
          annotate("text",label=Zero.outliers,y=max(nonzero.boxplot[,j])*1.1,x=1,size=3.5)+
          labs(title=Titlename)+
          theme_bw(20)
        
        ggsave(g,filename=paste(Titlename,".",Yvar.name,"_BOXPLOT",".png",sep=""))
      }
      ###empty data.frame for storage of ROC and box plot results###
      ROC.box.result<-matrix(0,nrow=nrow(t.test.result),ncol=4)
      ROC.box.result[t.test.indx==1,1]<-box.plot.scores
      ROC.box.result[t.test.indx==1,2]<-ROC.results[,1]
      ROC.box.result[t.test.indx==1,3]<-ROC.results[,2]
      ROC.box.result[t.test.indx==1,4]<-ROC.results[,3]
      
      colnames(ROC.box.result)<-c(paste(Yvar.name,".box.plot.score",sep=""),paste(Yvar.name,".lower.CI",sep=""),paste(Yvar.name,".AUC",sep=""),paste(Yvar.name,".higher.CI",sep=""))
     
      ###create subdirectory variable results###
      
      Variable.results<-cbind(XCMScolumns,t.test.result,ROC.box.result,Samples.comb)
            
      Variable.results<-Variable.results[which(t.test.indx==1),]
      
      write.csv(Variable.results,paste(Yvar.name,"results.csv",sep="_"),row.names=F)
      
      t.test.results<-rbind(t.test.results,t(t.test.result))
      
      ROC.box.results<-rbind(ROC.box.results,t(ROC.box.result))
      
      
  }
    
  }   
    
############continuous variables conditional####################  
} else if (max(Y[,1]>1)){
  ###save all parameters used in a dated .csv file for future reference###
  Parameters<-data.frame(X,Yvar,p.adjust.methods,box.plot.proportions=box.prop,
                         hclust.method,dist.m.method,pvalue.cutoff=pvalue,min.Fold.change=Fold.Ch,
                         Cluster.RT.tol=Clust.RT.tol,Mass.Accuracy.ClusterID=Clust.ppm)
  
  date<-Sys.time()
  date<-gsub("-",".",date)
  write.csv(Parameters,paste("Parameters",substr(date,1,10),".csv",sep=" "),row.names=FALSE)
    
  Y.X.combined<-cbind(Y,t(Samples))
  
  Y.X.raw<-cbind(Y,t(X.raw))
  
  t.test.results<-data.frame()
    
  ROC.box.results<-data.frame()
  
  for (k in 1:ncol(Y)){
    
    Yvar.name<-colnames(Y[k])
    
    dirname<-paste(wd,Yvar.name,sep="") #directory name
    
    dir.create(dirname) # create new folder for individual Y-variable data processing
    
    setwd(dirname) # set working directory to new folder
 
    high<-Y.X.combined[order(Y[,Yvar.name],decreasing=TRUE)[1:round(nrow(Y)*box.prop,digits=0)],]
    
    zero<-Y.X.combined[order(Y[,Yvar.name],decreasing=FALSE)[1:round(nrow(Y)*box.prop,digits=0)],]
    
    fold.change<-apply(Y.X.raw[order(Y[,Yvar.name],decreasing=TRUE)[1:round(nrow(Y)*box.prop,digits=0)],-c(1:ncol(Y))],2,mean)/
                 apply(Y.X.raw[order(Y[,Yvar.name],decreasing=FALSE)[1:round(nrow(Y)*box.prop,digits=0)],-c(1:ncol(Y))],2,mean)
      
      
    high[,k]<-1
    zero[,k]<-0
    t.test.df<-rbind(high,zero)
      
    ####apply t test and create data frame####
    
    t.test.result<-data.frame(t(sapply(t.test.df[-(1:ncol(Y))], function(x) unlist(t.test(x~t.test.df[,Yvar.name])[c("p.value","statistic")]))))
    
    ###p.value multiple testing adjustment
    Adjust.p.value<-p.adjust(a<-t.test.result[,1],method=p.adjust.method, n=nrow(t.test.result)) 
    
    ###number below p value threshold and fold change###
    
    t.test.indx<-(Adjust.p.value<=pvalue & fold.change > Fold.Ch)*1
    
    t.test.result<-data.frame(t.test.result,Adjust.p.value,fold.change,t.test.indx)
    
    colnames(t.test.result)<-c(paste(Yvar.name,"p.value",sep="."),paste(Yvar.name,"t.stat",sep="."),
                               paste("p.adjusted.",Yvar.name,sep=""),paste(Yvar.name,"fold.Change",sep="."),
                               paste(Yvar.name,"Signif.thresh",sep="."))
   
    zero.boxplot<-zero[,-c(1:ncol(Y))]
    nonzero.boxplot<-high[,-c(1:ncol(Y))]
    
    colnames(zero.boxplot)<-XCMScolumns$name
    colnames(nonzero.boxplot)<-XCMScolumns$name
    
    Samples.comb<-data.frame(t(zero.boxplot),t(nonzero.boxplot))
    
    ###BOX PLOT, MANHATTAN AND AREA UNDER THE CURVE###
   if (sum(t.test.indx)>1){
     zero.boxplot<-zero.boxplot[,t.test.indx==1]
     nonzero.boxplot<-nonzero.boxplot[,t.test.indx==1]
     Plot_names<-paste(XCMScolumns$XCMS_EIC,XCMScolumns$name,sep="_")[t.test.indx==1]
    box.plot.comb<-rbind(zero.boxplot,nonzero.boxplot)
     YBox<-c(rep("Y.low",length.out=nrow(zero.boxplot)),rep("Y.high",length.out=nrow(nonzero.boxplot)))
     
     ###empty vectors for box plot score and ROC results###
     box.plot.scores<-numeric()
     ROC.results<-data.frame()
     
     for (j in 1:ncol(zero.boxplot)){
    Titlename<-Plot_names[j]
    boxplot.df<-data.frame(Box=YBox,XObs_gLog_QC.LSC=as.numeric(box.plot.comb[,j]))  
    b1<-boxplot.stats(zero.boxplot[,j])
    b2<-boxplot.stats(nonzero.boxplot[,j])
    Above.zero.samples<-paste(round(sum.nonzero.above.zero<-sum(nonzero.boxplot[,j]>=b1$stats[[5]]),digits=0)," ","(",round(100*(sum.nonzero.above.zero/nrow(nonzero.boxplot)),digits=0),"%",")"," ",">1.5*IQR",sep="")
    Zero.outliers<-paste(zero.out<-length(which(zero.boxplot[,j]>b1$stats[[5]]))," ","(",round(100*(zero.out/nrow(zero.boxplot)),digits=0),"%",")"," ",">1.5*IQR",sep="") ###number of zero samples above 1.5*IQR threshold
    
    ###Box plot score calculation = percentage of non zero samples above zero group penalised by the proportion of samples above zero group whisker in zero group####
    box.plot.score<-((sum.nonzero.above.zero/length(nonzero.boxplot))-((zero.out/length(zero.boxplot))))
    box.plot.scores<-c(box.plot.scores,box.plot.score)
    
    ###ROC curve###
    png(paste(Titlename,Yvar.name,"ROC_curve",".png",sep="_"))#,width=1200,height=1200,res=275)
    ROC <- plot.roc(YBox,as.numeric(box.plot.comb[,j]),Box.plot.comb, main=paste(Titlename,sep=""), percent=TRUE, ci=TRUE, print.auc=TRUE) 
    dev.off()
    
    ci<-ROC$ci
    ROC.results<-rbind(ROC.results,ci)
    
    ####BOX and WHISKER plot#####
    g<-ggplot(boxplot.df,aes(x=Box,y=XObs_gLog_QC.LSC))+
      geom_boxplot(outlier.colour="red",outlier.size=0)+
      stat_summary(fun.y=mean, geom="point",shape=18, size=12,col="red",fill="red")+
      geom_jitter(position=position_jitter(w=0.15,h=0.15))+
      geom_hline(yintercept=b1$stats[[5]],colour="red",size=1)+
      annotate("text",label=Above.zero.samples,y=max(nonzero.boxplot[,j])*1.1,x=1,size=3.5)+
      annotate("text",label=Zero.outliers,y=max(nonzero.boxplot[,j])*1.1,x=2,size=3.5)+
      labs(title=Titlename)+
      theme_bw(20)
    
    ggsave(g,filename=paste(Titlename,".",Yvar.name,"_BOXPLOT",".png",sep=""))
      }
  ###empty data.frame for storage of ROC and box plot results###
  ROC.box.result<-matrix(0,nrow=nrow(t.test.result),ncol=4)
  ROC.box.result[t.test.indx==1,1]<-box.plot.scores
  ROC.box.result[t.test.indx==1,2]<-ROC.results[,1]
  ROC.box.result[t.test.indx==1,3]<-ROC.results[,2]
  ROC.box.result[t.test.indx==1,4]<-ROC.results[,3]
     
  colnames(ROC.box.result)<-c(paste(Yvar.name,".box.plot.score",sep=""),paste(Yvar.name,".lower.CI",sep=""),paste(Yvar.name,".AUC",sep=""),paste(Yvar.name,".higher.CI",sep=""))
   
     ###create subdirectory variable results###
     
     Variable.results<-cbind(XCMScolumns,t.test.result,ROC.box.result,Samples.comb)
     
     Variable.results<-Variable.results[which(t.test.indx==1),]
     
     write.csv(Variable.results,paste(Yvar.name,"results.csv",sep="_"),row.names=F)
     
     t.test.results<-rbind(t.test.results,t(t.test.result))
     
     ROC.box.results<-rbind(ROC.box.results,t(ROC.box.result))     
       
   }   
    
}
  
}

significant.markers<-cbind(XCMScolumns,t(t.test.results),t(ROC.box.results),Samples)
significant.markers<-significant.markers[which(apply(significant.markers[,grep("Signif.thresh",colnames(significant.markers))],1,sum)>0),]
setwd(wd) 
write.csv(significant.markers,"Features_Above_threshold.csv")

}
