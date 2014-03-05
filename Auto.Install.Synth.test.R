Auto.Install.Synth.test<-function(){
  
  ###use current working directory###
  Parent.wd<-paste(getwd(),"\\MetMSLine.Installation\\",sep="")
  ###source data directly from Github###
  require(devtools)
  require(httr)
  
  source_GitHubData <-function(url, sep = ",", header = TRUE, row.names=NULL)
  {
    request <- GET(url)
    stop_for_status(request)
    handle <- textConnection(content(request, as = 'text'))
    on.exit(close(handle))
    read.table(handle, sep = sep, header = header, row.names=row.names)
   }
  
  ###Address of Y synthetic data on GitHub###
  Y.UrlAddress<-"http://bit.ly/17VCYS6"
  
  ###Address of Metabolite synthetic data on GitHub###
  Metab.UrlAddress<-"http://bit.ly/IiSdc4"
  
  ###Address of X synthetic MS data on GitHub###
  X.UrlAddress<-"http://bit.ly/1a3vXJk"
  
  ###Download data and write to .csv files in working directory###
  write.csv(Synthetic_Y_data<-source_GitHubData(url = Y.UrlAddress,row.names=1),"Synthetic_Y_data.csv")
  ###Download Synthetic metabolite data and write to .csv in workinng directory###
  write.csv(Synthetic_Metabolite_DB<-source_GitHubData(url = Metab.UrlAddress,row.names=1),"Synthetic_Metabolite_DB.csv")
    
  ##create XCMS directory###
  XCMS_dir<-paste(Parent.wd,"XCMS\\",sep="")
  
  dir.create(XCMS_dir)
  
  setwd(XCMS_dir)
  ###write X synthetic matrix to XCMS folder from GitHub###
  write.csv(Synthetic_X_data<-source_GitHubData(url = X.UrlAddress, sep="\t"),row.names=FALSE,"Synthetic_X_data.csv")
    
  ###R scripts directory create###
  R.scripts.dir<-paste(Parent.wd,"MetMSLine.functions\\",sep="")
  dir.create(R.scripts.dir)
  
  ###Remove source_GitHubData function###
	rm(source_GitHubData)

  ###PreProc.QC.RLSC function data and save function###
  PreProc.QC.RLSC<-"http://bit.ly/1bl0h8H"
  eval(parse(file=PreProc.QC.RLSC))
  
  
  ###Perform PreProc.QC.LSC function on synthetic data###
  PreProc.QC.RLSC(X="Synthetic_X_data.csv",wd=Parent.wd,
                  XCMS_dir=XCMS_dir)
  
  setwd(R.scripts.dir)
  dump(c(lsf.str()), file="PreProc.QC.RLSC.R")
  rm(PreProc.QC.RLSC)
  
  ###Auto.PCA function data and save function###
  Auto.PCA<-"http://bit.ly/1bnm5QX"
  eval(parse(file=Auto.PCA))
  
  ###Perform Auto.PCA function on synthetic data###
  Auto.PCA(Yvar="Synthetic_Y_data.csv",wd=Parent.wd)
  
  setwd(R.scripts.dir)
  dump(c(lsf.str()), file="Auto.PCA.R")
  rm(Auto.PCA)
  
  ###Auto.MV.Regress function data and save function###
  Auto.MV.Regress<-"http://bit.ly/IgITWd"
  eval(parse(file=Auto.MV.Regress))
  
  ###Perform Auto.MV.Regress function on synthetic data###
  Auto.MV.Regress(wd=Parent.wd)
  
  setwd(R.scripts.dir)
  dump(c(lsf.str()), file="Auto.MV.Regress.R")
  rm(Auto.MV.Regress)
 
 
  ###Auto.MS.MS.match function data and save function###
  Auto.MS.MS.match<-"http://bit.ly/1fIvqES"
  eval(parse(file=Auto.MS.MS.match))
   
  setwd(R.scripts.dir)
  dump(c(lsf.str()), file="Auto.MS.MS.match.R")
  rm(Auto.MS.MS.match)
  
  ###DBAnnotate function data and save function###
  DBAnnotate<-"http://bit.ly/18nxCPU"
  eval(parse(file=DBAnnotate))
  
  
  ###Perform DBAnnotate function on synthetic data###
  DBAnnotate(database="Synthetic_Metabolite_DB.csv",wd=Parent.wd,unknowns.dir=paste(Parent.wd,"Auto.MV.Regress.results\\",sep=""))
  
  setwd(R.scripts.dir)
  dump(c(lsf.str()), file="DBAnnotate.R")
  rm(DBAnnotate)
    
 }



           
