Auto.MS.MS.match<-function(MSfeatures="Features_Above_threshold.csv",mode="negative",wd="D:\\R_data_processing\\STUDY NAME\\",  mzXML.dir="D:\\R_data_processing\\STUDY NAME\\MS_MS_mzXML\\",TICfilter=5000,minIntensity=500,delta=0.005,ret=5,Parent.tol=0.1,Fragment.tol=0.5){                                                
  ###minIntensity = minimum MS MS fragment intensities
  ##HMDB MS MS parent search tolerance. Set a broad tolerance here.
  ##mode for HMDB MS MS search                           
  ###HMDB MS MS fragment search tolerance.   
  require(readMzXmlData)
  
  wd<-paste(wd,"Auto.MV.Regress.results\\",sep="")
  
  ###empty data frame for collating MS MS matching results
  Precursor.Inter.frag.names<-data.frame()
  MS.MS.Matched.sign.features<-data.frame() 
  sign.features.MS.MS.matched<-data.frame()
  conjugates.mass.losses<-data.frame()
  
  
  setwd(wd)
  ####significant features######
  Sign.features<-read.csv(MSfeatures,header=T)
  Original.Sign.Features<-Sign.features
  
  
  setwd(mzXML.dir)
  files<-dir()
  Length.files<-length(files)
  files.seq<-seq(1,Length.files,1)
  file.name<-files[1]
  file.name<-substring(file.name,2)
  files<-paste(files.seq,file.name,sep="") ###OPTION OF CHANGING STARTING FILE NAME#####
  
  
  
  for (i in 1:length(files)) {
    Auto.MS.MS.file<-readMzXmlFile(files[i])
    
    Auto.MS.MS.file.name<-as.character(files[i])
    
    EIC_nums<-Sign.features[,"XCMS_EIC"]
    Med.RT<-Sign.features[,"rtmed"]
    Med.mz<-Sign.features[,"mzmed"]
    
    ###subsetting all lists containing a precursorMz###  
    MSlevel2.df<-data.frame()
    for (i in 1:length(Auto.MS.MS.file)) {
      MSlevel<-Auto.MS.MS.file[[i]][[2]]$msLevel
      MSlevel2.df<-rbind(MSlevel2.df,MSlevel)}
    
    Mslevel2.index<-which(MSlevel2.df==2)
    
    Auto_MS_MS.Mslevel2<-Auto.MS.MS.file[Mslevel2.index]#TIC intensity filtered
    
    
    
    Intensity.index<-data.frame()
    for (i in 1:length(Auto_MS_MS.Mslevel2)) {TIC<-Auto_MS_MS.Mslevel2[[i]][[2]]$totIonCurrent
                                              Intensity.index<-rbind(Intensity.index,TIC)}
    
    Intensity.filter.index<-which(Intensity.index>=TICfilter)
    Auto_MS_MS.Intensity.filtered<-Auto_MS_MS.Mslevel2[Intensity.filter.index]#TIC intensity filtered
    
    ###Precursor Mz list
     
    PrecursorMz.df<-data.frame()
    collisionEnergy.df<-data.frame()
    retentionTime.df<-data.frame()
    
    for (i in 1:length(Auto_MS_MS.Intensity.filtered)) {
      PrecursorMz<-Auto_MS_MS.Intensity.filtered[[i]][[2]]$precursorMz
      collisionEnergy<-Auto_MS_MS.Intensity.filtered[[i]][[2]]$collisionEnergy
      retentionTime<-Auto_MS_MS.Intensity.filtered[[i]][[2]]$retentionTime
      
      PrecursorMz.df<-rbind(PrecursorMz.df,PrecursorMz)
      collisionEnergy.df<-rbind(collisionEnergy.df,collisionEnergy)
      retentionTime.df<-rbind(retentionTime.df,retentionTime)}
    
    MS_MS_unique_ID<-seq(1,nrow(PrecursorMz.df),1)
    
    MS_MS_spectra.index<-cbind(MS_MS_unique_ID,PrecursorMz.df,collisionEnergy.df,retentionTime.df)
    colnames(MS_MS_spectra.index)<-c("MS_MS_unique_ID","PrecursorMz","collisionEnergy","retentionTime")
    
    ######search match function#####
    PrecursorMz<-as.numeric(MS_MS_spectra.index[,"PrecursorMz"])
    MS_MS_unique_ID<-as.numeric(MS_MS_spectra.index[,"MS_MS_unique_ID"])
    collisionEnergy<-as.numeric(MS_MS_spectra.index[,"collisionEnergy"])
    retentionTime<-as.numeric(MS_MS_spectra.index[,"retentionTime"])
    
   
    Match.results <- data.frame() #empty dataframe for DB search results
    
    for(i in 1:length(Med.mz)){
      max = Med.mz[i]+delta
      min = Med.mz[i]-delta
      for (k in 1:length(PrecursorMz)) {
        index <- which (PrecursorMz[k]<max & PrecursorMz[k]>min)
        
        if( length(index) > 0) {
          Matchmass<-PrecursorMz[k][index]
          ids <- MS_MS_unique_ID[k][index]
          Collision.energy<-collisionEnergy[k][index]
          MS_MS.RT<-retentionTime[k][index]
          
          identity <- cbind(EIC_nums[i],Med.mz[i],Med.RT[i],ids,Matchmass,Collision.energy,MS_MS.RT)
          Match.results <- rbind(Match.results,identity)
        }
      }
    }
    
    if (nrow(Match.results)!=0 & nrow(Match.results)>1){
      r.df <- as.matrix(Match.results)
      r.dfcolnames<-c("Sign.feat.EIC","Sign.feat.mz","Sign.feat.RT","MS_MS_unique_ID","PrecursorMz","collision.energy","MS_MS.RT")
      colnames(r.df)<-r.dfcolnames
      
      
      RT.MS_MSmin<-r.df[,"MS_MS.RT"]-ret
      RT.MS_MSmax<-r.df[,"MS_MS.RT"]+ret
      
      RT.match.index<-(r.df[,"Sign.feat.RT"]<RT.MS_MSmax & r.df[,"Sign.feat.RT"]>RT.MS_MSmin)
      
      
      sum.matches<-sum(RT.match.index)
      
      
      if (sum.matches!=0 & sum.matches>1){
        r.df<-r.df[RT.match.index,]
        MS_MS_matched.index<-r.df[,"MS_MS_unique_ID"]
        Matched.Auto.MS.MS<-Auto_MS_MS.Intensity.filtered[MS_MS_matched.index] ##Mz and RT matched list
        
        Intensity.results <- data.frame() #empty dataframe for intensity results
        
        ####max intensity fragment filtration####
        
        for (i in 1:length(Matched.Auto.MS.MS)){
          maxintensity <- max(Matched.Auto.MS.MS[[i]][[1]]$intensity)
          Intensity.results <-rbind(Intensity.results,maxintensity)#empty dataframe for intensity results
        }
        Intensity.fragments.index<-which(Intensity.results>=minIntensity)
        
        if (length(Intensity.fragments.index)!=0 & length(Intensity.fragments.index)>1){
          Intensity.results<-Intensity.results[Intensity.fragments.index,]
          Matched.Auto.MS.MS<-Matched.Auto.MS.MS[Intensity.fragments.index]#TIC intensity filtered
          r.df<-r.df[Intensity.fragments.index,]
          
          ############################################################################################################################################################################
          ############################################################################################################################################################################
          
          
          if(length(Matched.Auto.MS.MS)>1){
            Intensity.fragments <- data.frame()
            Mass.fragments<-data.frame()#empty dataframe for intensity results
            RI.fragments<-data.frame()##relative intensity of fragments compared to most intense
            HMDB.MS.MS.url<-data.frame()
            
            for (i in 1:length(Matched.Auto.MS.MS)){
              
              Frag.Intense<-Matched.Auto.MS.MS[[i]][[1]]$intensity
            
              Frag.mass <-Matched.Auto.MS.MS[[i]][[1]]$mass
             
              Collision.energy.plot<-Matched.Auto.MS.MS[[i]][[2]]$collisionEnergy
              
              File.origin<-rep(Auto.MS.MS.file.name,length(Matched.Auto.MS.MS))
              
              Frag.df<-as.data.frame(cbind(Frag.Intense,Frag.mass))
              Frag.df.order<-Frag.df[order(Frag.Intense,decreasing=TRUE),] #order by most intense
              Frag.df.top<-Frag.df.order[c(1,2,3,4,5,6),]
              RI.fragment<-(Frag.df.top[,1]/Frag.df.top[1,1])*100
              RI.fragment<-round(RI.fragment,digits=3)
              Frag.df.top<-cbind(Frag.df.top,RI.fragment)
              
              Mz.HMDB.search<-Frag.df.top[,2]
              Frag.df.mass.order<-Frag.df.top[order(Frag.df.top[,"Frag.mass"]),]
              Frag.df.top.intensity<-Frag.df.mass.order[,1]
              Frag.df.top.mass<-Frag.df.mass.order[,2]
              Frag.df.top.RI<-Frag.df.mass.order[,3]
              
              Precursor.name.plot<-Matched.Auto.MS.MS[[i]][[2]]$precursorMz
              Precursor.name.plot<-round(Precursor.name.plot,digits=2)
              
              #####HMDB MS MS search link creation ########
              
              HMDB.MS.MS.url.link<-as.data.frame(paste("http://www.hmdb.ca/spectra/ms_ms/search?utf8=%E2%9C%93&parent_ion_mass=",Precursor.name.plot,"&parent_ion_mass_tolerance=",Parent.tol,"&ionization_mode=",mode,"&collision_energy_level=","&peaks=",Mz.HMDB.search[1],"+",RI.fragment[1],"%0D%0A",Mz.HMDB.search[2],"+",RI.fragment[2],"%0D%0A",Mz.HMDB.search[3],"+",RI.fragment[3],"&mass_charge_tolerance=",Fragment.tol,"&commit=Search",sep=""))
              
              HMDB.MS.MS.url<-rbind(HMDB.MS.MS.url,HMDB.MS.MS.url.link)
              
              
              
              Intensity.fragments<-rbind(Intensity.fragments,Frag.df.top.intensity)
              Mass.fragments<-rbind(Mass.fragments,Frag.df.top.mass)
              RI.fragments<-rbind(RI.fragments,Frag.df.top.RI)
            }   
            
            
            Matched.MS.MS.fragment.df<-cbind(File.origin,HMDB.MS.MS.url,Mass.fragments,Intensity.fragments,RI.fragments)
            colnames(Matched.MS.MS.fragment.df)<-c("Auto_MS_MS_file","HMDB.MS.MS.url","Frag.1.mz","Frag.2.mz","Frag.3.mz","Frag.4.mz","Frag.5.mz","Frag.6.mz","Frag.1.Int","Frag.2.Int","Frag.3.Int","Frag.4.Int","Frag.5.Int","Frag.6.Int","Frag.1.RI","Frag.2.RI","Frag.3.RI","Frag.4.RI","Frag.5.RI","Frag.6.RI")
            Matched.MS.MS.fragment.df<-cbind(r.df,Matched.MS.MS.fragment.df)
            
            SignF.EICremIndex<-duplicated(Matched.MS.MS.fragment.df[,1])
            EIC.matched.sign.feat<-Matched.MS.MS.fragment.df[,1]
            EIC.matched.sign.feat<-EIC.matched.sign.feat[SignF.EICremIndex==FALSE] #####EIC removal from correlated features#####
            
            Sign.features.remIndex<-EIC_nums %in% EIC.matched.sign.feat
            Sign.features.matched<-Sign.features[Sign.features.remIndex==TRUE,]
            Sign.features<-Sign.features[Sign.features.remIndex==FALSE,] ###remove matched features from starting significant features
            
            #######Collating MS MS matched features######
            sign.features.MS.MS.matched<-rbind(sign.features.MS.MS.matched,Sign.features.matched)
  
            Matched.PrecursorMz<-as.matrix(Matched.MS.MS.fragment.df[,"PrecursorMz"])
            Matched.Frag.1.mz<-as.matrix(Matched.MS.MS.fragment.df[,"Frag.1.mz"])##first fragment
            Matched.Frag.2.mz<-as.matrix(Matched.MS.MS.fragment.df[,"Frag.2.mz"])#
            Matched.Frag.3.mz<-as.matrix(Matched.MS.MS.fragment.df[,"Frag.3.mz"])#
            Matched.Frag.4.mz<-as.matrix(Matched.MS.MS.fragment.df[,"Frag.4.mz"])#
            Matched.Frag.5.mz<-as.matrix(Matched.MS.MS.fragment.df[,"Frag.5.mz"])#
            Matched.Frag.6.mz<-as.matrix(Matched.MS.MS.fragment.df[,"Frag.6.mz"])#
            
            Matched.Frag.1.Int<-as.matrix(Matched.MS.MS.fragment.df[,"Frag.1.Int"])##first fragment
            Matched.Frag.2.Int<-as.matrix(Matched.MS.MS.fragment.df[,"Frag.2.Int"])#
            Matched.Frag.3.Int<-as.matrix(Matched.MS.MS.fragment.df[,"Frag.3.Int"])#
            Matched.Frag.4.Int<-as.matrix(Matched.MS.MS.fragment.df[,"Frag.4.Int"])#
            Matched.Frag.5.Int<-as.matrix(Matched.MS.MS.fragment.df[,"Frag.5.Int"])#
            Matched.Frag.6.Int<-as.matrix(Matched.MS.MS.fragment.df[,"Frag.6.Int"])#
            
            
            Fragments.matrix<-cbind(Matched.Frag.1.mz,Matched.Frag.2.mz,Matched.Frag.3.mz,Matched.Frag.4.mz,Matched.Frag.5.mz,Matched.Frag.6.mz)
            Intensity.matrix<-cbind(Matched.Frag.1.Int,Matched.Frag.2.Int,Matched.Frag.3.Int,Matched.Frag.4.Int,Matched.Frag.5.Int,Matched.Frag.6.Int)
            
            
            Inter.frag.differences<-cbind((Matched.Frag.2.mz-Matched.Frag.1.mz),(Matched.Frag.3.mz-Matched.Frag.2.mz),(Matched.Frag.4.mz-Matched.Frag.3.mz),(Matched.Frag.5.mz-Matched.Frag.4.mz),(Matched.Frag.6.mz-Matched.Frag.5.mz))
            Inter.frag.differences<-(apply(Inter.frag.differences,2,function(x) {ifelse(x<0,0,x)})) 
            
            Precursor.differences<-(as.vector(Matched.PrecursorMz)-Fragments.matrix)
            Precursor.differences<-(apply(Precursor.differences,2,function(x) {ifelse(x<0,0,x)})) 
            
            Precursor.Inter.frag.diff<-cbind(Precursor.differences,Inter.frag.differences)
            
            #######################FRAGMENT IDENTIFICATION############################################
            
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.diff,2,function(x){ifelse(x<(0.984015583+delta) & x>(0.984015583-delta),1,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(2.015650064+0.005) & x>(2.015650064-0.005),2,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(4.031300128+0.005) & x>(4.031300128-0.005),3,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(12.036385508+0.005) & x>(12.036385508-0.005),4,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(13.979264556+0.005) & x>(13.979264556-0.005),5,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(14.015650064+0.005) & x>(14.015650064-0.005),6,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(15.023475096+0.005) & x>(15.023475096-0.005),7,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(15.99491462+0.005) & x>(15.99491462-0.005),8,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(17.026549+0.005) & x>(17.026549-0.005),9,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(18.010565+0.005) & x>(18.010565-0.005),10,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(21.981945+0.005) & x>(21.981945-0.005),11,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(27.994915+0.005) & x>(27.994915-0.005),12,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(28.00614801+0.005) & x>(28.00614801-0.005),13,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(28.031300128+0.005) & x>(28.031300128-0.005),14,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(28.990163592+0.005) & x>(28.990163592-0.005),15,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(29.002739652+0.005) & x>(29.002739652-0.005),16,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(29.03912516+0.005) & x>(29.03912516-0.005),17,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(29.974179175+0.005) & x>(29.974179175-0.005),18,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(30.010564684+0.005) & x>(30.010564684-0.005),19,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(31.018498+0.005) & x>(31.018498-0.005),20,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(32.026214748+0.005) & x>(32.026214748-0.005),21,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(37.955882+0.005) & x>(37.955882-0.005),22,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(42.010565+0.005) & x>(42.010565-0.005),23,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(43.005813656+0.005) & x>(43.005813656-0.005),24,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(43.018389716+0.005) & x>(43.018389716-0.005),25,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(43.0547752+0.005) & x>(43.0547752-0.005),26,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(43.989829239+0.005) & x>(43.989829239-0.005),27,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(44.997654271+0.005) & x>(44.997654271-0.005),28,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(46.00548+0.005) & x>(46.00548-0.005),29,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(46.0419+0.005) & x>(46.0419-0.005),30,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(56.062600257+0.005) & x>(56.062600257-0.005),31,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(57.070425289+0.005) & x>(57.070425289-0.005),32,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(59.013304335+0.005) & x>(59.013304335-0.005),33,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(59.037113785+0.005) & x>(59.037113785-0.005),34,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(60.021129367+0.005) & x>(60.021129367-0.005),35,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(63.961904+0.005) & x>(63.961904-0.005),36,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(75.032029+0.005) & x>(75.032029-0.005),37,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(79.956814859+0.005) & x>(79.956814859-0.005),38,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(86.000395+0.005) & x>(86.000395-0.005),39,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(89.047679+0.005) & x>(89.047679-0.005),40,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(121.019753+0.005) & x>(121.019753-0.005),41,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(129.042594+0.005) & x>(129.042594-0.005),42,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(146.036779432+0.005) & x>(146.036779432-0.005),43,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(146.069142189+0.005) & x>(146.069142189-0.005),44,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(152.010958607+0.005) & x>(152.010958607-0.005),45,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(162.052825+0.005) & x>(162.052825-0.005),46,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(163.030318+0.005) & x>(163.030318-0.005),47,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(176.032087974+0.005) & x>(176.032087974-0.005),48,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(176.047344115+0.005) & x>(176.047344115-0.005),49,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(178.041213189+0.005) & x>(178.041213189-0.005),50,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(180.06339+0.005) & x>(180.06339-0.005),51,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(192.027+0.005) & x>(192.027-0.005),52,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(194.042655+0.005) & x>(194.042655-0.005),53,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(198.014035+0.005) & x>(198.014035-0.005),54,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(203.07937252+0.005) & x>(203.07937252-0.005),55,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(215.994605+0.005) & x>(215.994605-0.005),56,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(221.089937203+0.005) & x>(221.089937203-0.005),57,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(239.993994+0.005) & x>(239.993994-0.005),58,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(255.988909+0.005) & x>(255.988909-0.005),59,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(258.004564+0.005) & x>(258.004564-0.005),60,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(273.096087+0.005) & x>(273.096087-0.005),61,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(273.999479+0.005) & x>(273.999479-0.005),62,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(275.111737+0.005) & x>(275.111737-0.005),63,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(307.083811+0.005) & x>(307.083811-0.005),64,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(352.06418+0.005) & x>(352.06418-0.005),65,x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(388.08531+0.005) & x>(388.08531-0.005),66,x)}))
            
            Precursor.Inter.frag.name<-round(Precursor.Inter.frag.name,digits=4)
            
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==1,"[+O-NH3]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==2,"[Reduction]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==3,"[-H4]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==4,"[+O-C2H4]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==5,"[+H2-O]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==6,"[Demethyl]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==7,"[CH3]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==8,"[O/+O-S]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==9,"[NH3]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==10,"[H20]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==11,"[Na]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==12,"[CO]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==13,"[N2]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==14,"[C2H4] ",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==15,"[+H-NO] ",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==16,"[CHO]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==17,"[C2H5] ",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==18,"[+H2-O2]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==19,"[CH2O] ",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==20,"[OCH3]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==21,"[MeOH]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==22,"[K]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==23,"[AcKetene]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==24,"[HCNO]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==25,"[CH3CO]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==26,"[C3H7]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==27,"[CO2]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==28,"[formate/+H-NO2]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==29,"[HCOOH]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==30,"[EtOH]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==31,"[C4H8]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==32,"[C4H9]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==33,"[acetate]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==34,"[Acetamide]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==35,"[CH3COOH]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==36,"[SO2]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==37,"[Gly]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==38,"[Sulf]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==39,"[Malonyl]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==40,"[Cys conj-ala]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==41,"[Cys]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==42,"[GSH AnhydroGlu]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==43,"[Coumaroyl loss]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==44,"[Ala-Gly/ GSH Glu loss]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==45,"[Galloyl loss]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==46,"[AnhydroGlucose]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==47,"[N-AcCys]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==48,"[AnhydroGluc]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==49,"[Feruloyl loss]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==50,"[Cys-Gly loss]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==51,"[Gluc loss]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==52,"[Hydroxylation+Gluc]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==53,"[Gluc]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==54,"[AnhydroGluc + Na]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==55,"[AnhydroGlucNAc]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==56,"[Gluc + Na]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==57,"[GlcNAc loss]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==58,"[AnhydroGluc + SO2]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==59,"[AnhydroGluc + SO3]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==60,"[Gluc + SO2]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==61,"[GSH GluAlaGly-2H]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==62,"[Gluc + SO3]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==63,"[GSH-GluAlaGly]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==64,"[GSH]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==65,"[DiAnhydroGlu]",x)}))
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x==66,"[DiGlu]",x)}))
            
            Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){paste("-",x,sep="")}))
            colnames(Precursor.Inter.frag.name)<-c("Precursor.Frag.diff.1","Precursor.Frag.diff.2","Precursor.Frag.diff.3","Precursor.Frag.diff.4","Precursor.Frag.diff.5","Precursor.Frag.diff.6","Frag.mz.2-Frag.mz.1","Frag.mz.3-Frag.mz.2","Frag.mz.4-Frag.mz.3","Frag.mz.5-Frag.mz.4","Frag.mz.6-Frag.mz.5")
            #################################END INTER-FRAGMENT IDENTIFICATION#####################################################
            
            conjugates.mass.loss<-cbind(Fragments.matrix,(as.vector(Matched.PrecursorMz)-Fragments.matrix))
                      
            ######################PLOT CREATION LOOP############
            EIC.name.plots<-r.df[,1]
            XCMS.mz<-round(r.df[,2],digits=0)
            XCMS.rt<-round(r.df[,3],digits=0)
            MS.MS.unique.ID.plot<-r.df[,4]
            
            if(length(Matched.Auto.MS.MS)>1){
              
              Plot.file.url<-data.frame()
              
              for (i in 1:length(Matched.Auto.MS.MS)){
                
                Frag.Intense<-Matched.Auto.MS.MS[[i]][[1]]$intensity
                # Frag.Intense<-Matched.Auto.MS.MS[[1]][[1]]$intensity
                Frag.mass <-Matched.Auto.MS.MS[[i]][[1]]$mass
                #  Frag.mass <-Matched.Auto.MS.MS[[1]][[1]]$mass
                Collision.energy.plot<-Matched.Auto.MS.MS[[i]][[2]]$collisionEnergy
                
                Top.mass.points<-Fragments.matrix[i,]
                Top.Intensity.points<- Intensity.matrix[i,]
                
                
                Fragment.diff.plot<-Precursor.Inter.frag.name[i,]
                
                Inter.frag.position<-Top.mass.points[1:5]
                Inter.frag.height<-Top.Intensity.points[1:5]
                Inter.frag.labels<-Fragment.diff.plot[7:11]
                Parent.neutral.losses<-Fragment.diff.plot[1:6]
                
                Precursor.name.plot<-Matched.Auto.MS.MS[[i]][[2]]$precursorMz
                Precursor.name.plot<-round(Precursor.name.plot,digits=2)
                
                maxXlimit<-Precursor.name.plot+250
                
                EIC.name<-EIC.name.plots[i]
                
                MS.MS.unique.ID<-MS.MS.unique.ID.plot[i]
                
                Intensity.plot<-round(Frag.Intense,digits=4)
                Mass.plot<-round(Frag.mass,digits=4)
                MS.MS.plot.lim.index<-which(Mass.plot<maxXlimit)
                Mass.plot<-Mass.plot[MS.MS.plot.lim.index]
                Intensity.plot<-Intensity.plot[MS.MS.plot.lim.index]
                XCMS.feature<-paste("M",XCMS.mz[i],"T",XCMS.rt[i],sep="")
                
                plot.file.name<-paste(Auto.MS.MS.file.name,MS.MS.unique.ID,"XCMS_EIC",EIC.name,XCMS.feature,"eV",Collision.energy.plot,sep="_")
                plot.file.name<-paste(plot.file.name,".png",sep="")
                
                Plot.url<-as.data.frame(paste(mzXML.dir,plot.file.name,sep="")) ##url for corresponding MS MS file 
                Plot.file.url<-rbind(Plot.file.url,Plot.url)
                
                
                
                png(plot.file.name,width=1200,height=1200,res=275)
                plot(Mass.plot,Intensity.plot, main=paste("XCMS_EIC",EIC.name,XCMS.feature,"eV", Collision.energy.plot),sub=paste("Precursor",Precursor.name.plot,Auto.MS.MS.file.name),xlab="m/z", ylab="MS_intensity",type="h",xlim=c(0,maxXlimit),ylim=c(0,(max(Intensity.plot)*1.25)))
                points(Top.mass.points,Top.Intensity.points,type="p",pch=23,bg="red")
                points(Precursor.name.plot,40,type="p",pch=23,bg="blue")
                text(Top.mass.points,(Top.Intensity.points+(max(Top.Intensity.points)*0.064)),format(Top.mass.points,digits=6),pos=3,cex=0.5)
                text(Inter.frag.position,Inter.frag.height,Inter.frag.labels,pos=3,col="red",cex=0.5)
                text(Top.mass.points,(Top.Intensity.points+(max(Top.Intensity.points)*0.032)),Parent.neutral.losses,pos=3,col="blue",cex=0.5)
                dev.off()
                
              }
            }
            Matched.MS.MS.fragment.df<-cbind(Plot.file.url,Matched.MS.MS.fragment.df)
            colnames(Matched.MS.MS.fragment.df)[1]<-"MS.MS.plot.file.url"
            Precursor.Inter.frag.names<-rbind(Precursor.Inter.frag.names,Precursor.Inter.frag.name) 
            MS.MS.Matched.sign.features<-rbind(MS.MS.Matched.sign.features,Matched.MS.MS.fragment.df)
            conjugates.mass.losses<-rbind(conjugates.mass.losses,conjugates.mass.loss)
            
          }
        }
      }
    }
  }
  #}  
  MS.MS.Matched.sign.features<-as.data.frame(MS.MS.Matched.sign.features)
  unique.Sign.features<-unique(MS.MS.Matched.sign.features[,"Sign.feat.EIC"]) ###unique feature matches
  
  
  Potential.gluc<-(apply(conjugates.mass.losses,2,function(x) {x<(176.032088+delta) & x>(176.032088-delta)}))+0 
  Potential.sulf<-(apply(conjugates.mass.losses,2,function(x) {x<(79.95681456+delta) & x>(79.95681456-delta)}))
  Potential.digluc<-(apply(conjugates.mass.losses,2,function(x) {x<((176.032088*2)+delta) & x>((176.032088*2)-delta)}))+0 
  Potential.disulf<-(apply(conjugates.mass.losses,2,function(x) {x<((79.95681456*2)+delta) & x>((79.95681456*2)-delta)}))+0 
  Potential.gluc.sulf<-(apply(conjugates.mass.losses,2,function(x) {x<((176.032088+79.95681456)+delta) & x>((176.032088+79.95681456)-delta)}))+0 
  
  Glucuronide<-ifelse(rowSums(Potential.gluc)>=1,"Glucuronide",0)
  Sulfate<-ifelse(rowSums(Potential.sulf)>=1,"Sulfate",0)
  DiGlucuronide<-ifelse(rowSums(Potential.digluc)>=1,"DiGlucuronide",0)
  DiSulfate<-ifelse(rowSums(Potential.disulf)>=1,"DiSulfate",0)
  GlucuronideSulfate<-ifelse(rowSums(Potential.gluc.sulf)>=1,"GlucuronideSulfate",0)
  
  Conjugate.columns<-data.frame()
  
  if (is.character(Glucuronide)==TRUE) {
    Glucuronide<-as.data.frame(Glucuronide)
    Conjugate.columns<-rbind(Conjugate.columns,Glucuronide)
  }
  
  if (is.character(Sulfate)==TRUE) {
    Sulfate<-as.data.frame(Sulfate)
    Conjugate.columns<-cbind(Conjugate.columns,Sulfate)
  }
  
  if (is.character(DiGlucuronide)==TRUE) {
    DiGlucuronide<-as.data.frame(DiGlucuronide)
    Conjugate.columns<-cbind(Conjugate.columns,DiGlucuronide)
  }
  
  if (is.character(DiSulfate)==TRUE) {
    DiSulfate<-as.data.frame(DiSulfate)
    Conjugate.columns<-cbind(Conjugate.columns,DiSulfate)
  }
  
  if (is.character(GlucuronideSulfate)==TRUE) {
    GlucuronideSulfate<-as.data.frame(GlucuronideSulfate)
    Conjugate.columns<-cbind(Conjugate.columns,GlucuronideSulfate)
  }
  
  if (ncol(Conjugate.columns)!=0){
    MS.MS.Matched.sign.features<-cbind(MS.MS.Matched.sign.features,Conjugate.columns)
  }
  
  Original.EICs<-Original.Sign.Features[,"XCMS_EIC"] ###column name subject to change
  
  MS_MS.matched.dummy<-as.matrix(Original.EICs %in% unique.Sign.features)+0 ###dummy matrix matched sign features
  colnames(MS_MS.matched.dummy)<-"Auto.MS.MS.matched"
  Original.Sign.Features<-cbind(Original.Sign.Features,MS_MS.matched.dummy)
    
  ####conjugate and mass loss matching####
    
  ###Index matching values###
  Matched.Original.features.index<-match(MS.MS.Matched.sign.features[,"Sign.feat.EIC"],Original.Sign.Features[,"XCMS_EIC"])
  MS.MS.Matched.sign.features<-cbind(Original.Sign.Features[Matched.Original.features.index,],MS.MS.Matched.sign.features,Precursor.Inter.frag.names)
  
  
 
  
  
  setwd(wd)
  write.csv(Original.Sign.Features,MSfeatures,row.names=FALSE)
  
  setwd(mzXML.dir)
  write.csv(Sign.features,"unmatched_sign_features.csv",row.names=FALSE)
  write.csv(MS.MS.Matched.sign.features,"MS_MS_Matched_sign_features.csv",row.names=FALSE)
  
}

###END###