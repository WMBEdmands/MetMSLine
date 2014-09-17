Auto.MS.MS.match<-function(MSfeatures="Features_Above_threshold.csv",  mode="Negative", TICfilter=5000,Precursor.ppm=10,Frag.ppm=20,ret=5, Parent.tol=0.1, Fragment.tol=0.5,wd="~/STUDY NAME/",mzXML.dir="~/STUDY NAME/MS_MS_mzXML/"){                                                  ###minIntensity = minimum MS MS fragment intensities
    ###minIntensity = minimum MS MS fragment intensities
    ##HMDB MS MS parent search tolerance. Set a broad tolerance here.
    ##mode for HMDB MS MS search                           
    ###HMDB MS MS fragment search tolerance.   
    require(readMzXmlData)
    
    ###empty data frame for collating MS MS matching results
    Precursor.Inter.frag.names<-data.frame()
    MS.MS.Matched.MS1.features<-data.frame() 
    MS1.features.MS.MS.matched<-data.frame()
    conjugates.mass.losses<-data.frame()
    
    
    setwd(wd)
    ####significant features######
    
    message("Reading MS1 feature table...PLEASE WAIT")#,quote=F)
    flush.console()
    
    MS1.features<-read.csv(MSfeatures,header=T)
    Original.MS1.features<-MS1.features
    
    
    message("...Done")#,quote=F)
    flush.console()
    
    setwd(mzXML.dir)


    ###save all parameters used in a dated .csv file for future reference###
	Parameters<-data.frame(MSfeatures,mode, working.dir=wd,MS.MS.file.dir=mzXML.dir,TICfilter,
	Precursor.ppm,Frag.ppm,RT.tol=ret,HMDB.Parent.delta=Parent.tol,HMDB.Parent.delta=Fragment.tol)

	date<-Sys.time()
	date<-gsub("-",".",date)
	write.csv(Parameters,paste("Parameters",substr(date,1,10),".csv",sep=" "),row.names=FALSE)

    ###identify all mzXML files in raw-data directory###
    files = list.files(pattern = "*.mzXML",full.names=F)
	
	message(paste(length(files)," MS2 .mzXML files were detected",sep=""))#,quote=F)
	flush.console()

  EIC.column.name<-colnames(MS1.features)[1]
	mzmed.column.name<-colnames(MS1.features)[2]
	RTmed.column.name<-colnames(MS1.features)[3]  
    
    for (i in 1:length(files)) {
      
      Auto.MS.MS.file.name<-as.character(files[i])
      
      message(paste("Loading ",Auto.MS.MS.file.name,"...PLEASE WAIT",sep=""))#,quote=F)
      flush.console()
      
      Auto.MS.MS.file<-readMzXmlFile(files[i])
    
      message("..Done")#,quote=F)
      flush.console()
      
      EIC_nums<-MS1.features[,EIC.column.name]
      Med.mz<-MS1.features[,mzmed.column.name]
      Med.RT<-MS1.features[,RTmed.column.name]
      
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
      MS_MS_Scan.num.df<-data.frame()
      
      for (i in 1:length(Auto_MS_MS.Intensity.filtered)) {
        PrecursorMz<-Auto_MS_MS.Intensity.filtered[[i]][[2]]$precursorMz
        collisionEnergy<-Auto_MS_MS.Intensity.filtered[[i]][[2]]$collisionEnergy
        retentionTime<-Auto_MS_MS.Intensity.filtered[[i]][[2]]$retentionTime
        MS_MS_Scan.num<-Auto_MS_MS.Intensity.filtered[[i]][[2]]$num
        
        PrecursorMz.df<-rbind(PrecursorMz.df,PrecursorMz)
        collisionEnergy.df<-rbind(collisionEnergy.df,collisionEnergy)
        retentionTime.df<-rbind(retentionTime.df,retentionTime)
        MS_MS_Scan.num.df<-rbind(MS_MS_Scan.num.df,MS_MS_Scan.num)
      }
      
      Unique.ID<-seq(1,length(Auto_MS_MS.Intensity.filtered),1)
      
      MS_MS_spectra.index<-cbind(Unique.ID,MS_MS_Scan.num.df,PrecursorMz.df,collisionEnergy.df,retentionTime.df)
      colnames(MS_MS_spectra.index)<-c("Unique.ID","MS_MS_Scan.num","PrecursorMz","collisionEnergy","retentionTime")
      
      ######search match function#####
      PrecursorMz<-as.numeric(MS_MS_spectra.index[,"PrecursorMz"])
      Unique.ID<-as.numeric(MS_MS_spectra.index[,"Unique.ID"])
      MS_MS_Scan.num<-as.numeric(MS_MS_spectra.index[,"MS_MS_Scan.num"])
      collisionEnergy<-as.numeric(MS_MS_spectra.index[,"collisionEnergy"])
      retentionTime<-as.numeric(MS_MS_spectra.index[,"retentionTime"])
      
      message("Matching MS1 features to MS2 spectra precursors...PLEASE WAIT")#,quote=F)
      flush.console()
      
      pb<-txtProgressBar(min=0,max=length(Med.mz),style=3)#,width=300)#title="Scatter plot progress bar"
      
      Match.results <- data.frame() #empty dataframe for DB search results
      
      for(i in 1:length(Med.mz)){
   
        Sys.sleep(0.1)
        setTxtProgressBar(pb,i)
        # setTkProgressBar(pb,i,label=paste( round(i/nrow(RawfoldCV)*100, 0),"% done"))
        flush.console()
        
        max = Med.mz[i]+((Precursor.ppm/1000000)*Med.mz[i])
        min = Med.mz[i]-((Precursor.ppm/1000000)*Med.mz[i])
        for (k in 1:length(PrecursorMz)) {
          index <- which (PrecursorMz[k]<max & PrecursorMz[k]>min)
          
          if( length(index) > 0) {
            Matchmass<-PrecursorMz[k][index]
            Scan.no <- MS_MS_Scan.num[k][index]
            Collision.energy<-collisionEnergy[k][index]
            MS_MS.RT<-retentionTime[k][index]
            IDs<-Unique.ID[k][index]
            
            
            identity <- cbind(EIC_nums[i],Med.mz[i],Med.RT[i],Scan.no,IDs,Matchmass,Collision.energy,MS_MS.RT)
            Match.results <- rbind(Match.results,identity)
          }
        }
      }
      
      if (nrow(Match.results)!=0 & nrow(Match.results)>1){
        r.df <- as.matrix(Match.results)
        r.dfcolnames<-c("MS1.feat.ID","MS1.feat.mz","MS1.feat.RT","MS_MS_Scan.num","Precursor.no","PrecursorMz","collision.energy","MS_MS.RT")
        colnames(r.df)<-r.dfcolnames
        
        
        RT.MS_MSmin<-r.df[,"MS_MS.RT"]-ret
        RT.MS_MSmax<-r.df[,"MS_MS.RT"]+ret
        
        RT.match.index<-(r.df[,"MS1.feat.RT"]<RT.MS_MSmax & r.df[,"MS1.feat.RT"]>RT.MS_MSmin)
        
        
        sum.matches<-sum(RT.match.index)
        
        
        if (sum.matches!=0 & sum.matches>1){
          r.df<-r.df[RT.match.index,]
          MS_MS_matched.index<-r.df[,"Precursor.no"]
          Matched.Auto.MS.MS<-Auto_MS_MS.Intensity.filtered[MS_MS_matched.index] ##Mz and RT matched list
          
          Intensity.results <- data.frame() #empty dataframe for intensity results
          
          ####max intensity fragment filtration####
          
          for (i in 1:length(Matched.Auto.MS.MS)){
            maxintensity <- max(Matched.Auto.MS.MS[[i]][[1]]$intensity)
            Intensity.results <-rbind(Intensity.results,maxintensity)#empty dataframe for intensity results
          }
          Intensity.fragments.index<-which(Intensity.results>=50)
          
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
              
              MS1.features.remIndex<-EIC_nums %in% EIC.matched.sign.feat
              MS1.features.matched<-MS1.features[MS1.features.remIndex==TRUE,]
              MS1.features<-MS1.features[MS1.features.remIndex==FALSE,] ###remove matched features from starting significant features
              
              #######Collating MS MS matched features######
              MS1.features.MS.MS.matched<-rbind(MS1.features.MS.MS.matched,MS1.features.matched)
              
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
              
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.diff,2,function(x){ifelse(x<(a<-0.984015583)+((Frag.ppm/1000000)*a) & x>(b<-0.984015583)-((Frag.ppm/1000000)*b),1,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-2.015650064)+((Frag.ppm/1000000)*a) & x>(b<-2.015650064)-((Frag.ppm/1000000)*b),2,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-4.031300128)+((Frag.ppm/1000000)*a) & x>(b<-4.031300128)-((Frag.ppm/1000000)*b),3,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-12.036385508)+((Frag.ppm/1000000)*a) & x>(b<-12.036385508)-((Frag.ppm/1000000)*b),4,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-13.979264556)+((Frag.ppm/1000000)*a) & x>(b<-13.979264556)-((Frag.ppm/1000000)*b),5,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-14.015650064)+((Frag.ppm/1000000)*a) & x>(b<-14.015650064)-((Frag.ppm/1000000)*b),6,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-15.023475096)+((Frag.ppm/1000000)*a) & x>(b<-15.023475096)-((Frag.ppm/1000000)*b),7,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-15.99491462)+((Frag.ppm/1000000)*a) & x>(b<-15.99491462)-((Frag.ppm/1000000)*b),8,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-17.026549)+((Frag.ppm/1000000)*a) & x>(b<-17.026549)-((Frag.ppm/1000000)*b),9,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-18.010565)+((Frag.ppm/1000000)*a) & x>(b<-18.010565)-((Frag.ppm/1000000)*b),10,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-21.981945)+((Frag.ppm/1000000)*a) & x>(b<-21.981945)-((Frag.ppm/1000000)*b),11,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-27.994915)+((Frag.ppm/1000000)*a) & x>(b<-27.994915)-((Frag.ppm/1000000)*b),12,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-28.00614801)+((Frag.ppm/1000000)*a) & x>(b<-28.00614801)-((Frag.ppm/1000000)*b),13,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-28.031300128)+((Frag.ppm/1000000)*a) & x>(b<-28.031300128)-((Frag.ppm/1000000)*b),14,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-28.990163592)+((Frag.ppm/1000000)*a) & x>(b<-28.990163592)-((Frag.ppm/1000000)*b),15,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-29.002739652)+((Frag.ppm/1000000)*a) & x>(b<-29.002739652)-((Frag.ppm/1000000)*b),16,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-29.03912516)+((Frag.ppm/1000000)*a) & x>(b<-29.03912516)-((Frag.ppm/1000000)*b),17,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-29.974179175)+((Frag.ppm/1000000)*a) & x>(b<-29.974179175)-((Frag.ppm/1000000)*b),18,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-30.010564684)+((Frag.ppm/1000000)*a) & x>(b<-30.010564684)-((Frag.ppm/1000000)*b),19,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-31.018498)+((Frag.ppm/1000000)*a) & x>(b<-31.018498)-((Frag.ppm/1000000)*b),20,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-32.026214748)+((Frag.ppm/1000000)*a) & x>(b<-32.026214748)-((Frag.ppm/1000000)*b),21,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-37.955882)+((Frag.ppm/1000000)*a) & x>(b<-37.955882)-((Frag.ppm/1000000)*b),22,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-42.010565)+((Frag.ppm/1000000)*a) & x>(b<-42.010565)-((Frag.ppm/1000000)*b),23,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-43.005813656)+((Frag.ppm/1000000)*a) & x>(b<-43.005813656)-((Frag.ppm/1000000)*b),24,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-43.018389716)+((Frag.ppm/1000000)*a) & x>(b<-43.018389716)-((Frag.ppm/1000000)*b),25,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-43.0547752)+((Frag.ppm/1000000)*a) & x>(b<-43.0547752)-((Frag.ppm/1000000)*b),26,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-43.989829239)+((Frag.ppm/1000000)*a) & x>(b<-43.989829239)-((Frag.ppm/1000000)*b),27,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-44.997654271)+((Frag.ppm/1000000)*a) & x>(b<-44.997654271)-((Frag.ppm/1000000)*b),28,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-46.00548)+((Frag.ppm/1000000)*a) & x>(b<-46.00548)-((Frag.ppm/1000000)*b),29,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-46.0419)+((Frag.ppm/1000000)*a) & x>(b<-46.0419)-((Frag.ppm/1000000)*b),30,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-56.062600257)+((Frag.ppm/1000000)*a) & x>(b<-56.062600257)-((Frag.ppm/1000000)*b),31,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-57.070425289)+((Frag.ppm/1000000)*a) & x>(b<-57.070425289)-((Frag.ppm/1000000)*b),32,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-59.013304335)+((Frag.ppm/1000000)*a) & x>(b<-59.013304335)-((Frag.ppm/1000000)*b),33,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-59.037113785)+((Frag.ppm/1000000)*a) & x>(b<-59.037113785)-((Frag.ppm/1000000)*b),34,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-60.021129367)+((Frag.ppm/1000000)*a) & x>(b<-60.021129367)-((Frag.ppm/1000000)*b),35,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-63.961904)+((Frag.ppm/1000000)*a) & x>(b<-63.961904)-((Frag.ppm/1000000)*b),36,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-75.032029)+((Frag.ppm/1000000)*a) & x>(b<-75.032029)-((Frag.ppm/1000000)*b),37,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-79.956814859)+((Frag.ppm/1000000)*a) & x>(b<-79.956814859)-((Frag.ppm/1000000)*b),38,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-86.000395)+((Frag.ppm/1000000)*a) & x>(b<-86.000395)-((Frag.ppm/1000000)*b),39,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-89.047679)+((Frag.ppm/1000000)*a) & x>(b<-89.047679)-((Frag.ppm/1000000)*b),40,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-121.019753)+((Frag.ppm/1000000)*a) & x>(b<-121.019753)-((Frag.ppm/1000000)*b),41,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-129.042594)+((Frag.ppm/1000000)*a) & x>(b<-129.042594)-((Frag.ppm/1000000)*b),42,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-146.036779432)+((Frag.ppm/1000000)*a) & x>(b<-146.036779432)-((Frag.ppm/1000000)*b),43,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-146.069142189)+((Frag.ppm/1000000)*a) & x>(b<-146.069142189)-((Frag.ppm/1000000)*b),44,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-152.010958607)+((Frag.ppm/1000000)*a) & x>(b<-152.010958607)-((Frag.ppm/1000000)*b),45,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-162.052825)+((Frag.ppm/1000000)*a) & x>(b<-162.052825)-((Frag.ppm/1000000)*b),46,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-163.030318)+((Frag.ppm/1000000)*a) & x>(b<-163.030318)-((Frag.ppm/1000000)*b),47,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-176.032087974)+((Frag.ppm/1000000)*a) & x>(b<-176.032087974)-((Frag.ppm/1000000)*b),48,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-176.047344115)+((Frag.ppm/1000000)*a) & x>(b<-176.047344115)-((Frag.ppm/1000000)*b),49,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-178.041213189)+((Frag.ppm/1000000)*a) & x>(b<-178.041213189)-((Frag.ppm/1000000)*b),50,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-180.06339)+((Frag.ppm/1000000)*a) & x>(b<-180.06339)-((Frag.ppm/1000000)*b),51,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-192.027)+((Frag.ppm/1000000)*a) & x>(b<-192.027)-((Frag.ppm/1000000)*b),52,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-194.042655)+((Frag.ppm/1000000)*a) & x>(b<-194.042655)-((Frag.ppm/1000000)*b),53,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-198.014035)+((Frag.ppm/1000000)*a) & x>(b<-198.014035)-((Frag.ppm/1000000)*b),54,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-203.07937252)+((Frag.ppm/1000000)*a) & x>(b<-203.07937252)-((Frag.ppm/1000000)*b),55,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-215.994605)+((Frag.ppm/1000000)*a) & x>(b<-215.994605)-((Frag.ppm/1000000)*b),56,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-221.089937203)+((Frag.ppm/1000000)*a) & x>(b<-221.089937203)-((Frag.ppm/1000000)*b),57,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-239.993994)+((Frag.ppm/1000000)*a) & x>(b<-239.993994)-((Frag.ppm/1000000)*b),58,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-255.988909)+((Frag.ppm/1000000)*a) & x>(b<-255.988909)-((Frag.ppm/1000000)*b),59,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-258.004564)+((Frag.ppm/1000000)*a) & x>(b<-258.004564)-((Frag.ppm/1000000)*b),60,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-273.096087)+((Frag.ppm/1000000)*a) & x>(b<-273.096087)-((Frag.ppm/1000000)*b),61,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-273.999479)+((Frag.ppm/1000000)*a) & x>(b<-273.999479)-((Frag.ppm/1000000)*b),62,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-275.111737)+((Frag.ppm/1000000)*a) & x>(b<-275.111737)-((Frag.ppm/1000000)*b),63,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-307.083811)+((Frag.ppm/1000000)*a) & x>(b<-307.083811)-((Frag.ppm/1000000)*b),64,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-352.06418)+((Frag.ppm/1000000)*a) & x>(b<-352.06418)-((Frag.ppm/1000000)*b),65,x)}))
              Precursor.Inter.frag.name<-(apply(Precursor.Inter.frag.name,2,function(x){ifelse(x<(a<-388.08531)+((Frag.ppm/1000000)*a) & x>(b<-388.08531)-((Frag.ppm/1000000)*b),66,x)}))
              
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
              MS.MS.Scan.num.plot<-r.df[,4]
              
              if(length(Matched.Auto.MS.MS)>1){
                
                Plot.file.url<-data.frame()
                
                message(paste("Generating plots ",Auto.MS.MS.file.name,"...PLEASE WAIT",sep=""))#,quote=F)
                flush.console()
                
                pb<-txtProgressBar(min=0,max=length(Matched.Auto.MS.MS),style=3)#,width=300)#title="Scatter plot progress bar"
          
                for (i in 1:length(Matched.Auto.MS.MS)){
                  
                  Sys.sleep(0.1)
                  setTxtProgressBar(pb,i)
                  # setTkProgressBar(pb,i,label=paste( round(i/nrow(RawfoldCV)*100, 0),"% done"))
                  flush.console()
                  
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
                  
                  MS.MS.Scan.num<-MS.MS.Scan.num.plot[i]
                  
                  Intensity.plot<-round(Frag.Intense,digits=4)
                  Mass.plot<-round(Frag.mass,digits=4)
                  MS.MS.plot.lim.index<-which(Mass.plot<maxXlimit)
                  Mass.plot<-Mass.plot[MS.MS.plot.lim.index]
                  Intensity.plot<-Intensity.plot[MS.MS.plot.lim.index]
                  XCMS.feature<-paste("M",XCMS.mz[i],"T",XCMS.rt[i],sep="")
                  
                  plot.file.name<-paste(Auto.MS.MS.file.name,MS.MS.Scan.num,"ID",EIC.name,XCMS.feature,"eV",Collision.energy.plot,sep="_")
                  plot.file.name<-paste(plot.file.name,".png",sep="")
                  
                  Plot.url<-as.data.frame(paste(mzXML.dir,plot.file.name,sep="")) ##url for corresponding MS MS file 
                  Plot.file.url<-rbind(Plot.file.url,Plot.url)
                  
                  
                  
                  png(plot.file.name,width=1200,height=1200,res=275)
                  plot(Mass.plot,Intensity.plot, main=paste("ID",EIC.name,XCMS.feature,"eV", Collision.energy.plot),sub=paste("Precursor",Precursor.name.plot,Auto.MS.MS.file.name),xlab="m/z", ylab="MS_intensity",type="h",xlim=c(0,maxXlimit),ylim=c(0,(max(Intensity.plot)*1.25)))
                  points(Top.mass.points,Top.Intensity.points,type="p",pch=23,bg="red")
                  points(Precursor.name.plot,40,type="p",pch=23,bg="blue")
                  text(Top.mass.points,(Top.Intensity.points+(max(Top.Intensity.points)*0.064)),format(Top.mass.points,digits=6),pos=3,cex=0.5)
                  text(Inter.frag.position,Inter.frag.height,Inter.frag.labels,pos=3,col="red",cex=0.5)
                  text(Top.mass.points,(Top.Intensity.points+(max(Top.Intensity.points)*0.032)),Parent.neutral.losses,pos=3,col="blue",cex=0.5)
                  dev.off()
                  
                }
                
                message("...Done")#,quote=F)
                flush.console()
                
              }
              Matched.MS.MS.fragment.df<-cbind(Plot.file.url,Matched.MS.MS.fragment.df)
              colnames(Matched.MS.MS.fragment.df)[1]<-"MS.MS.plot.file.url"
              Precursor.Inter.frag.names<-rbind(Precursor.Inter.frag.names,Precursor.Inter.frag.name) 
              MS.MS.Matched.MS1.features<-rbind(MS.MS.Matched.MS1.features,Matched.MS.MS.fragment.df)
              conjugates.mass.losses<-rbind(conjugates.mass.losses,conjugates.mass.loss)
              
            }
          }
        }
      }
    }
    
    MS.MS.Matched.MS1.features<-as.data.frame(MS.MS.Matched.MS1.features)
    unique.MS1.features<-unique(MS.MS.Matched.MS1.features[,"MS1.feat.ID"]) ###unique feature matches
    
    
    #Potential.gluc<-(apply(conjugates.mass.losses,2,function(x) {x<(a<-176.032088)+((Frag.ppm/1000000)*a) & x>(b<-176.032088)-((Frag.ppm/1000000)*b)}))+0 
    #Potential.sulf<-(apply(conjugates.mass.losses,2,function(x) {x<(a<-79.95681456)+((Frag.ppm/1000000)*a) & x>(b<-79.95681456)-((Frag.ppm/1000000)*b)}))+0
    #Potential.digluc<-(apply(conjugates.mass.losses,2,function(x) {x<((a<-176.032088*2))+((Frag.ppm/1000000)*a) & x>(b<-(176.032088*2))-((Frag.ppm/1000000)*b)}))+0 
    #Potential.disulf<-(apply(conjugates.mass.losses,2,function(x) {x<((a<-79.95681456*2))+((Frag.ppm/1000000)*a) & x>(b<-(79.95681456*2))-((Frag.ppm/1000000)*b)}))+0 
    #Potential.gluc.sulf<-(apply(conjugates.mass.losses,2,function(x) {x<(a<-(176.032088+79.95681456))+((Frag.ppm/1000000)*a) & x>(b<-(176.032088+79.95681456))-((Frag.ppm/1000000)*b)}))+0 
    
    #Glucuronide<-ifelse(rowSums(Potential.gluc)>=1,"Glucuronide",0)
    #Sulfate<-ifelse(rowSums(Potential.sulf)>=1,"Sulfate",0)
    #DiGlucuronide<-ifelse(rowSums(Potential.digluc)>=1,"DiGlucuronide",0)
    #DiSulfate<-ifelse(rowSums(Potential.disulf)>=1,"DiSulfate",0)
    #GlucuronideSulfate<-ifelse(rowSums(Potential.gluc.sulf)>=1,"GlucuronideSulfate",0)
    
    #Conjugate.columns<-data.frame()
    
    #if (is.character(Glucuronide)==TRUE) {
     # Glucuronide<-as.data.frame(Glucuronide)
      #Conjugate.columns<-rbind(Conjugate.columns,Glucuronide)
    #}
    
    #if (is.character(Sulfate)==TRUE) {
     # Sulfate<-as.data.frame(Sulfate)
      #Conjugate.columns<-cbind(Conjugate.columns,Sulfate)
    #}
    
   # if (is.character(DiGlucuronide)==TRUE) {
    #  DiGlucuronide<-as.data.frame(DiGlucuronide)
     # Conjugate.columns<-cbind(Conjugate.columns,DiGlucuronide)
    #}
    
    #if (is.character(DiSulfate)==TRUE) {
     # DiSulfate<-as.data.frame(DiSulfate)
      #Conjugate.columns<-cbind(Conjugate.columns,DiSulfate)
    #}
    
    #if (is.character(GlucuronideSulfate)==TRUE) {
     # GlucuronideSulfate<-as.data.frame(GlucuronideSulfate)
      #Conjugate.columns<-cbind(Conjugate.columns,GlucuronideSulfate)
    #}
    
    #if (ncol(Conjugate.columns)!=0){
     # MS.MS.Matched.MS1.features<-cbind(MS.MS.Matched.MS1.features,Conjugate.columns)
    #}
    
    Original.EICs<-Original.MS1.features[,EIC.column.name] ###column name subject to change
    
    MS_MS.matched.dummy<-as.matrix(Original.EICs %in% unique.MS1.features)+0 ###dummy matrix matched sign features
    colnames(MS_MS.matched.dummy)<-"Auto.MS.MS.matched"
    Original.MS1.features<-cbind(Original.MS1.features,MS_MS.matched.dummy)
    
    ####conjugate and mass loss matching####
    
    ###Index matching values###
    Matched.Original.features.index<-match(MS.MS.Matched.MS1.features[,"MS1.feat.ID"],Original.MS1.features[,EIC.column.name])
    MS.MS.Matched.MS1.features<-cbind(Original.MS1.features[Matched.Original.features.index,],MS.MS.Matched.MS1.features,Precursor.Inter.frag.names)
    
    
    
    
    
    setwd(wd)
    write.csv(Original.MS1.features,MSfeatures,row.names=FALSE)
    
    setwd(mzXML.dir)
    write.csv(MS1.features,"unmatched_sign_features.csv",row.names=FALSE)
    write.csv(MS.MS.Matched.MS1.features,"MS_MS_Matched_sign_features.csv",row.names=FALSE)
    
  }

  ###END###
