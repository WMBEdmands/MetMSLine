DBAnnotate<-function(X="Features_above_threshold.csv",database="metabolite_DB.csv",mode="negative",conjugates="NO",MassAcc=10,unknowns.dir="~/STUDY NAME/",wd="~/STUDY NAME/"){
  # input XCMS diff report and tentative metabolite list with monoisotopic masses, also require
  # ionisation mode, anticipated acceptable mass accuracy (dependent on mass spectrometer performance) and retention time window for adduct determination returns txt file in working
  setwd(unknowns.dir)
  
  message("Reading MS1 feature table...PLEASE WAIT")#,quote=F)
  flush.console()
  
  sample<-read.csv(X,header=T)
  
  message("...Done")#,quote=F)
  flush.console()
  
  EIC.column.name<-colnames(sample)[1]
  
  #sample<-sample[,-1]
  setwd(wd)
  MetaboliteData<-data.frame(read.csv(database,header=T))
  setwd(unknowns.dir)
  
  addColumns<-length(MetaboliteData)-2
   
  ParentIDno<-seq(1,length(MetaboliteData[,1]),1) #Parent unique ID number
  MetaboliteData<-cbind(ParentIDno,MetaboliteData) #Parent ID number concatenated before conjugate and adduct calculation

  if (conjugates=="YES")
  {
  conjugates=c("Gluc","Sulf","DiGluc","DiSulf","GlucSulf","NAcCys","Glyc")
  conjugatesID<-c("Gluc","Sulf","DiGluc","DiSulf","GlucSulf","NAcCys","Glyc") #all conjugates more need to be added here when more conjugates are considered
  conjugatesmatch<-as.numeric(which(conjugatesID%in%conjugates)) # matches user input conjugates with complete ID vector
  
  #Phase 2 conjugation table creation:: according to the character vector conjugates provided in the function.
  Conjugate_name_label<-c("_glucuronide","_sulfate","_diglucuronide","_diSulfate","_glucuronidesulfate","_Nacetylcysteine","_Glycine")
  Conjugate_column<-c("Glucuronide","Sulfate","DiGlucuronide","DiSulfate","GlucuronideSulfate","Nacetylcysteine","Glycine")
  Conjugate_name_labels<-Conjugate_name_label[conjugatesmatch] # selecting conjugate text to concatenate
  Conjugate_column<-Conjugate_column[conjugatesmatch]
  
  # Monomass conjugate calculation and matching
  Gluc=176.032088 
  Sulf=79.95681456
  DiGluc=Gluc*2
  DiSulf=Sulf*2  
  GlucSulf=Gluc+Sulf
  NAcCys=163.030318
  Glyc=74.024753

  Mono_mass_conjugate<-c(Gluc,Sulf,DiGluc,DiSulf,GlucSulf,NAcCys,Glyc) #monoisotopic mass for conjugate
  Conjugate_mono_masses<-Mono_mass_conjugate[conjugatesmatch] #match mono masses to be utilised
  
  ConjugateCharVector<-rep(Conjugate_column, each=length(ParentIDno))#character vector for concatenation
  ConjugateLabelVector<-rep( Conjugate_name_labels, each=length(ParentIDno))#character vector for concatenation
  ParentCharVector<-rep("Parent",each=length(ParentIDno)) # Parent names
  Conjugate_type<-append(ParentCharVector,ConjugateCharVector)
  ConjRep<-as.numeric(length(conjugates)) # vector to create ConjugateMatrix
  ConjugateMatrix<-as.matrix(apply(MetaboliteData,2,rep,times=ConjRep))#conjugatematrix
  ConjugateMatrix[,2]<-as.character(paste(ConjugateMatrix[,2],ConjugateLabelVector,sep=""))
  
  
  ConjugateMultVector<-rep(Conjugate_mono_masses, each=length(ParentIDno)) # create vector for conjugate matrix multiplication
  ConjugateMass<-as.numeric(ConjugateMatrix[,3])#add conjugate masses
  ConjugateMass<-ConjugateMass+ConjugateMultVector
  ConjugateMatrix[,3]<-ConjugateMass
  MetaboliteData<-rbind(MetaboliteData,ConjugateMatrix) # row bind conjugate calc with original matrix
  
  MetabIDno<-seq(1,length(MetaboliteData[,1]),1) # sequence for unique metabolite ID
  MetaboliteData<-cbind(MetabIDno,Conjugate_type,MetaboliteData)# concat unique ID number
  
  colnames(MetaboliteData)[5]<-"Monoisotopic_mass"
  
  MetaboliteData<-aggregate.data.frame(MetaboliteData[-6],by=list(MetaboliteData$Monoisotopic_mass),paste,collapse=";")
  
  colnames(MetaboliteData)[1]<-"Monoisotopic_mass"
  
  
  } else if (conjugates=="NO"){
    ParentCharVector<-rep("Parent",each=length(ParentIDno)) # Parent names
    MetabIDno<-seq(1,length(MetaboliteData[,1]),1) # sequence for unique metabolite ID
    MetaboliteData<-cbind(MetabIDno,ParentCharVector,MetaboliteData)# concat unique ID number
  
    colnames(MetaboliteData)[2]<-"Conjugate_type"
        
    colnames(MetaboliteData)[5]<-"Monoisotopic_mass"
    
    MetaboliteData<-aggregate.data.frame(MetaboliteData[-6],by=list(MetaboliteData$Monoisotopic_mass),paste,collapse=";")
    
    colnames(MetaboliteData)[1]<-"Monoisotopic_mass"
  }
  
  
  # calculation of ESI adducts commonly found in negative and positive modes and generation of the search fields
    
  if (mode=="positive"){
    
      
    Monomass<-as.numeric(MetaboliteData[,1])
    M_H=c(Monomass+1.007276)  
    M_3H  =c(Monomass/3 +1.007276)  
    M_2H_Na =c(Monomass/3 + 8.334590 )  
    M_H_2Na =c(Monomass/3 +15.7661904)	
    M_3Na =c(Monomass/3 +22.989218)	
    M_2H=c(Monomass/2 +1.007276)	
    M_H_NH4=c(Monomass/2 +9.520550)	
    M_H_Na=c(Monomass/2 +11.998247)	
    M_H_K=c(Monomass/2 +19.985217)	
    M_2Na=c(Monomass/2 +22.989218)	
    M_NH4=c(Monomass +18.033823)
    M_Na=c(Monomass +22.989218)	
    M_CH3OH_H=c(Monomass +33.033489)
    M_K=c(Monomass +38.963158)	
    M_2Na_H=c(Monomass +44.971160)
    M_2K_H=c(Monomass +76.919040)	
    M2_H=c(2*Monomass +1.007276)	
    M2_NH4=c(2*Monomass +83.060373)	
    M2_Na=c(2*Monomass +22.989218)	
    M2_3H2O_2H=c(2*Monomass +28.02312)	
    M2_K=c(2*Monomass +38.963158)	
  
    DB<-data.frame(MetaboliteData, M_H=M_H,  M_3H  =M_3H  ,  M_2H_Na =M_2H_Na ,	M_H_2Na =M_H_2Na ,	M_3Na =M_3Na ,	M_2H=M_2H,	M_H_NH4=M_H_NH4,	M_H_Na=M_H_Na,	M_H_K=M_H_K,	M_2Na=M_2Na,	M_NH4=M_NH4,	M_Na=M_Na,	M_CH3OH_H=M_CH3OH_H,	M_K=M_K,	M_2Na_H=M_2Na_H,	M_2K_H=M_2K_H,	M2_H=M2_H,	M2_NH4=M2_NH4,	M2_Na=M2_Na,	M2_3H2O_2H=M2_3H2O_2H,	M2_K=M2_K)
    
     
    A<-21 #number of adduct columns
    B<-addColumns+6 #first adduct column for matrix subsetting
    C<-A+addColumns+5 #last column of new dataframe (4 =input DB file original database ID, unique identifier, name and parent mass)
    
    
  } else if (mode=="negative") {
    
    
    Monomass<-as.numeric(MetaboliteData[,1])
    M_H=c(Monomass-1.007276)  
    M_H2O_H=c(Monomass- 19.01839)  
    M_FA_H=c(Monomass + 44.998201 )	
    M_Na_2H=c(Monomass + 20.974666 )	
    M_K_2H=c(Monomass + 36.948606 )	
    M_Hac_H =c(Monomass + 59.013851 )	
    M2_FA_H=c(2*Monomass + 44.998201 )	
    M_2H=c(Monomass/2 - 1.007276 )	
    M_3H=c(Monomass/3 - 1.007276)	
    M2_H=c(2*Monomass - 1.007276 )	
    M3_H=c(3*Monomass - 1.007276 )	
    M2_Hac_H =c(2*Monomass + 59.013851 )
    
       
    DB<-data.frame(MetaboliteData,M_H=M_H,  M_H2O_H=M_H2O_H,  M_FA_H=M_FA_H,  M_Na_2H=M_Na_2H,	M_K_2H=M_K_2H,	M_Hac_H =M_Hac_H ,	M2_FA_H=M2_FA_H,	M_2H=M_2H,	M_3H=M_3H,	M2_H=M2_H,	M3_H=M3_H,	M2_Hac_H =M2_Hac_H)
    
    A<-12 #number of adduct columns
    B<-addColumns+6 #first adduct column for matrix subsetting
    C<-A+addColumns+5 #last column of new dataframe (4 =input DB file original database ID, unique identifier, name and parent mass)
  
    
  } 
   colnames(DB)[1]<-"Expected_ion"
 
  ##DATABASE CREATED##
  ################################################################################################################################
  DBMat <- as.matrix(DB[,B:C])
    
  ions <-sample[,2] 
  
  names <- as.character(DB[,5]) #compound names
  ParentIDno<-as.character(DB[,4]) #Parent ID number
  MetabIDno<-as.character(DB[,2])# Metabolite unique ID number
  
  adductnames<-colnames(DB)[B:C] #adduct names
  
  EICnames <- sample[,1] #return EIC numbers from Diff Report
  medRET<-sample[,3] #return median retention time
 
  message(paste("matching",nrow(sample),"MS1 features to database entries...PLEASE WAIT"))
  flush.console()
  
  pb<-txtProgressBar(min=0,max=length(ions),style=3)#,width=300)#title="Scatter plot progress bar"
  
  results <- data.frame() #empty dataframe for DB search results
  
  for(i in 1:length(ions)){
  
    Sys.sleep(0.1)
    setTxtProgressBar(pb,i)
    # setTkProgressBar(pb,i,label=paste( round(i/nrow(RawfoldCV)*100, 0),"% done"))
    flush.console()
        
    max = ions[i]+((MassAcc/1000000)*ions[i])
    min = ions[i]-((MassAcc/1000000)*ions[i])
    for (k in 1:A) {
      index <- which (DBMat[,k]<max & DBMat[,k]>min)
      
      
      if( length(index) > 0) {
        ExpMass<-DBMat[,k][index]
        metab_no<-MetabIDno[index]
        parent_no<-ParentIDno[index]
        ids <- names[index]
       
        
        identity <- cbind(EICnames[i],ions[i],medRET[i],adductnames[k],metab_no,parent_no,ids,ExpMass)#RSD_raw[i],RSD_corr[i],AverageIntQCs_raw[i],AverageIntQCs_corr[i])
        results <- rbind(results,identity)
      }
    }
  }
  
  
  message("...Done")
  flush.console()
    
  if(nrow(results)>1){
    
  r.df <- as.data.frame(results)
    
  r.dfcolnames<-c("MS1_feature_ID","m.z","medRT","adducts","MetabIDno","Parent_no","tent_assignment","ExpMass")

  colnames(r.df)<-r.dfcolnames 
   
  
  
  MetaboliteData<-as.data.frame(MetaboliteData)
  
  Metab_merge<-merge(r.df,MetaboliteData,sort=FALSE)
  
  Metab_merge<-as.matrix(Metab_merge)
  
  
  
  DeltaMass<-round(as.numeric(Metab_merge[,"m.z"])-as.numeric(Metab_merge[,"ExpMass"]),digits=4)
  
  
  ppm<-round((DeltaMass/as.numeric(Metab_merge[,"m.z"]))*1000000,digits=2)
  
  ppm<-ifelse(ppm<0,ppm*-1,ppm)#change negative ppm to positive
  
  ppm_belowMassAcc<-ifelse (MassAcc<ppm,0,1)#ppm below mass accuracy
  
  Index_ppm_below<-as.numeric(which(ppm_belowMassAcc==1)) #index of accurate assignments
  
  Metab_merge<-cbind(Metab_merge,DeltaMass,ppm)#data frame of search results
  
  aligned<-as.data.frame(Metab_merge[Index_ppm_below,])
  
  if (ncol(aligned)>1){
    
    aligned<-aggregate.data.frame(aligned,by=list(aligned$MS1_feature_ID,aligned$m.z,aligned$medRT),paste,collapse=";")
    
    aligned<-aligned[,-c(5,6,7)]
    
    colnames(aligned)[1]<-"MS1_feature_ID";colnames(aligned)[2]<-"m.z";colnames(aligned)[3]<-"medRT"
    
    index_belowdelta<-match(aligned[,"MS1_feature_ID"],sample[,1]) # sample index for results below delta
    
    aligned<-cbind(aligned,sample[index_belowdelta,])# extract rows from XCMS data below delta and create dataframe (much wider search than below ppm)
    
    unmatched<-sample[-index_belowdelta,]
    
    aligned<-subset(aligned, select=-c(MS1_feature_ID,MetabIDno,Monoisotopic_mass,ParentIDno,name))
    
    ###Place the EIC column name in the first column position
    aligned<-cbind(aligned[,EIC.column.name],aligned[,which(colnames(aligned)!=EIC.column.name)])
    
    write.csv (aligned, file="Allresults_belowMassTol.csv",row.names=FALSE)
    
    write.csv (unmatched,file="unAnnotated.csv",row.names=FALSE)
    
    
  }else if (ncol(aligned)<=1){
    message(paste("No Matches Below",MassAcc,"ppm mass accuracy"))
    flush.console()
  } 
  
  } else if (nrow(results)<1){
    message(paste("No Matches Below",MassAcc,"ppm mass accuracy"))
    flush.console()
  } 
  
}
###END###