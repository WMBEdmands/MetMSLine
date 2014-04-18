DBAnnotate<-function(X="Coffee_24hrDR_DM_results.csv",database="coffee_FooDB_filtered.csv",mode="positive",conjugates=c("Gluc","Sulf","DiGluc","DiSulf","GlucSulf","NAcCys","Glyc"),MassAcc=10,wd="Z:\\Raw_data\\EPIC_POS_Ur_n481_04_02_14\\Auto.MV.Discriminant.results\\24hrDR\\Coffee_24hrDR_DM\\",unknowns.dir="Z:\\Raw_data\\EPIC_POS_Ur_n481_04_02_14\\Auto.MV.Discriminant.results\\24hrDR\\Coffee_24hrDR_DM\\"){
  # input XCMS diff report and tentative metabolite list with monoisotopic masses, also require
  # ionisation mode, anticipated acceptable mass accuracy (dependent on mass spectrometer performance) and retention time window for adduct determination returns txt file in working
  setwd(unknowns.dir)
  
  sample<-read.csv(X,header=T)
  #sample<-sample[,-1]
  setwd(wd)
  MetaboliteData<-data.frame(read.csv(database,header=T))
  setwd(unknowns.dir)
  
  addColumns<-length(MetaboliteData)-2
   
  ParentIDno<-seq(1,length(MetaboliteData[,1]),1) #Parent unique ID number
  MetaboliteData<-cbind(ParentIDno,MetaboliteData) #Parent ID number concatenated before conjugate and adduct calculation
  
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
  Glyc=75.032029
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
  
  # calculation of ESI adducts commonly found in negative and positive modes and generation of the search fields
    
  if (mode=="positive"){
    
      
    Monomass<-as.numeric(MetaboliteData[,1])
    M_H=c(Monomass+1.007276)  
    
  
    DB<-data.frame(MetaboliteData, M_H=M_H)
    
     
    A<-1 #number of adduct columns
    B<-addColumns+6 #first adduct column for matrix subsetting
    C<-A+addColumns+5 #last column of new dataframe (4 =input DB file original database ID, unique identifier, name and parent mass)
    
    
  } else if (mode=="negative") {
    
    
    Monomass<-as.numeric(MetaboliteData[,1])
    M_H=c(Monomass-1.007276)  
       
       
    DB<-data.frame(MetaboliteData,M_H=M_H)
    
    A<-1 #number of adduct columns
    B<-addColumns+6 #first adduct column for matrix subsetting
    C<-A+addColumns+5 #last column of new dataframe (4 =input DB file original database ID, unique identifier, name and parent mass)
  
    
  } 
   colnames(DB)[1]<-"Expected_ion"
 
  ##DATABASE CREATED##
  ################################################################################################################################
  DBMat <- as.matrix(DB[,B:C])
    
  ions <-sample[,"mzmed"] 
  
  names <- as.character(DB[,5]) #compound names
  ParentIDno<-as.character(DB[,4]) #Parent ID number
  MetabIDno<-as.character(DB[,2])# Metabolite unique ID number
  
  adductnames<-colnames(DB)[B:C] #adduct names
  
  EICnames <- sample[,"XCMS_EIC"] #return EIC numbers from Diff Report
  medRET<-sample[,"rtmed"] #return median retention time
 
  
  results <- data.frame() #empty dataframe for DB search results
  
  for(i in 1:length(ions)){
    max = ions[i]+0.1
    min = ions[i]-0.1
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
  
  if(nrow(results)>1){
    
    r.df <- as.data.frame(results)
  
  
  
  r.dfcolnames<-c("XCMS_EIC","XCMSmzmed","medRT","adducts","MetabIDno","Parent_no","tent_assignment","ExpMass")

  colnames(r.df)<-r.dfcolnames 
   
  
  
  MetaboliteData<-as.data.frame(MetaboliteData)
  
  Metab_merge<-merge(r.df,MetaboliteData,sort=FALSE)
  
  Metab_merge<-as.matrix(Metab_merge)
  
  
  
  DeltaMass<-round(as.numeric(Metab_merge[,"XCMSmzmed"])-as.numeric(Metab_merge[,"ExpMass"]),digits=4)
  
  
  ppm<-round((DeltaMass/as.numeric(Metab_merge[,"XCMSmzmed"]))*1000000,digits=2)
  
  ppm<-ifelse(ppm<0,ppm*-1,ppm)#change negative ppm to positive
  
  ppm_belowMassAcc<-ifelse (MassAcc<ppm,0,1)#ppm below mass accuracy
  
  Index_ppm_below<-as.numeric(which(ppm_belowMassAcc==1)) #index of accurate assignments
  
  Metab_merge<-cbind(Metab_merge,DeltaMass,ppm)#data frame of search results
  
  XCMSaligned<-as.data.frame(Metab_merge[Index_ppm_below,])
  
  if (ncol(XCMSaligned)>1){
    
  XCMSaligned<-aggregate.data.frame(XCMSaligned,by=list(XCMSaligned$XCMS_EIC,XCMSaligned$XCMSmzmed,XCMSaligned$medRT),paste,collapse=";")
  
  XCMSaligned<-XCMSaligned[,-c(5,6,7)]
  
  colnames(XCMSaligned)[1]<-"XCMS_EIC";colnames(XCMSaligned)[2]<-"XCMSmzmed";colnames(XCMSaligned)[3]<-"medRT"
  
  XCMSindex_belowdelta<-match(XCMSaligned[,"XCMS_EIC"],sample[,"XCMS_EIC"]) # sample index for results below delta
  
  XCMSaligned<-cbind(XCMSaligned,sample[XCMSindex_belowdelta,])# extract rows from XCMS data below delta and create dataframe (much wider search than below ppm)
  
  unmatched<-sample[-XCMSindex_belowdelta,]
 
  XCMSaligned<-subset(XCMSaligned, select=-c(XCMS_EIC,MetabIDno,Monoisotopic_mass,ParentIDno,name))
  
  write.csv (XCMSaligned, file="Allresults_belowMassTol.csv",row.names=FALSE)

  write.csv (unmatched,file="unAnnotated.csv",row.names=FALSE)
  
  
  }else if (ncol(XCMSaligned)<=1){
    print("No Matches Below ppm value")
  } 
  
  } else if (nrow(results)<1){
    print("No Matches Below delta 0.01 mass difference")
    
   } 

}

###END###