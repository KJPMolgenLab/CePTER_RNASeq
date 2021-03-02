ModuleTest<-function(x,y,z,mypval){
  coef1<<-z[[1]]
  coef2<<-z[[2]]
   #Selecting samples
  Rapa_samples<-SampleInfo[which(SampleInfo$Rapa==x & SampleInfo$Diff_bin==y & (SampleInfo$Type%in%z)),] #Metadata
  Rapa_samples$Type<-factor(Rapa_samples$Type)
  Rapa_samples$CellLine<-factor(Rapa_samples$CellLine)
  Rapa_samples$Rapa<-factor(Rapa_samples$Rapa)
  Rapa_samples$Diff_bin<-factor(Rapa_samples$Diff_bin)
  
  Rapa_samples$SampleID<-rownames(Rapa_samples)
  RNA_Rapa<-MEs[which(rownames(MEs) %in% Rapa_samples$SampleID),]
  
  group<-factor(Rapa_samples$Type) #Cell type as group
  Cell_Line<-factor(Rapa_samples$CellLine)
  
  
  #####################################################################################
  ####################### (1) Cell Types as groups ####################################
  #####################################################################################
  
  ################################
  # Creating design matrix
  ################################
  design = model.matrix( ~0 + Type, Rapa_samples)
  colnames(design)<- gsub("Type","",colnames(design))
  
  
  ################################ 
  #using lmfit with random effect
  ################################
  
  dupcor <- duplicateCorrelation(t(RNA_Rapa), design, block=Rapa_samples$CellLine)
  
  #This model can handle single random effect and allows its magnitude to be the same across modules/genes.
  fitDupCor <- lmFit(t(RNA_Rapa), design, block=Rapa_samples$CellLine, correlation=dupcor$consensus)
  colnames(fitDupCor$design)<- gsub("Type","",colnames(fitDupCor$design))
  colnames(fitDupCor$coefficients)<- gsub("Type","",colnames(fitDupCor$coefficients))
  
  
  #construct the contrast matrix 
  cont.matrix<-makeContrasts(paste(coef1,"-",coef2,sep=""),levels=make.names(colnames(fitDupCor$design)))
  cont.matrix 
  fit_v1<- contrasts.fit(fitDupCor, cont.matrix)
  fit<- eBayes(fit_v1)
  
  
  
  #get significant DE genes only (adjusted p-value < mypval)
  Output<-topTable(fit,coef=paste(coef1,"-",coef2,sep=""),n=dim(fit)[1])
  Results_limma_Filtered <- Output[which(Output$adj.P.Val<=mypval),]
  SignificantModList<-paste("SignificantModules_",y,"_Rapa_",x,"_",paste(coef1,"-",coef2,sep=""),"RandomCellLine",".txt",sep="")
  write.table(Results_limma_Filtered,file=SignificantModList,sep="\t")
}