ModuleTest_IndividualCellLines<-function(x,y,z,mypval,w){
  coef1<<-z[[1]]
  coef2<<-z[[2]]
  #Selecting samples
  Rapa_samples_CL<-SampleInfo[which(SampleInfo$CellLine==w),] #9
  Rapa_samples<-Rapa_samples_CL[which(Rapa_samples_CL$Rapa==x & Rapa_samples_CL$Diff_bin==y  & (Rapa_samples_CL$Type%in%z)),] #9#Metadata
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
  #using lmfit
  ################################
  
  design = model.matrix( ~Type, Rapa_samples)
  fit <- lmFit(t(RNA_Rapa),design)
  fit <- eBayes(fit)
  options(digits=3)
  
  
  #get significant DE genes only (adjusted p-value < mypval)
  Output<-topTable(fit,coef="TypeNTC",n=dim(fit)[1])
  Results_limma_Filtered <- Output[which(Output$adj.P.Val<=mypval),]
  SignificantModList<-paste("SignificantModules_",y,"_Rapa_",x,"_",paste(coef1,"-",coef2,sep=""),w,"CellLine",".txt",sep="")
  write.table(Results_limma_Filtered,file=SignificantModList,sep="\t")
}