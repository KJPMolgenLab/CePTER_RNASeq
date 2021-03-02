#GenesTest_IndividualCellLines<-function(Rapa,Diff,CoefToTest,Type_Inc,mypval,CL){
GenesTest_IndividualCellLines<-function(x,y,z,mypval,w){
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
  RNA_Rapa<-Countdata_cpm[,Rapa_samples$SampleID] #RNAseq data
  
  RNA_Rapa<-Countdata_cpm[,which(colnames(Countdata_cpm) %in% Rapa_samples$SampleID)] #RNAs

  group<-factor(Rapa_samples$Type) #Cell type as group
  CellLine<-factor(Rapa_samples$CellLine)
  
#############################
#
# Data Normalization
#
############################

  #Converting counts to DGElist
  DGE_object<- DGEList(RNA_Rapa) #Also does not accept log counts
  #Preprocessing
  Norm <- calcNormFactors(DGE_object)
  dim(Norm)

###############
# Model
###############
################################
# Creating design matrix
################################


  design = model.matrix( ~Type, Rapa_samples)
  fit <- lmFit(RNA_Rapa,design)
  fit <- eBayes(fit)
  options(digits=3)

  Output<-topTable(fit,coef="TypeNTC",n=dim(fit)[1])
  Results_limma_Filtered <- Output[which(Output$adj.P.Val<=mypval),]
  Results_limma_Filtered$entrezgene<-rownames(Results_limma_Filtered)
  Results_limma_Filtered_v1<-merge(GeneNames,Results_limma_Filtered,by="entrezgene")
  SignificantGeneList<-paste("SignificantGenes_",y,"_Rapa_",x,"_",paste(coef1,"-",coef2,sep=""),w,"_CellLine",".txt",sep="")
  write.table(Results_limma_Filtered_v1 ,file=SignificantGeneList,sep="\t",row.names = TRUE)
}
