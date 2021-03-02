GenesTest_AllCellLines<-function(x,y,z,mypval){
  coef1<<-z[[1]]
  coef2<<-z[[2]]
#Selecting samples
  Rapa_samples<-SampleInfo[which(SampleInfo$Rapa==x & SampleInfo$Diff_bin==y & (SampleInfo$Type%in%z)),] #Metadata
  Rapa_samples$Type<-factor(Rapa_samples$Type)
  Rapa_samples$CellLine<-factor(Rapa_samples$CellLine)
  Rapa_samples$Rapa<-factor(Rapa_samples$Rapa)
  Rapa_samples$Diff_bin<-factor(Rapa_samples$Diff_bin)
  
  
  Rapa_samples$SampleID<-rownames(Rapa_samples)
  RNA_Rapa<-Countdata_cpm[,Rapa_samples$SampleID] #RNAseq data

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

  param = SnowParam(4, "SOCK", progressbar=TRUE)
  register(param)

  # The variable to be tested must be a fixed effect
  form <- ~ Type + (1|CellLine) 

  # estimate weights using linear mixed model of dream
  vobjDream = voomWithDreamWeights(Norm, form, Rapa_samples)

  # Fit the dream model on each gene
  # By default, uses the Satterthwaite approximation for the hypothesis test
  fitmm = dream( vobjDream, form, Rapa_samples )
  #Sys.sleep(280)

  # get significant DE genes only (adjusted p-value < mypval)
  Output<-topTable(fitmm,coef="TypeNTC",n=dim(fit)[1])
  Results_limma_Filtered <- Output[which(Output$adj.P.Val<=mypval),]
  SignificantGeneList<-paste("SignificantGenes_",y,"_Rapa_",x,"_",paste(coef1,"-",coef2,sep=""),"RandomCellLine",".txt",sep="")
  write.table(Results_limma_Filtered ,file=SignificantGeneList,sep="\t",row.names = TRUE)
}
