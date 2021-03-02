#Calling libraries
require("ggpubr")
require("gplots")
require("devtools")

#Load the saved data
load("WGCNA_Analysis.RData")
load("CompleteTranscriptome_Analysis.RData")

#Calling WGCNA matrix
AllSamples_MEs<-MEs[which(rownames(MEs) %in% SampleInfo$SampleID),]
AllSamples_MEs$SampleID<-rownames(AllSamples_MEs)
AllSamples_MEs_v1<-merge(AllSamples_MEs,SampleInfo,by="SampleID")


ModuleNames<-names(AllSamples_MEs_v1[ , grepl("ME",names(AllSamples_MEs_v1))])

Plot_lists<-list()
for(i in 1:length(ModuleNames))
  {
  Module=ModuleNames[i]
  plots<-ggbarplot(AllSamples_MEs_v1, x = "CellLine",title=paste("AllSamples_GroupedByCellLine_",Module,sep=""),sort.by.groups=FALSE,y = Module,
            add = "mean_se", fill = "Type",palette = c("red","blue","green"),position = position_dodge(0.8))
  Plot_lists[[i]]<-plots
   }

for(i in 1:length(ModuleNames))
{
  Module=ModuleNames[i]
  filename=paste("AllSamples_GroupedByCellLine_",Module,sep="")
  png(paste(filename,".png"))
  print(Plot_lists[[i]])
  dev.off()
}

#Selecting samples
Rapa="n"
Diff="FALSE"


Rapa_samples<-SampleInfo[which(SampleInfo$Rapa==Rapa & SampleInfo$Diff_bin==Diff),] 
RNA_Rapa<-MEs[which(rownames(MEs) %in% Rapa_samples$SampleID),]
RNA_Rapa<-MEs[which(rownames(MEs) %in% Rapa_samples$SampleID),]
RNA_Rapa$SampleID<-rownames(RNA_Rapa)
RNA_Rapa_v1<-merge(RNA_Rapa,SampleInfo,by="SampleID")
RNA_Rapa<-RNA_Rapa_v1


BarPlot_lists<-list()
for(i in 1:length(ModuleNames))
{
  Module=ModuleNames[i]
  barplots<-ggbarplot(RNA_Rapa, x = "CellLine", y = Module, main=paste("Differentiation:",Diff,"_RAPA:",Rapa,"_GroupedByCellLine_",Module,sep=""),
          fill = "Type",add = "mean_se",color="Type",palette = c("red","blue","green"),position = position_dodge(0.9))
  BarPlot_lists[[i]]<-barplots
}

for(i in 1:length(ModuleNames))
{
  Module=colnames(RNA_Rapa[2:21])[i]
  filename=paste("Barplot_Diff",Diff,"_RAPA",Rapa,"_GroupedByCellLine_",Module,sep="")
  png(paste(filename,".png"))
  print(BarPlot_lists[[i]])
  dev.off()
}


############################################################
#
# Plotting 1000 Genes based on variance
#
############################################################


SampleInfo$DifferentiationStatus<-ifelse(SampleInfo$Diff_bin=="TRUE","DIFF","PROLIF")
SampleInfo$RapaStatus<-ifelse(SampleInfo$Rapa=="y","YRAP","NRAP")
SampleInfo$UniqIDs<-substr(SampleInfo$SampleID, start = 12, stop = 17)
SampleInfo$UniqueName<-paste(SampleInfo$UniqIDs,"_",SampleInfo$CellLine,"_",SampleInfo$DifferentiationStatus,"_",SampleInfo$RapaStatus,"_",SampleInfo$Type,sep="")
rownames(SampleInfo)<-paste(SampleInfo$UniqIDs,"_",SampleInfo$CellLine,"_",SampleInfo$DifferentiationStatus,"_",SampleInfo$RapaStatus,"_",SampleInfo$Type,sep="")


#Renaming the sample IDS for heatmap
for(i in 1:length(colnames(Countdata))){
  for(j in 1:length(SampleInfo$SampleID)){
    colnames(Countdata)[i]<-ifelse(substr(colnames(Countdata[i]),start = 12, stop = 17) == SampleInfo$UniqIDs[j],SampleInfo$UniqueName[j],colnames(Countdata)[i])
  }
}

Countdata$entrezgene<-rownames(Countdata)
Countdata_v<-merge(Countdata, GeneNames,by="entrezgene")
rownames(Countdata_v)<-make.unique(Countdata_v$hgnc_symbol,sep=".")



CountdataNames<-names(Countdata_v[ , grepl("11",names(Countdata_v))])
Countdata_v<-Countdata_v[,CountdataNames]

DGE_object<- DGEList(Countdata_v)

y<-DGE_object
group.col=c("green","red","yellow")[SampleInfo$Type]


logcounts <- cpm(y,log=TRUE)
var_genes <- apply(logcounts, 1, var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:1000]
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

Mcolors=sample(c("blue","red","green"), length(SampleInfo$CellLine), replace = TRUE, prob = NULL)
Ncolors=sample(c("green","red","yellow"), length(SampleInfo$Type), replace = TRUE, prob = NULL)
Rcolors=sample(c("purple","orange"), length(SampleInfo$Rapa), replace = TRUE, prob = NULL)
DColors=sample(c("grey","black"), length(SampleInfo$Diff_bin), replace = TRUE, prob = NULL)

clab=cbind(Mcolors,Ncolors,Rcolors,DColors)
colnames(clab)=c("CellLine","CellType","Rapa","Diff")

# Plot the heatmap
png("Heatmap3.png",width=1500, height=1500)
heatmap.3(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="1000 most variable genes",ColSideColors=clab,ColSideColorsSize=4,scale="row",margins=c(16,16))
legend("topright",legend=c("D244","D62","ReN","","2.1","2.2",
                           "NTC","","Rapa-Treated","Rapa-NonTreated","","Non-Diff.","Diff."),
       fill=c("blue","red","green","white","green","red","yellow","white","purple",
              "orange","white","grey","black"),  border=FALSE, bty="n", y.intersp = 0.7, cex=1)
dev.off()


pdf("Heatmap3.pdf",height=12,width = 12)
heatmap.3(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="1000 most variable genes",ColSideColors=clab,ColSideColorsSize=4,scale="row",margins=c(18,18))
legend("topright",legend=c("D244","D62","ReN","","2.1","2.2",
                           "NTC","","Rapa-Treated","Rapa-NonTreated","","Non-Diff.","Diff."),
       fill=c("blue","red","green","white","green","red","yellow","white","purple",
              "orange","white","grey","black"),  border=FALSE, bty="n", y.intersp = 0.7, cex=1)
dev.off()

###################################################
#
# Plotting highly expressed genes w.r.t each group
#
####################################################

Rapamycin=c("n","y")
Differentiation=c("FALSE","TRUE")
Type_In<-list(c("NTC","2.2"),c("NTC","2.1"))
CL=c("ReN","D62","D244")
mypval=0.05

for(r in 1:length(Rapamycin)){
  for(d in 1:length(Differentiation)){
    #for(Tp in 1:length(Type_In)){
      for(Cl in 1:length(CL))
        GenesTest_IndividualCellLines(Rapamycin[[r]],Differentiation[[d]],mypval,CL[[Cl]])
   # }
  }
}

#Selecting samples with Rapa==n, Differentiation=DIFF and Cell Line=D62

Plot_cellLine<-function(x,y,mypval,z){
  Rapa_samples_cellLine<-SampleInfo[which(SampleInfo$CellLine==z),]
  Rapa_samples<-Rapa_samples_cellLine[which(Rapa_samples_cellLine$Rapa==x & Rapa_samples_cellLine$Diff_bin==y),] #9
  RNA_Rapa<-MEs[which(rownames(MEs) %in% Rapa_samples$SampleID),]
  group<-factor(Rapa_samples$Type) #Cell type as group
  
  ################################## Define colours for heatmap based on Cell Type ############
  
  mypalette <- brewer.pal(11,"RdYlBu")
  morecols <- colorRampPalette(mypalette)
  # Set up colour vector for celltype variable
  col.cell <- c("purple","orange","red")[Rapa_samples$Type]
  
  ############################## Plotting highly expressed genes (counts) ########################
  
  select = order(rowMeans(RNA_Rapa), decreasing=TRUE)
  highexprgenes_counts<-RNA_Rapa[select,]
  highexprgenes_counts_t<-t(highexprgenes_counts)
  colnames(highexprgenes_counts_t)<- group
  
  filename<-paste("Heatmap_",y,"_Rapa_",x,"_",z,"_CellLine",sep="")
  png(paste(filename,".png",sep=""))
  heatmap.2(data.matrix(highexprgenes_counts_t),col=rev(morecols(50)),trace="none", main="Modules expression",ColSideColors=col.cell,scale="row")
  dev.off()
  
  
  ################################### Plotting variance of the genes ########################
  var_genes <- apply(RNA_Rapa, 2, var)
  head(var_genes)
  select_var <- names(sort(var_genes, decreasing=TRUE))
  
  head(select_var)
  # Subset logcounts matrix
  highly_variable <- data.frame(RNA_Rapa[,select_var])
  dim(highly_variable)
  
  filename<-paste("Arranged_Heatmap_",y,"_Rapa_",x,"_",z,"_CellLine",sep="")
  png(paste(filename,".png",sep=""))
  heatmap.2(data.matrix(t(highly_variable)),col=rev(morecols(50)),trace="none", main="Variable Modules across samples",ColSideColors=col.cell,scale="row")
  dev.off()
  
  
  ####################################### PCA plot ##############################################
  
  data_for_PCA <- highly_variable_lcpm
  dim(data_for_PCA)
  
  ## calculate MDS (matrix of dissimilarities)
  mds <- cmdscale(dist(data_for_PCA), k=3, eig=TRUE)  
  
  eig_pc <- mds$eig * 100 / sum(mds$eig)
  # plot the PCA
  filename<-paste("PCA_PropExplainedVariance_",y,"_Rapa_",x,"_",z,"_CellLine",sep="")
  png(paste(filename,".png",sep=""))
   barplot(eig_pc,
          las=1,
          xlab="Dimensions", 
          ylab="Proportion of explained variance (%)", y.axis=NULL,
          col="darkgrey")
  dev.off()
  
  
  ## calculate MDS
  mds <- cmdscale(dist(data_for_PCA)) # Performs MDS analysis 
  filename<-paste("PCA_Dim1vsDim2_",y,"_Rapa_",x,"_",z,"_CellLine",sep="")
  png(paste(filename,".png",sep=""))
  plot(mds[,1], -mds[,2], type="n", xlab="Dimension 1", ylab="Dimension 2", main="")
  text(mds[,1], -mds[,2], rownames(mds), cex=0.8) 
  dev.off()
}



####################################################################
#
# GO analysis
#
####################################################################


library(GOstats)

SignificantGenesfiles<-list.files(path = "/media/afsheenyousaf/My Passport/Mattson_Analysis/NewRawFiles/20201104_P2020-121-LEX/TranscriptomeAnalysis/GenesAnalysis/", pattern = ".txt$", full.names = TRUE)


MF<-list()

for(i in 1:length(SignificantGenesfiles)){
  entrez_genes<-list()
  MF[i]<-list(read.table(SignificantGenesfiles[i]))
  entrezgeneids <- as.character(rownames(MF[[i]]))
  length(entrezgeneids)
  universeids <- rownames(Countdata_filt)
  length(universeids)

hgCutoff <- 0.05
params <- new("GOHyperGParams",annotation="org.Hs.eg",geneIds=entrezgeneids,
              universeGeneIds=universeids,ontology="BP",pvalueCutoff=hgCutoff,conditional=TRUE,testDirection="over")
#  Run the test
hg <- hyperGTest(params)
# Get the p-values of the test
hg.pv <- pvalues(hg)
## Adjust p-values for multiple test (FDR)
hg.pv.fdr <- p.adjust(hg.pv,'fdr')
## select the GO terms with adjusted p-value less than the cut off
sigGO.ID <- names(hg.pv[hg.pv < hgCutoff])


# get table from HyperG test result
df <- summary(hg)
# keep only significant GO terms in the table
GOannot.table <- df[df[,1] %in% sigGO.ID,]
head(GOannot.table)

# Create text report of the significantly over-represented GO terms
ParName<-sub(".*/", "", SignificantGenesfiles[i])
NameID<-sub(".txt.*", "", ParName)
filename=paste("GOannot_",NameID,sep="")
write.table(GOannot.table,file=filename,sep="\t",row.names=F)


#Generating Biological pathway Plots

#GO plots
Genes<-entrezgeneids
UG<-universeids

pheno=paste("GOannot_",NameID,sep="")
#1) Gprofiler
gostres <- gost(query = as.vector(Genes), 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "custom", custom_bg = strtoi(UG), 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)


pdf(paste(pheno,"_.pdf",sep=""),height=4,width=5)
p<-gostplot(gostres, capped = FALSE, interactive = FALSE)
publish_gostplot(p,width = NA, height = NA, filename = NULL )
dev.off()



#2) GO Plots
GOoutput = gostres$result
if(identical(GOoutput$term_id,logical(0))==FALSE) {
  pdf(paste(pheno,"_GO_Plot.pdf",sep=""),height=4,width=5);
  par(oma=c(3,6,3,6));
  par(mgp=c(1, 0, 0))
  GO = GOoutput[(GOoutput$source=="GO:BP"|GOoutput$source=="GO:BP"|GOoutput$source=="GO:CC"),];
  if(length(GO$term_id)>0){
    bp = barplot(GO[order(GO$p_value,decreasing = FALSE),][,3][1:10],
                 main=paste("GO Ontology"),
                 horiz=TRUE,
                 yaxt='n',col="#D95F02",
                 xlab="p_value",
                 cex.main=0.5,cex.axis=0.3,cex.lab=0.5,tck=-0.05);
    axis(2,at=bp,labels=GO[order(GO$p_value,decreasing = FALSE),][,11][1:10],
         tick=FALSE,las=2,cex.axis=0.3);
    abline(v=0.05,col="red",lwd=2,lty=1);
    dev.off()
  }else {
    cat("Not enough enriched terms for GO plot")
  }
  
  #3) KEGG plot   
  pdf(paste(pheno,"KEGGPlot.pdf",sep=""),height=4,width=5);
  par(oma=c(3,6,3,6));
  par(mgp=c(1, 0, 0))
  kegg = GOoutput[(GOoutput$source=="KEGG"),];
  if(length(kegg$term_id)>0){
    KG = barplot(kegg[order(kegg$p_value,decreasing = FALSE),][,3][1:10],
                 main=paste("KEGG pathways"),
                 horiz=TRUE,
                 yaxt='n',col="#D95F02",
                 xlab="p_value",
                 cex.main=0.5,cex.axis=0.3,cex.lab=0.5,tck=-0.05);
    axis(2,at=KG,labels=kegg[order(kegg$p_value,decreasing = FALSE),][,11][1:10],
         tick=FALSE,las=2,cex.axis=0.3);
    abline(v=0.05,col="red",lwd=2,lty=1);
    dev.off()
  } else {
    cat("Not enough enriched terms for KEGG plot")
  }
}


}




