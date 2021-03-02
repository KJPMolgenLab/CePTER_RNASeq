###############################################
#					                                  	#
# 	01) Mapping And Alignment	              	#
#					                                    #
###############################################


################################
#				                       #
#	Libraries to call	           #
#				                       #
################################

package<-c("limma", "RColorBrewer","flashClust","textshape","Rtsne","ggpubr")
if(length(setdiff(package, rownames(installed.packages()))) > 0)	{
  install.packages(setdiff(package, rownames(installed.packages())),INSTALL_opts = c('--no-lock'),ask=FALSE)  
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("biomaRt","Rsubread","edgeR","WGCNA","DESeq2","variancePartition","BiocParallel"))


require("Rsubread")
require("biomaRt")
require("edgeR")
require("limma")
require("RColorBrewer")
require("WGCNA")
require("DESeq2")
require("flashClust")
require("variancePartition")
require("BiocParallel")
require("textshape")
require("Rtsne")
require("ggpubr")
require("GOstats")
require("gprofiler2")

################################
#				                       #
#	Sourcing the functions       #
#				                       #
################################


source("TranscriptomeAnalysis/Codes/ModuleTest.R") # Function to identify modules differentially regulated
################################
#				                       #
#	Preprocessing & QC           #
#				                       #
################################

#QC was performed in FastQC 
#Filtering was performed in trimmomatic (trimmomatic_processing.sh)

#Reference base was downloaded from UCSC
#http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/

#Build Index
#buildindex(basename="Analysis_Jones",reference="/home/afsheenyousaf/bowtie_output/HumanRefGenome/hg38.fa.gz")

#Read in the fastq files
#fastq.files<-list.files(path="/media/afsheenyousaf/My Passport/Mattson_Analysis/NewRawFiles/20201104_P2020-121-LEX/FastqFiles_Extracted",pattern=".fq$",full.names=TRUE)

#Aligning the reads
#align(index="Analysis_Jones",readfile1=fastq.files)

#Reading in Bam files
#bam.files <- list.files(path = "/media/afsheenyousaf/My Passport/Mattson_Analysis/NewRawFiles/20201104_P2020-121-LEX/FastqFiles_Extracted", pattern = ".BAM$", full.names = TRUE)

#Getting the feature counts
#fc <- featureCounts(bam.files, annot.inbuilt="hg38",countMultiMappingReads=TRUE)

###############################################################################################################################################
###############################################################################################################################################
#                                                                                                                                             #
# Fastqc files were read in, mapped, aligned and counts were called. All these objects are stored in  "01_MappingAndAlignment_trim.RData"     #
#                                                                                                                                             #
###############################################################################################################################################
###############################################################################################################################################

#loading the data
load("01_MappingAndAlignment_trim.RData")

#Deleting extra information from column names
Countdata <- data.frame(fc$counts)
colnames(Countdata)<-sub("_.*", "", colnames(Countdata))
write.table(Countdata,file="CountData_UnAnnotated_Unfiltered.txt",sep="\t", quote=F,append=F)

################
#  SampleInfo  #
################

#Creating factors for groups
SampleInfo<-read.table("/media/afsheenyousaf/My Passport/Mattson_Analysis/TranscriptomeAnalysis_pheno.txt",sep="\t",header=TRUE)
rownames(SampleInfo)<-SampleInfo$SampleID
SampleInfo$Differentiation.1<-as.factor(SampleInfo$Differentiation.1)
SampleInfo$Diff_bin<-as.factor(ifelse(grepl("Diff", SampleInfo$Differentiation.1), "TRUE", "FALSE"))
SampleInfo$KO_bin<-as.factor(ifelse(grepl("NTC", SampleInfo$Type), "FALSE", "TRUE"))
SampleInfo$CellLine<-as.factor(SampleInfo$CellLine)
SampleInfo$Type<-as.factor(SampleInfo$Type)
SampleInfo$Rapa<-as.factor(SampleInfo$Rapa)
write.table(SampleInfo,"SampleInfo_Final.txt",sep="\t",row.names=FALSE,quote=FALSE)

SampleInfo <- SampleInfo %>% column_to_rownames("SampleID")
nSampleInfo = nrow(SampleInfo)
nSampleData = ncol(Countdata)
nSampleGenes = nrow(Countdata)
Samples = rownames(SampleInfo)[rownames(SampleInfo) %in% colnames(Countdata)]
noverlap = length(Samples)
Countdata = Countdata[,Samples]
SampleInfo = SampleInfo[Samples,]
nReads = colSums(Countdata)
avreads = round(mean(nReads),3)
sdreads = round(sd(nReads),3)



#Quality Check Plots

#MDS plot
vars=apply(Countdata,1,var)
index = order(vars, decreasing = T)[]
forcluster=log(t(Countdata[index,]+1))
d = dist(forcluster)
fit = cmdscale(d,eig=TRUE, k=5)
pairplot <- function(coltoshow){
  pairs(fit$points, pch=16, col=SampleInfo[,coltoshow], main = paste0("PCA ", coltoshow))
  labels=levels(SampleInfo[,coltoshow])
  legend("topleft", legend = labels, col=c(1:length(labels)), xpd=T, pch=16, bty = "n", cex=0.7)
}
pairplot("Type")
pairplot("CellLine")
pairplot("Diff_bin")
pairplot("KO_bin")
pairplot("Rapa")

#TSNE plot
tsne_model_1 = Rtsne(as.matrix(forcluster), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=3)
tsnepairplot <- function(coltoshow){
  pairs(tsne_model_1$Y, pch=16, col=SampleInfo[,coltoshow], main = paste0("tsne ", coltoshow))
  labels=levels(SampleInfo[,coltoshow])
  legend("topleft", legend = labels, col=c(1:length(labels)), xpd=T, pch=16, bty = "n", cex=0.7)
}
tsnepairplot("Type")
tsnepairplot("CellLine")
tsnepairplot("Diff_bin")
tsnepairplot("KO_bin")
tsnepairplot("Rapa")




###############################################
#						                                  #
# 		02) WGCNA			                          #
#						                                  #
###############################################



colData=data.frame(row.names = colnames(Countdata))
coldata <- data.frame(SampleInfo[,c("Type","CellLine","Rapa","Diff_bin","KO_bin")])
rownames(coldata)<-rownames(SampleInfo)
all(rownames(coldata) %in% colnames(Countdata))# TRUE
all(rownames(coldata) == colnames(Countdata))#FALSE
Countdata <- Countdata[, rownames(coldata)]
all(rownames(coldata) == colnames(Countdata))

#################
#   Filtering	 #
#################

#Keep genes with least 5 count-per-million reads (cpm) in at least 3 samples

Expression <- rowSums(cpm(Countdata)>5) >= 3 #Since we hv three replicates
table(Expression)
#FALSE  TRUE 
#14093 14302 


#Converting the counts to logcounts
Countdata_filt<-Countdata[Expression,]
Countdata_log<-cpm(Countdata_filt,log=TRUE)

Dataset_pre= DESeqDataSetFromMatrix(countData = Countdata_filt,colData = coldata,design=~ Type + CellLine + Diff_bin + Rapa + Type*Rapa)


#Estimate size factors and dispersions
dds2 = DESeq(Dataset_pre)

#Perform normalization from the fitted dispersion-mean relation(s) and then transforms the count data (normalized by division by the size factor), suggested by WGNA manuals
vsd = getVarianceStabilizedData(dds2)

WGCNA_matrix <- t(vsd) #Need to transform for further calculations


s = abs(bicor(WGCNA_matrix)) #biweight mid-correlation
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(WGCNA_matrix, powerVector = powers, verbose = 5)

png("Softpower.png")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence'));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col='red'); abline(h=0.90,col='red')
dev.off()


png("ScaleFreeToplogy.png")
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.80,col="red")
dev.off()



#Mean connectivity
png("MeanConnectivity.png")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

softPower = 6; # The point where the curve flattens
#calclute the adjacency matrix
adj= adjacency(WGCNA_matrix,type = "unsigned", power = softPower);

#Converting adjacency matrix into To so that the noice could be reduced
TOM=TOMsimilarityFromExpr(WGCNA_matrix,networkType = "unsigned", TOMType = "unsigned", power = softPower);

SubGeneNames<-colnames(WGCNA_matrix)

colnames(TOM) =rownames(TOM) =SubGeneNames
dissTOM=1-TOM

#hierarchical clustering
geneTree = flashClust(as.dist(dissTOM),method="average");

#plot the resulting clustering tree (dendrogram)
png("ClusteringTree.png")
plot(geneTree, xlab="", sub="",cex=0.3);
dev.off()

#Set the minimum module size
minModuleSize = 30;

#Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,method="tree", minClusterSize = minModuleSize);

#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)

#Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

png("Dendro.png")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

#set the diagonal of the dissimilarity to NA 
diag(dissTOM) = NA;

#Visualize the Tom plot. Raise the dissimilarity matrix to the power of 4 to bring out the module structure
png("TOMplot.png")
TOMplot(dissTOM, geneTree, as.character(dynamicColors))
dev.off()


for (color in dynamicColors){
  module=SubGeneNames[which(dynamicColors==color)]
  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}

module.order <- unlist(tapply(1:ncol(WGCNA_matrix),as.factor(dynamicColors),I))
m<-t(t(WGCNA_matrix[,module.order])/apply(WGCNA_matrix[,module.order],2,max))

png("Moduleplot.png")
heatmap(t(m),zlim=c(0,1),col=gray.colors(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",RowSideColors=dynamicColors[module.order])
dev.off()

#calculate eigengenes
MEList = moduleEigengenes(WGCNA_matrix, colors = dynamicColors)
MEs = MEList$eigengenes 

png("EigenGeneNetwork")
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

# Call an automatic merging function
merge = mergeCloseModules(WGCNA_matrix, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
png("Plots_before_afterMerging.png")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
#Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs; # We selected these MEs

# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file="WGCNAAnalysis_merged.RData")

for (color in moduleColors){
  module=SubGeneNames[which(moduleColors==color)]
  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}


save.image(file = "WGCNA_Analysis.RData")



#Fetch module genes

ensembl_version = "https://dec2016.archive.ensembl.org"

Modulefiles <- list.files(path = "/media/afsheenyousaf/My Passport/Mattson_Analysis/NewRawFiles/20201104_P2020-121-LEX/TranscriptomeAnalysis/WGCNA/Merged_Modules/Entrez_Modules/", pattern = ".txt$", full.names = TRUE)


#Name the modules
#paste("Module",Module[i])<-read.table(Modulefiles[1])
MF<-list()

for(i in 1:length(Modulefiles)){
  entrez_genes<-list()
  MF[i]<-list(read.table(Modulefiles[i]))
  x<-sub(".*_", "", Modulefiles[i]) 
  y<-sub(".txt.*", "", x)
  colnames(MF[[i]])<-paste("Module_",y,sep ="")
  
  # get the Ensembl annotation for human genome
  mart<- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL",host=ensembl_version))
  
  # get entrez gene IDs from limma output table
  entrez_genes[[i]] <- MF[[i]]
  length(entrez_genes[[i]][,1])
  
  # interrogate the BioMart database to get gene symbol and description for these genes 
  GeneNames<- getBM(
    filters= "entrezgene",
    attributes= c("entrezgene","hgnc_symbol","description"),
    values= entrez_genes[[i]],
    mart= mart)
  
  #Name the genes
  filename = paste("GeneNames_",colnames(MF[[i]]),".txt",sep="")
  write.table(GeneNames,file=filename,row.names=F,quote=F,sep="\t")
  
}


###############################################
#					                                  	#
#     02) Modules up/down regulated upon KO	  #
#    of DEPDC5					                      #
#						                                  #
###############################################

source("TranscriptomeAnalysis/Codes/ModuleTest.R")

Rapamycin=c("n","y")
Differentiation=c("FALSE","TRUE")
Type_In<-list(c("NTC","2.2"),c("NTC","2.1"))
mypval=0.05

for(r in 1:length(Rapamycin)){
  for(d in 1:length(Differentiation)){
    for(Tp in 1:length(Type_In)){
            ModuleTest(Rapamycin[[r]],Differentiation[[d]],Type_In[[Tp]],mypval)
                                }
  }
}

#################################################
#
# Modules DE w.r.t individual cell lines
#
#################################################


source("TranscriptomeAnalysis/Codes/ModuleTest_IndividualCellLines.R")
Rapamycin=c("n","y")
Differentiation=c("FALSE","TRUE")
Type_In<-list(c("NTC","2.2"),c("NTC","2.1"))
CL=c("ReN","D62","D244")
mypval=0.05

for(r in 1:length(Rapamycin)){
  for(d in 1:length(Differentiation)){
    for(Tp in 1:length(Type_In)){
      for(Cl in 1:length(CL))
        ModuleTest_IndividualCellLines(Rapamycin[[r]],Differentiation[[d]],Type_In[[Tp]],mypval,CL[[Cl]])
    }
  }
}



##########################################################################################################################################################################################################
#
#                                                             Gene Analysis
#
##########################################################################################################################################################################################################


#############################
#                           #
#   Mapping to HGNC symbols #
#                           #
#############################

Countdata_Annotated<-Countdata_cpm

Countdata_Annotated$entrezgene<-rownames(Countdata_cpm)

ensembl_version = "https://dec2016.archive.ensembl.org"

# get the Ensembl annotation for human genome
mart<- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL",host=ensembl_version))

# get entrez gene IDs
entrez_genes <- as.character(rownames(Countdata_Annotated))
length(entrez_genes)

# interrogate the BioMart database to get gene symbol and description for these genes 
GeneNames <- getBM(filters= "entrezgene",attributes= c("entrezgene","hgnc_symbol","description"),values= entrez_genes,mart= mart)

Countdata_Genes<-merge(Countdata_Annotated, GeneNames,by="entrezgene",all.x=TRUE)
rownames(Countdata_Genes)<-make.names(Countdata_Genes$hgnc_symbol,unique = TRUE)
write.table(Countdata_Genes,file="CountData_Annotated_Filtered.txt",sep="\t", quote=F,append=F)



###############################################
#						                                  #
# 01) Genes up/down regulated, Cell Line      #
#     as random effect	                      #
#                                     	      #
#						                                  #
###############################################



source("TranscriptomeAnalysis/Codes/GenesTest_AllCellLines.R")
Rapamycin=c("n","y")
Differentiation=c("FALSE","TRUE")
Type_In<-list(c("NTC","2.2"),c("NTC","2.1"))
mypval=0.05

for(r in 1:length(Rapamycin)){
  for(d in 1:length(Differentiation)){
    for(Tp in 1:length(Type_In)){
      GenesTest_AllCellLines(Rapamycin[[r]],Differentiation[[d]],Type_In[[Tp]],mypval)
    }
  }
}



################################################
#
# Performing in each cell type
#
#################################################
Rapa="n"
Diff="TRUE"
CoefToTest="NTC-2.2"
Type_Inc<-c("NTC","2.2")
mypval=0.05


source("TranscriptomeAnalysis/Codes/GenesTest_IndividualCellLines.R")
Rapamycin=c("n","y")
Differentiation=c("FALSE","TRUE")
Type_In<-list(c("NTC","2.2"),c("NTC","2.1"))
CL=c("ReN","D62","D244")
mypval=0.05

for(r in 1:length(Rapamycin)){
  for(d in 1:length(Differentiation)){
    for(Tp in 1:length(Type_In)){
      for(Cl in 1:length(CL))
      GenesTest_IndividualCellLines(Rapamycin[[r]],Differentiation[[d]],Type_In[[Tp]],mypval,CL[[Cl]])
    }
  }
}



save.image(file = "CompleteTranscriptome_Analysis.RData")

