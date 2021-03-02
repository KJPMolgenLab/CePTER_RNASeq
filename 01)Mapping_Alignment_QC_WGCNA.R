###############################################
#						#
# 	01) Mapping And Alignment		#
#						#
###############################################

fastqcfiles = "M:/afsheen/Mattson_TranscriptomeAnalysis/FastQC_results/"
refhg38 = "M:/chiocchetti/ReferenceGenomes/hg38.fa.gz"

################################
#				#
#	Libraries to call	#
#				#
################################

package<-c("edgeR","limma", "RColorBrewer")

if(length(setdiff(package, rownames(installed.packages()))) > 0)	{
	install.packages(setdiff(package, rownames(installed.packages())),INSTALL_opts = c('--no-lock'),ask=FALSE)  
                                                				}
                                                              										
if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")
	BiocManager::install(c("biomaRt", "Rsubread"))
	
require("Rsubread")
require("biomaRt")
require("edgeR")
require("limma")
require("RColorBrewer")
require("WGCNA")
require("DESeq2")
require("flashClust")
################################
#				#
#	Preprocessing & QC	#
#				#
################################

#QC was performed in FastQC 
#Filtering was performed in trimmomatic


#Read filter
#for filename in *.fastq
#do
        # first, make the base by removing fastq.gz
#        base=$(basename $filename .fastq)
#        echo $base

#       trimmomatic SE ${base}.fastq \
#                ${base}.qc.fq \
#                ILLUMINACLIP:/home/afsheenyousaf/anaconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/TruSeq3-SE.fa:2:0:15 \
#                LEADING:3 TRAILING:3 \
#                SLIDINGWINDOW:4:2 \
#                MINLEN:36
#done              

#Reference base was downloaded from UCSC

#http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/


#Build Index
buildindex(basename="Analysis_Jones",reference=refhg38)

#Read in the fastq files
fastq.files<-list.files(path="/media/afsheenyousaf/My Passport/Mattson_Analysis/NewRawFiles/20201104_P2020-121-LEX/FastqFiles_Extracted",pattern=".fq$",full.names=TRUE)


#Aligning the reads
align(index="Analysis_Jones",readfile1=fastq.files)


#Reading in Bam files
bam.files <- list.files(path = "/media/afsheenyousaf/My Passport/Mattson_Analysis/NewRawFiles/20201104_P2020-121-LEX/FastqFiles_Extracted", pattern = ".BAM$", full.names = TRUE)


#Getting the feature counts
fc <- featureCounts(bam.files, annot.inbuilt="hg38",countMultiMappingReads=TRUE)


#Deleting extra information from column names
Countdata <- data.frame(fc$counts)
colnames(Countdata)<-sub("_.*", "", colnames(Countdata))





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



#################
#   Filtering	 #
#################

#Keep genes with least 1 count-per-million reads (cpm) in at least 3 samples


Expression <- rowSums(Countdata>5) >= 3 #Since we hv three replicates

table(Expression)
#FALSE  TRUE 
#10037 18358 



DGE_object<- DGEList(Countdata[Expression,])

logcounts <- cpm(DGE_object,log=TRUE)

#Quality Check Plots

#library Sizes
png("Plot_LibrarySizes.png")
barplot(DGE_object$samples$lib.size,names=colnames(DGE_object),las=2)
title("Barplot of library sizes")
dev.off()


#Check distribution of samples
png("SampleDistribution.png")
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
dev.off()


#MDS plot

#levels(SampleInfo$Type)
col.Type <- c("green","red","yellow")[SampleInfo$Type]
pch.celll<-c(21,22,23)[SampleInfo$CellLine]
col.Diff<- c(4,11)[SampleInfo$Diff_bin]

pchl<-c(21,22,23) #cell line
pchD<-c(4,11)# differentiation




png("MDSplot_main_Type_Diff.png")
plotMDS(DGE_object,pch=col.Diff,col=col.Type,cex=1.5)
legend("topleft",fill=c("green","red","yellow"),legend=levels(SampleInfo$Type))
legend("bottomleft",c("FALSE","TRUE"),pch=pchD, legend=c("Non-Diff","Diff"))
dev.off()



png("MDSplot_main_Line_Diff.png")
col.CL<-c("blue","orange","violet")[SampleInfo$Type]
plotMDS(DGE_object,pch=col.Diff,col=col.CL,cex=1.5)
legend("topleft",fill=c("blue","orange","violet"),legend=levels(SampleInfo$CellLine))
legend("bottomleft",c("FALSE","TRUE"),pch=pchD, legend=c("Non-Diff","Diff"))
dev.off()



png("MDSplot_CellLine.png")
levels(SampleInfo$CellLine)
col.cellline <- c("blue","red","dark green")[SampleInfo$CellLine]
plotMDS(DGE_object,col=col.cellline)
legend("topleft",fill=c("blue","red","dark green"),legend=levels(SampleInfo$CellLine),cex=0.8)
title("Cell_Line")
dev.off()



png("RAPA_mds.png")
col.Rapa<- c("purple","orange")[SampleInfo$Rapa]
plotMDS(DGE_object,col=col.Rapa,pch=col.Diff)
legend("topright",c("yes","no"),fill=c("purple","orange"),legend=c("Yes","No"))
legend("bottomleft",c("FALSE","TRUE"),pch=pchD,legend=c("Non-Diff","Diff"))
title("Rapa treatment")
dev.off()



png("Differentiation_mds.png")
col.Diff<- c("red","blue")[SampleInfo$Diff_bin]
plotMDS(DGE_object,col=col.Diff)
legend("topleft",c("FALSE","TRUE"),fill=c("red","blue"), legend=c("Non-Diff","Diff"))
title("Diff vs Non-diff")
dev.off()



#Generate matrixplot

png("matrixplot.png")
my_cols <- c("green", "red", "blue")  #NTC,2.1,2.2
my_cols_Levels <- c("orange", "brown", "violet")  #244,62,ReN
my_cols_Diff <- c("grey","black")  #Nondiff,Diff
my_cols_Rapa<- c("purple","orange")  #FALSE,TRUE
my_cols_KO<-c("pink","black")#FALSE,TRUE

TA_pca <- prcomp(t(Countdata))
TA_pca.proportionvariances <- ((TA_pca$sdev^2) / (sum(TA_pca$sdev^2)))*100



png("PC1-PC5_cellline.png")
par(cex=1.0, cex.axis=0.8, cex.main=0.8)
pairs(TA_pca$x[,1:5], col=my_cols_Levels[SampleInfo$CellLine], main="Principal components_Cell Line", pch = c(21),  cex = 1, bg = c("orange", "brown", "violet"))
dev.off()

png("PC1-PC5_Diff.png")
par(cex=1.0, cex.axis=0.8, cex.main=0.8)
pairs(TA_pca$x[,1:5], col=my_cols_Diff[SampleInfo$Diff_bin], main="Principal components_Diff", pch = c(21),  cex = 1, bg = c("grey","black"))
dev.off()

png("PC1-PC5_Rapa.png")
par(cex=1.0, cex.axis=0.8, cex.main=0.8)
pairs(TA_pca$x[,1:5], col=my_cols_Rapa[SampleInfo$Rapa], main="Principal components_Rapa", pch = c(21),  cex = 1, bg = c("purple","orange"))
dev.off()

png("PC1-PC5_KO.png")
par(cex=1.0, cex.axis=0.8, cex.main=0.8)
pairs(TA_pca$x[,1:5], col=my_cols_KO[SampleInfo$KO_bin], main="Principal components_KO", pch = c(21),  cex = 1, bg = c("pink","black"))
dev.off()



###############################################
#						#
# 		02) WGCNA			#
#						#
###############################################



colData=data.frame(row.names = colnames(Countdata))
coldata <- SampleInfo[,c("Type","Rapa")]
coldata$Type<- factor(coldata$Type)
coldata$Rapa <- factor(coldata$Rapa)
all(rownames(coldata) %in% colnames(Countdata))# TRUE
all(rownames(coldata) == colnames(Countdata))#FALSE
Countdata <- Countdata[, rownames(coldata)]
all(rownames(coldata) == colnames(Countdata))


Dataset_pre= DESeqDataSetFromMatrix(countData = Countdata, colData=coldata,design=~Type+Rapa)


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
MEs = MEList$eigengenes # We selected these MEs
png("EigenGeneNetwork")
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))
dev.off()

save.image(file = "01_MappingAndAlignment_trim.RData")

