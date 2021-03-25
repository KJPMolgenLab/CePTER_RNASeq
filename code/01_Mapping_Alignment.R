###############################################
#						#
# 	01) Mapping And Alignment		#
#						#
###############################################

mode = "normal"  # test only creates subset, "reanalyze" will calculate all BAM  files  again

fastqcfiles = "/scratch/fuchs/agmisc/chiocchetti/RNASeq_Data/Cepter/Rawfiles_resent/"
refhgfolder = "/scratch/fuchs/agmisc/chiocchetti/ReferenceGenomes/"
reffile = "hg38.fa.gz"
output="/scratch/fuchs/agmisc/chiocchetti/RNASeq_Data/Cepter/Output/"

dir.create(output)


require("parallel")
require("Rsubread")


################################
#			                   #
#	Preprocessing & QC	       #
#				               #
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


ncores = detectCores()


##Build Rerenze genome index

setwd(refhgfolder)
if(! file.exists("hg38.reads")){
    buildindex(basename=paste0(refhgfolder,"hg38"),reference=reffile)
}

setwd(output)
##Read in the fastq files
fastq.files<-list.files(path=fastqcfiles, pattern=".fq$",full.names=TRUE)

if (mode=="test")
    fastq.files=fastq.files[1:5]

                            
if (mode == "reanalyze")
    align(index=paste0(refhgfolder,"hg38"),readfile1=fastq.files, nthreads = ncores)

if (mode=="normal"){
    chckfiles=paste0(fastq.files,".subread.BAM.summary")
    filestomake=setdiff(chckfiles, list.files(fastqcfiles,full.names = TRUE))
    filestomake=gsub(".subread.BAM.summary","",filestomake)
    align(index=paste0(refhgfolder,"hg38"),readfile1=filestomake, nthreads = ncores)
}


##Reading in Bam files
bam.files <- list.files(path =fastqcfiles, pattern = ".BAM$", full.names = TRUE)

#Getting the feature counts
fc <- featureCounts(bam.files, annot.inbuilt="hg38",countMultiMappingReads=TRUE,
                    nthreads = ncores)

#Deleting extra information from column names
Countdata <- data.frame(fc$counts)
colnames(Countdata)<-sub("_.*", "", colnames(Countdata))

save(Countdata, file=paste0(output,"Countmatrix.RData"))
