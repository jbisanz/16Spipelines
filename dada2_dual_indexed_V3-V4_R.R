
PROJECTNAME=Sys.getenv("PROJECTNAME")
WORKDIR=Sys.getenv("WORKDIR")
READS=Sys.getenv("READS")
SAMPLESHEET<-Sys.getenv("SAMPLESHEET")
TRIM=Sys.getenv("TRIM")
FILT=Sys.getenv("FILT")
INTERIM=Sys.getenv("INTERIM")
CUTADAPT=Sys.getenv("CUTADAPT")
print(paste("Analysis started at", date()))


#------------------------------------------------

## 1.1 Library Import 

require(dada2) #version 1.5.2
require(MicrobeR) #version 0.1.1
require(R.utils) # version 2.5.0
require(ggplot2) # version 2.2.1
require(plyr) # version 1.8.4
require(dplyr) # version 0.5.0
require(data.table) # version 1.10.4
require(reshape2) # version 4.5.6
require(plotly) # version 4.5.6
require(ape) # version 4.1
require(stringr) 

dir.create(INTERIM)
setwd(WORKDIR)
print(sessionInfo())


#------------------------------------------------

## 1.2 Raw Read Depth

#get metadata from illumina sample sheet; need to skip the first X lines - check the sample sheet!
metadata<-read.table(SAMPLESHEET, header=T, sep=',', stringsAsFactors = F, skip = 20, comment.char="", colClasses = "character")

#get fastq file names into new dataframe
reads<-data.frame(FASTQs=list.files(READS, include.dirs = FALSE, recursive = TRUE), stringsAsFactors = F)

#get the last 15 chars from file name, which contains whether reads are fwd or rev 
reads$Read<-sapply(reads$FASTQs, function(x){str_sub(x, -15)})

#get 2 characters indicating R1 (fwd) or R2 (rev)
reads$Read<-str_sub(reads$Read, 1, 2)

#create sample id by removing extraneous parts of file name 
reads$Sample_ID<-gsub("-[0-9]+/..+fastq\\.gz", "", reads$FASTQs)

#reshape data frame to have R1, R2 columns and populate with the fastq filenames
reads<-dcast(reads, Sample_ID~Read, value.var="FASTQs")

#count the number of lines in the fastq file; divide by 4 to get number of seq records
reads$UnFilt.NRead<-sapply(reads$R1, function(x){countLines(paste0(READS,x))[1]/4})

#join the metadata dataframe and the newly generated reads dataframe
metadata<-metadata %>% left_join(reads, by=c("Sample_ID"))
metadata<-data.table(metadata)

print("Summary stats per plate:")
print(ddply(metadata, c("Sample_Plate"), summarize, 
      n.actual=sum(UnFilt.NRead>=1, na.rm=T),
      n.expected=length(UnFilt.NRead),
      mean=mean(UnFilt.NRead, na.rm=T),
      sd=sd(UnFilt.NRead, na.rm=T),
      sem=sd(UnFilt.NRead, na.rm=T)/sqrt(length(UnFilt.NRead)),
      median=median(UnFilt.NRead, na.rm=T),
      min=min(UnFilt.NRead, na.rm=T),
      max=max(UnFilt.NRead, na.rm=T)
      ))

print("The following samples were not returned by Illumina Demultiplexing, check their barcodes:")
Nice.Table(metadata[is.na(R1)])


#------------------------------------------------

## 1.3 Adapter Removal

#remove the samples that were not returned by illumina demultiplexing
metadata<-metadata[!is.na(R1)]

#run cutadapt for each sample and save output in trim directory
if(!dir.exists(TRIM)){
dir.create(TRIM)
for (i in 1:nrow(metadata)){  
x<-metadata[i,]

#write results to log, delineating by sample id
system(paste("echo \"==============================", x$Sample_ID, "==============================\" >> trimlog.txt"))

#discard reads that do not have the primer seq
system(paste0(CUTADAPT, " --discard-untrimmed --pair-filter=both --error-rate=0.1 --times=2",

              #trim the primer sequences used in the first round PCR; g for fwd reads, G for rev reads
              " -g ^GTGCCAGCMGCCGCGGTAA -G ^GGACTACHVGGGTWTCTAAT",
              " --output=",paste0(TRIM,"/",x$Sample_ID,".R1.fastq.gz"),
              " --paired-output=",paste0(TRIM,"/",x$Sample_ID,".R2.fastq.gz"),
              " ",WORKDIR,"/",READS, x$R1,
              " ",WORKDIR,"/",READS, x$R2,
              " >> trimlog.txt"
              ), intern = F)
}
}

#add another column to metadata df with number of reads remaining after trimming primer seqs
metadata$Trimmed.NRead<-sapply(paste0(TRIM,"/",metadata$Sample_ID,".R1.fastq.gz"), function(x){countLines(x)[1]/4})
print(metadata)


#------------------------------------------------

## 1.4 Prefilter samples with low read numbers

#set the read number threshold by which to filter
threshold=5000

#remove this chunk from if statement, if filtering by read numbers is desired
if(FALSE){
paste0("The following samples had less than ", threshold, " reads and have been discarded:")
print(metadata[Trimmed.NRead<threshold])
metadata<-metadata[Trimmed.NRead>=threshold]
}


#------------------------------------------------

## 1.5 Read Quality




#------------------------------------------------
