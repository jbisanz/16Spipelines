
#import environment variables from the shell
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

#load the packages
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

#set up environment
dir.create(INTERIM)
setwd(WORKDIR)
print(sessionInfo())


#------------------------------------------------

## 1.2 Raw Read Depth

#get metadata from illumina sample sheet; need to skip the first X lines - double check the sample sheet!
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
print(metadata[is.na(R1)])


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

              #trim the primer sequences used in the first-round PCR; g for fwd reads, G for rev reads
              " -g ^CCTACGGGNGGCWGCAG -G ^GACTACHVGGGTATCTAATCC",
              " --output=",paste0(TRIM,"/",x$Sample_ID,".R1.fastq.gz"),
              " --paired-output=",paste0(TRIM,"/",x$Sample_ID,".R2.fastq.gz"),
              " ",WORKDIR,"/",READS, x$R1,
              " ",WORKDIR,"/",READS, x$R2,
              " >> trimlog.txt"
              ), intern = F)
}
}

#add another column to metadata dataframe with number of reads remaining after trimming primer seqs
metadata$Trimmed.NRead<-sapply(paste0(TRIM,"/",metadata$Sample_ID,".R1.fastq.gz"), function(x){countLines(x)[1]/4})

print("Summary of read number post trimming")
print(metadata)


#------------------------------------------------

## 1.4 Prefilter samples with low read numbers

#set the read number threshold by which to filter
threshold=5000

#remove this chunk from if statement, if filtering by read numbers is desired
if(FALSE){
print(paste("The following samples had less than", threshold, "reads and have been discarded:"))
print(metadata[Trimmed.NRead < threshold])
metadata<-metadata[Trimmed.NRead >= threshold]
}


#------------------------------------------------

## 1.5 Read Quality

#sample the first several files for plotting quality profile
subset_num = 6
if(nrow(metadata)<6){
    subset_num = nrow(metadata)
}
samps<-sample(metadata$Sample_ID, subset_num, replace = F)

#plot quality profile for fwd and rev reads
#plotQualityProfile(paste0(TRIM, samps, ".R1.fastq.gz")) + ggtitle("Forward Qualities")
#plotQualityProfile(paste0(TRIM, samps, ".R2.fastq.gz")) + ggtitle("Reverse Qualities")


#------------------------------------------------

## 1.6 Filtering and Trimming

#set the fwd and rev lengths to trim to; this depends on amplicon size, primer trimming, and read length!
fwd_len = 270
rev_len = 220

#make pre- and post-filter filenames and save output in filt dir
metadata$R1.trim<-paste0(TRIM,"/",metadata$Sample_ID,".R1.fastq.gz")
metadata$R2.trim<-paste0(TRIM,"/",metadata$Sample_ID,".R2.fastq.gz")
metadata$R1.filt<-paste0(FILT,"/",metadata$Sample_ID,".filtered.R1.fastq.gz")
metadata$R2.filt<-paste0(FILT,"/",metadata$Sample_ID,".filtered.R2.fastq.gz")
if(!dir.exists(FILT)){
dir.create(FILT)

#call filter function to trim lengths, filer out poor quality bases, and remove phix
filt.sum<-filterAndTrim(metadata$R1.trim, metadata$R1.filt, 
                        metadata$R2.trim, metadata$R2.filt, 
                        truncLen=c(fwd_len, rev_len),
                        maxN=0,
                        maxEE=c(2,2),
                        truncQ=2,
                        rm.phix=TRUE,
                        compress=TRUE,
                        multithread=TRUE
                        )
}


#------------------------------------------------

## 1.7 Error Profile Learning

#if they do not exist, generate the error profiles (using the first 2.5M reads) 
if(!file.exists(paste0(INTERIM,"/R2.errprofile.RDS"))){
R1.errprofile<-learnErrors(metadata$R1.filt, multithread=TRUE, nreads = 2.5e6)
R2.errprofile<-learnErrors(metadata$R2.filt, multithread=TRUE, nreads = 2.5e6)

saveRDS(R1.errprofile,paste0(INTERIM,"/R1.errprofile.RDS"))
saveRDS(R2.errprofile,paste0(INTERIM,"/R2.errprofile.RDS"))

#use the error profiles if they have been previously generated and saved
}else{
  R1.errprofile<-readRDS(paste0(INTERIM,"/R1.errprofile.RDS"))
  R2.errprofile<-readRDS(paste0(INTERIM,"/R2.errprofile.RDS"))
}

#plot the fwd and rev error profiles
#plotErrors(R1.errprofile, nominalQ=TRUE) + ggtitle("Forward Error Profile") + theme_bw()
#plotErrors(R2.errprofile, nominalQ=TRUE) + ggtitle("Reverse Error Profile") + theme_bw()


#------------------------------------------------

## 1.8 Dereplication

#if they do not exist, generate the deplicated reads 
if(!file.exists(paste0(INTERIM,"/derepRs.RDS"))){
derepFs <- derepFastq(metadata$R1.filt, verbose=TRUE)
derepRs <- derepFastq(metadata$R2.filt, verbose=TRUE)
names(derepFs)<-gsub("\\...+","", names(derepFs))
names(derepRs)<-gsub("\\...+","", names(derepRs))

saveRDS(derepFs, paste0(INTERIM,"/derepFs.RDS"))
saveRDS(derepRs, paste0(INTERIM,"/derepRs.RDS"))

#use the dereplicated reads if they have been previously generated and saved
} else {
  derepFs<-readRDS(paste0(INTERIM,"/derepFs.RDS"))
  derepRs<-readRDS(paste0(INTERIM,"/derepRs.RDS"))
}


#------------------------------------------------

## 1.9 Denoising

#if they do not exist, generate the denoised reads using the error profiles
if(!file.exists(paste0(INTERIM,"/dadaRs.RDS"))){
dadaFs <- dada(derepFs, err=R1.errprofile, multithread=TRUE)
dadaRs <- dada(derepRs, err=R2.errprofile, multithread=TRUE)

saveRDS(dadaFs, paste0(INTERIM,"/dadaFs.RDS"))
saveRDS(dadaRs, paste0(INTERIM,"/dadaRs.RDS"))

#use the denoised reads if they have been previously generated and saved
} else {
  dadaFs<-readRDS(paste0(INTERIM,"/dadaFs.RDS"))
  dadaRs<-readRDS(paste0(INTERIM,"/dadaRs.RDS"))
}


#------------------------------------------------

## 1.10 Sequence Overlap 

#if they do not exist, generate the merged fwd and rev reads from denoised reads
if(!file.exists(paste0(INTERIM,"/mergers.RDS"))){
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
saveRDS(mergers, paste0(INTERIM,"/mergers.RDS"))

#use the denoised reads if they have been previously generated and saved
} else {
mergers<-readRDS(paste0(INTERIM,"/mergers.RDS"))
}


#------------------------------------------------

## 1.11 Sequence Table

#make the sequence variant table from merged reads
SVs.w.chim<-makeSequenceTable(mergers)
print(paste("There are", ncol(SVs.w.chim),"SVs in", nrow(SVs.w.chim), "Samples"))


#------------------------------------------------

## 1.12 Chimera Removal

#remove chimeras from the sequence variant table
SVtable<- removeBimeraDenovo(SVs.w.chim, method="consensus", multithread=TRUE, verbose=TRUE)
print(paste("After chimera removal, there are", ncol(SVtable),"SVs in", nrow(SVtable), "Samples with the following size distribution:"))
table(nchar(getSequences(SVtable)))

#------------------------------------------------

## 1.13 <i>In Silico</i> Size Selection

#set expected length, allowing for 20% variation; note that primer sequence have been trimmed, so expected length = amplicon size minus primer lengths (V3-V4: 470 - 40)
ExpectedSize = 430
AllowVar = 0.2
lower_lim = ExpectedSize*(1-AllowVar)
upper_lim = ExpectedSize*(1+AllowVar)

#do the in silico size selection using upper and lower limits for expected length
SizeSelect<-SVtable[,nchar(colnames(SVtable)) %in% seq(lower_lim, upper_lim)]
print(paste("Expected fragment size after trimming primers is", ExpectedSize, "but allow for size variation of", AllowVar*100, "%"))
print(paste("Removed", ncol(SVtable)-ncol(SizeSelect), "Sequence Variants of", ncol(SVtable), "that were smaller than", lower_lim, "or greater than", upper_lim))
SVtable<-SizeSelect
print(paste("Post size selection, the size distribution is now:"))
table(nchar(getSequences(SVtable)))


#------------------------------------------------

## 1.14 Read loss tracking

#get stats from each processing step
metadata$Merged.NRead<-rowSums(SVs.w.chim)[metadata$Sample_ID]
metadata$Final.NRead<-rowSums(SVtable)[metadata$Sample_ID]
readplot<-melt(metadata, id.vars=c("Sample_ID","Sample_Plate"), measure.vars=c("UnFilt.NRead", "Trimmed.NRead", "Merged.NRead","Final.NRead"), variable.name="Step", value.name="ReadNumber")

#plot read loss/retention at each step
#ggplot(readplot, aes(x=Step, y=ReadNumber, group=Sample_ID, color=Sample_Plate)) + geom_violin(aes(group=Step), alpha=0.2) + geom_line(alpha=0.5) + theme_bw()


#------------------------------------------------

## 1.15 Taxonomy Assignment

#if they do not exist, generate the taxonomic assignments
if(!file.exists(paste0(INTERIM,"/taxonomy.RDS"))){

#use SILVA database
taxonomy <- assignTaxonomy(SVtable, "/netapp/home/jbisanz/dada2_training_sets/silva_nr_v123_train_set.fa.gz", multithread=TRUE)
taxonomy <- addSpecies(taxonomy, "/netapp/home/jbisanz/dada2_training_sets/silva_species_assignment_v123.fa.gz", allowMultiple=TRUE)

#use Greengenes database
ggtaxonomy <- assignTaxonomy(SVtable, "/netapp/home/jbisanz/dada2_training_sets/gg_13_8_train_set_97.fa.gz", multithread=TRUE)

saveRDS(taxonomy, paste0(INTERIM,"/taxonomy.RDS"))
saveRDS(ggtaxonomy, paste0(INTERIM,"/greengenes_taxonomy.RDS"))

#use the taxonomic assignments if they have been previously generated and saved
} else {
  taxonomy<-readRDS(paste0(INTERIM,"/taxonomy.RDS"))
  ggtaxonomy<-readRDS(paste0(INTERIM,"/greengenes_taxonomy.RDS"))
}

#get stats on taxonmic assignment
taxonomy<-as.data.frame(taxonomy)
frac<-apply(taxonomy, 2, function(x) {100*sum(!is.na(x))/length(x)})
frac<-data.frame(Level=names(frac), Assigned=frac)
frac$Level<-factor(frac$Level, levels=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

#plot fraction of sequence variants assigned at given taxonomic level
#ggplot(frac, aes(x=Level, y=Assigned)) + geom_bar(stat="identity") + ggtitle("Fraction SVs assigned Taxonomy") + theme_bw() + ylab("% Assigned")

#transpose SV table for later use and save
SVtable<-t(SVtable)
saveRDS(SVtable,paste0(INTERIM,"/SVtable.RDS"))


#------------------------------------------------

## Completion

print(paste("Analysis complete at", date()))
write.table(SVtable, "SVtable.txt", sep='\t', quote=F, col.names=NA)
write.table(taxonomy, "SilvaTaxonomy.txt", sep='\t', quote=F, col.names=NA)
write.table(ggtaxonomy, "GGTaxonomy.txt", sep='\t', quote=F, col.names=NA)