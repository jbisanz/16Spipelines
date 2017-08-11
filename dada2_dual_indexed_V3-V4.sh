#!/bin/bash
#$ -S /bin/bash                     
#$ -o dada2.stdout.txt
#$ -e dada2.stderr.txt
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=4G
#$ -l arch=linux-x64
#$ -l netapp=10G,scratch=10G
#$ -l h_rt=96:0:0
#$ -pe smp 16

echo "##################################################################################"
echo "$(date)	Running DADA2 Pipeline v1.0 for demultiplexed sequences from Basespace on node $(hostname)"
echo "$(date)	Assuming overlapping V3-V4 reads with primer sequences still embedded"

#Enter your settings below which will be exported to the environment for import to R
#Pipeline has been updated to version 1.5.2 of dada2 installed from github

export PROJECTNAME="Test"
export WORKDIR="."                                  #working directory
export READS="reads/"                               #name of read folder from Basespace folder
export SAMPLESHEET="sample_sheet.csv"               #illumina sample sheet from sequencer, the result of MakeSampleSheet.R
export TRIM="TrimmedSeqs/"                          #name of folder for trimmed sequences from cutadapt
export FILT="FilteredSeqs/"                         #name of folder for filtered sequences from dada2
export INTERIM="RDS/"                               #name of folder for staging temp files
export CUTADAPT="/netapp/home/jbisanz/qiime_for_cluster_1.9.1/dependencies/Python-2.7.11-build/bin/cutadapt" #location of cutadapt install

#########################################################################
source /netapp/home/jbisanz/.bash_profile
#########################################################################

echo "$(date)	Handing over to R"
export R_LIBS="/netapp/home/jbisanz/R/x86_64-pc-linux-gnu-library/3.3/"
export OMP_NUM_THREADS=16
module load MRO

###
if [ -e ${PROJECTNAME}.Rmd ]; then
	echo "Removing previous Rmd File: ${PROJECTNAME}.Rmd"
	rm ${PROJECTNAME}.Rmd
fi

cat >> ${PROJECTNAME}.Rmd<<'EOF'
---
title: '`r paste("Project:", Sys.getenv("PROJECTNAME"))`'
subtitle: "DADA2 Pipeline for Read Processing and QC of V3-V4 amplicon sequence data"
author:
    - "J. Bisanz, original (https://github.com/jbisanz/16Spipelines)"
    - "K. Lam, updated (https://github.com/itskathylam/16Spipelines)"
date: '`r format(Sys.time(), "%Y-%m-%d %H:%M")`'
output: 
  html_document:
    code_folding: hide 
    theme: spacelab
    highlight: monochrome
    fig_width: 11
    fig_height: 8.5
---

#Sections:

## 1. [Data Import, Denoising and QC](#S1)

* 1.1 [Library Import](#S1.1)
* 1.2 [Raw Read Depth](#S1.2)
* 1.3 [Adapter Removal](#S1.3)
* 1.4 [Prefilter on reads](#S1.4)
* 1.5 [Read Quality](#S1.5)
* 1.6 [Filtering](#S1.6)
* 1.7 [Error Profile Learning](#S1.7)
* 1.8 [Dereplication](#S1.8)
* 1.9 [Denoising](#S1.9)
* 1.10 [Sequence Overlap](#S1.10)
* 1.11 [Sequence Table](#S1.11)
* 1.12 [Chimera Removal](#S1.12)
* 1.13 [<i>In Silico</i> Size Selection](#S1.13)
* 1.14 [Read loss tracking](#S1.14)
* 1.15 [Taxonomy Assignment](#S1.15)

## 2. [Standards and Controls](#S2)

* 2.1 [Community Standards](#S2.1)
* 2.2 [Negative Controls and GF Mice](#S2.2)

## 3. [Completion](#S3)


***


# 1. Data Import, Denoising and QC {#S1}

```{r setup, include=TRUE, message=F, warning=F}

#import environment variables from the shell
PROJECTNAME=Sys.getenv("PROJECTNAME")
WORKDIR=Sys.getenv("WORKDIR")
READS=Sys.getenv("READS")
SAMPLESHEET<-Sys.getenv("SAMPLESHEET")
TRIM=Sys.getenv("TRIM")
FILT=Sys.getenv("FILT")
INTERIM=Sys.getenv("INTERIM")
CUTADAPT=Sys.getenv("CUTADAPT")
knitr::opts_chunk$set(echo = TRUE, message=TRUE, warning=FALSE, tidy=TRUE, cache=FALSE, fig.width=10, fig.height=7.5)
print(paste("Analysis started at", date()))
```

***
## 1.1 Library Import {#S1.1}
```{r LibraryImport, message=FALSE, warning=FALSE}

#load the packages
require(dada2) #version 1.5.2
require(MicrobeR) #version 0.1.1
require(R.utils) # version 2.5.0
require(stringr) # version 1.2.0
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
setwd(WORKDIR) #working directory
print(sessionInfo())
```

***
## 1.2 Raw Read Depth {#S1.2}

Get sequencing read depth from raw data

```{r MetadataAndDepth}

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
Nice.Table(ddply(metadata, c("Sample_Plate"), summarize, 
      n.actual=sum(UnFilt.NRead>=1, na.rm=T),
      n.expected=length(UnFilt.NRead),
      mean=mean(UnFilt.NRead, na.rm=T),
      sd=sd(UnFilt.NRead, na.rm=T),
      sem=sd(UnFilt.NRead, na.rm=T)/sqrt(length(UnFilt.NRead)),
      median=median(UnFilt.NRead, na.rm=T),
      min=min(UnFilt.NRead, na.rm=T),
      max=max(UnFilt.NRead, na.rm=T)
      ))

ggplot(metadata[!is.na(UnFilt.NRead)], aes(UnFilt.NRead, fill=Sample_Plate)) + geom_density(stat="bin", alpha=0.3, binwidth=5000) + theme_bw()

print("The following samples were not returned by Illumina Demultiplexing, check their barcodes:")
Nice.Table(metadata[is.na(R1)])
```

***
## 1.3 Adapter Removal {#S1.3}

Now removing adapters and being sure to only keep samples that had reads returned from Illumina demultiplexing.

```{r AdapterTrim}

#remove the samples that were not returned by illumina demultiplexing
metadata<-metadata[!is.na(R1)]

#run cutadapt for each sample and save output in trim directory
if(!dir.exists(TRIM)){
dir.create(TRIM)
for (i in 1:nrow(metadata)){  
x<-metadata[i,]

#write results to log, delineating by sample id
system(paste("echo \"==============================", x$Sample_ID, "==============================\" >> trimlog.txt"))

#discard reads that do not have the primer seq; 'pair-fiter=any' will discard the pair if either do not have primer
system(paste0(CUTADAPT, " --discard-untrimmed --pair-filter=any --error-rate=0.1 --times=2",

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
```

***
## 1.4 Prefilter on reads {#S1.4}

Optional: removing all samples with less than threshold number of reads before going forward.

```{r Prefilter}

#set the read number threshold by which to filter
threshold=5000

#remove this chunk from if statement, if filtering by read numbers is desired
if(FALSE){
print(paste("The following samples had less than", threshold, "reads and have been discarded:"))
print(metadata[Trimmed.NRead < threshold])
metadata<-metadata[Trimmed.NRead >= threshold]
}
```

***
## 1.5 Read Quality {#S1.5}

Plot quality for first several samples

```{r QualityPlot}

#sample the first several files for plotting quality profile
subset_num = 6
if(nrow(metadata)<6){
    subset_num = nrow(metadata)
}
samps<-sample(metadata$Sample_ID, subset_num, replace = F)

#plot quality profile for fwd and rev reads
plotQualityProfile(paste0(TRIM, samps, ".R1.fastq.gz")) + ggtitle("Forward Qualities")
plotQualityProfile(paste0(TRIM, samps, ".R2.fastq.gz")) + ggtitle("Reverse Qualities")
```

***
## 1.6 Filtering {#S1.6}

Trim the forward and reverse reads. The length will depend on amplicon size, primer trimming, and read length!

```{r FilterAndTrim}

#set the fwd and rev lengths to trim to; 
fwd_len = 270
rev_len = 220
print(paste("Forward trim length set to", fwd_len, "and reverse trim length set to", rev_len))

#make pre- and post-filter filenames and save output in filt dir
metadata$R1.trim<-paste0(TRIM,"/",metadata$Sample_ID,".R1.fastq.gz")
metadata$R2.trim<-paste0(TRIM,"/",metadata$Sample_ID,".R2.fastq.gz")
metadata$R1.filt<-paste0(FILT,"/",metadata$Sample_ID,".filtered.R1.fastq.gz")
metadata$R2.filt<-paste0(FILT,"/",metadata$Sample_ID,".filtered.R2.fastq.gz")
if(!dir.exists(FILT)){
dir.create(FILT)

#call filter function to trim lengths, filter out poor quality bases, and remove phix
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
```

***
## 1.7 Error Profile Learning {#S1.7}

Learn error profiles for forward and reverse reads

```{r ErrorLearning}

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
plotErrors(R1.errprofile, nominalQ=TRUE) + ggtitle("Forward Error Profile") + theme_bw()
plotErrors(R2.errprofile, nominalQ=TRUE) + ggtitle("Reverse Error Profile") + theme_bw()
```

***
## 1.8 Dereplication {#S1.8}

Dereplicate forward and reverse reads

```{r Dereplication}

#if they do not exist, generate the dereplicated reads 
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
```

***
## 1.9 Denoising {#S1.9}

Denoise forward and reverse reads, using the respective error profiles

```{r Denoising}

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
```

***
## 1.10 Sequence Overlap {#S1.10}

Merge forward and reverse reads

```{r Overlap}

#if they do not exist, generate the merged fwd and rev reads from denoised reads
if(!file.exists(paste0(INTERIM,"/mergers.RDS"))){
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
saveRDS(mergers, paste0(INTERIM,"/mergers.RDS"))

#use the merged reads if they have been previously generated and saved
} else {
mergers<-readRDS(paste0(INTERIM,"/mergers.RDS"))
}
```

***
## 1.11 Sequence Table {#S1.11}

Make table of sequence variants from the merged reads

```{r TableBuild}

#make the sequence variant table from merged reads
SVs.w.chim<-makeSequenceTable(mergers)
print(paste("There are", ncol(SVs.w.chim),"SVs in", nrow(SVs.w.chim), "Samples"))
```

***
## 1.12 Chimera Removal {#S1.12}

Removing chimeras

```{r Chimeras}

#remove chimeras from the sequence variant table
SVtable<- removeBimeraDenovo(SVs.w.chim, method="consensus", multithread=TRUE, verbose=TRUE)
print(paste("After chimera removal, there are", ncol(SVtable),"SVs in", nrow(SVtable), "Samples with the following size distribution:"))
table(nchar(getSequences(SVtable)))
```

***
## 1.13 <i>In Silico</i> Size Selection {#S1.13}

In silico Size select, check on future runs if this is still necessary. Alowing for expected size ± 20%

```{r InSilicoSize}

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
```

***
## 1.14 Read loss tracking {#S1.14}

Track loss of reads at each processing step

```{r ReadLossTracking}

#get stats from each processing step
metadata$Merged.NRead<-rowSums(SVs.w.chim)[metadata$Sample_ID]
metadata$Final.NRead<-rowSums(SVtable)[metadata$Sample_ID]
readplot<-melt(metadata, id.vars=c("Sample_ID","Sample_Plate"), measure.vars=c("UnFilt.NRead", "Trimmed.NRead", "Merged.NRead","Final.NRead"), variable.name="Step", value.name="ReadNumber")

#plot read loss/retention at each step
ggplot(readplot, aes(x=Step, y=ReadNumber, group=Sample_ID, color=Sample_Plate)) + geom_violin(aes(group=Step), alpha=0.2) + geom_line(alpha=0.5) + theme_bw()
```

***
## 1.15 Taxonomy Assignment {#S1.15}

Assign taxonomy to sequence variants using SILVA and Greengenes databases

```{r TaxAssignment}
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
ggplot(frac, aes(x=Level, y=Assigned)) + geom_bar(stat="identity") + ggtitle("Fraction SVs assigned Taxonomy") + theme_bw() + ylab("% Assigned")
```

Now transpose SV table for future use

```{r Transpose}

#transpose SV table and save
SVtable<-t(SVtable)
saveRDS(SVtable,paste0(INTERIM,"/SVtable.RDS"))
```

***

# 2. Standards and Controls {#S2}

***
## 2.1 Community Standards {#S2.1}

Check standards at the genus level. I have assumed that samples which are the community standards have "STANDARD" somewhere in their name.

```{r CommunityStandards}
if(sum(colnames(SVtable) %like% "STANDARD")>0){
standardstab<-SVtable[,colnames(SVtable) %like% "STANDARD"]
Microbiome.Barplot(Summarize.Taxa(standardstab, taxonomy)[[6]])
Nice.Table(merge(standardstab[rowSums(standardstab)>2,], taxonomy, all.x=T, all.y=F, by="row.names"))
} else {print("No STANDARDs found, skipping this...")}
```

***
## 2.2 Negative Controls and GF Mice {#S2.2}

Check Controls, assuming names contain "EXTRACTION", "NTC" and/or "GF".

```{r NegativeControls}
if(any(colnames(SVtable) %like% "EXTRACTION", colnames(SVtable) %like% "NTC", colnames(SVtable) %like% "GF")){
constab<-as.data.frame(SVtable[,grep("EXTRACTION|GF|NTC", colnames(SVtable))])
#Microbiome.Barplot(Summarize.Taxa(constab, taxonomy)[[6]])
Nice.Table(merge(constab[rowSums(constab)>2,], taxonomy, all.x=T, all.y=F, by="row.names"))
} else {print("No EXTRACTIONs, NTCs, or GFs found, skipping this...")}
```

***
# 3. Completion {#S3}

Write SV table and taxonomy assignments to .txt files

```{r Completion}
print(paste("Analysis complete at", date()))
write.table(SVtable, "SVtable.txt", sep='\t', quote=F, col.names=NA)
write.table(taxonomy, "TaxonomySilva.txt", sep='\t', quote=F, col.names=NA)
write.table(ggtaxonomy, "TaxonomyGG.txt", sep='\t', quote=F, col.names=NA)
```

Processing complete, use SVtable.RDS and taxonomy.RDS for futher analysis!

EOF

###
Rscript -e "rmarkdown::render('${PROJECTNAME}.Rmd')"

