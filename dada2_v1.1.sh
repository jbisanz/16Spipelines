#!/bin/bash
#$ -S /bin/bash
#$ -o dada2.outputlog.txt
#$ -e dada2.errorlog.txt
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=4G
#$ -l arch=linux-x64
#$ -l netapp=10G,scratch=10G
#$ -l h_rt=96:0:0
#$ -pe smp 24

echo "##################################################################################"
echo "$(date)	Running DADA2 Pipeline v1.1 for demultiplexed sequences from Basespace on node $(hostname) using $NSLOTS CPUs"
#Enter your settings below which will be exported to the environment for import to R
#Pipeline has been updated to version 1.70 of dada2 installed from github

export PATH=/netapp/home/jbisanz/bin/:$PATH

export PL_PROJECTNAME="MHGEP_Microbiota"
export PL_USERNAME="J. Bisanz"

export PL_WORKDIR="/turnbaugh/qb3share/jbisanz/MHGEP/" #The location where accessable from the QB3 cluster
export PL_READS=$PL_WORKDIR/MGHEP-39797785/ #name of read folder from Basespace 
export PL_SAMPLESHEET=$PL_WORKDIR/MGHEP_samplesheet.csv #illumina sample sheet from sequencer, the result of MakeSampleSheet.R

export PL_FTRUNC=200 #Trim for forward read, this should be fairly optimized for V4 reads *This is after the primer is removed
export PL_RTRUNC=180 #Trim for reverse read, this should be fairly optimized for V4 reads  *This is after the primer is removed

export PL_AMPLICONLENGTH=253 #expected bp for V4
export PL_LENGTHVAR=0.03 #Fraction to allow length variation in (0.03 is 3%)

export PL_FOR_PRIMER="GTGCCAGCMGCCGCGGTAA" #515F, replace if different primer
export PL_REV_PRIMER="GGACTACHVGGGTWTCTAAT" #806R, replace if different primer

export CLEANONEXIT="NO" #Should intermediate folders of trimmed and filtered reads be deleted on completion? YES/NO

############################################################################################
#Do not edit below unless you want to customize pipeline
############################################################################################

echo "$(date)	Handing over to R"
export R_LIBS="/netapp/home/jbisanz/R/x86_64-pc-linux-gnu-library/3.3/"
export OMP_NUM_THREADS=$NSLOTS
module load MRO


###
if [ -e ${PL_PROJECTNAME}.Rmd ]; then
	echo "Removing previous Rmd File: ${PROJECTNAME}.Rmd"
	rm ${PL_PROJECTNAME}.Rmd
fi

cat >> ${PL_PROJECTNAME}.Rmd<<'EOF'
---
title: '`r paste0("Project:", Sys.getenv("PL_PROJECTNAME"),": Read processing and QC V1.1.nov9.2017")`'
author: '`r Sys.getenv("PL_USERNAME")`'
date: '`r format(Sys.time(), "%Y-%m-%d %H:%M")`'
output: 
  html_document:
    code_folding: hide
    theme: spacelab
    highlight: monochrome
    fig_width: 11
    fig_height: 8.5
---

***

#Sections:

## 1. [Data Import, Denoising and QC](#S1)

* 1.1 [Variable Import](#S1.1)
* 1.2 [Data Import](#S1.2)
* 1.3 [Adapter Removal](#S1.3)
* 1.4 [Prefilter on reads](#S1.4)
* 1.5 [Read Quality](#S1.5)
* 1.6 [Filtering](#S1.6)
* 1.7 [Error Profile Learning](#S1.7)
* 1.8 [Dereplication](#S1.8)
* 1.9 [Denoising](#S1.9)
* 1.10 [Sequence Overlap](#S1.10)
* 1.11 [Sequence Table and In Silico Size Selection](#S1.11)
* 1.12 [Chimera Removal](#S1.12)
* 1.13 [Read loss tracking](#S1.13)
* 1.14 [Taxonomy Assignment](#S1.14)
* 1.15 [Alpha Diversity Metrics](#S1.15)
* 1.16 [Phylogenetic Tree Building](#S1.16)
* 1.17 [Beta Diversity Metrics](#S1.17)

## 2. [Standards and Controls](#S2)

* 2.1 [Community Standards](#S2.1)
* 2.2 [Negative Controls and GF Mice](#S2.2)

***


# 1. Variable Import


```{r setup, include=TRUE, message=F, warning=F}
PROJECTNAME=Sys.getenv("PL_PROJECTNAME")
WORKDIR=Sys.getenv("PL_WORKDIR")
READS=Sys.getenv("PL_READS")
SAMPLESHEET<-Sys.getenv("PL_SAMPLESHEET")
FTRUNC=Sys.getenv("PL_FTRUNC")
RTRUNC=Sys.getenv("PL_RTRUNC")
AMPLICONLENGTH=Sys.getenv("PL_AMPLICONLENGTH")
LENGTHVAR=Sys.getenv("PL_LENGTHVAR")
FOR_PRIMER=Sys.getenv("PL_FOR_PRIMER")
REV_PRIMER=Sys.getenv("PL_REV_PRIMER")


knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE, tidy=TRUE, cache=FALSE, fig.width=10, fig.height=7.5)
print(paste("Analysis started at", date()))
```

### User Settings:

```{r}
MicrobeR::Nice.Table(
data.frame(
Variable=c("PROJECTNAME","WORKDIR","READS","SAMPLESHEET","FTRUNC","RTRUNC","AMPLICONLENGTH","LENGTHVAR","FOR_PRIMER","REV_PRIMER"),
Value=c(PROJECTNAME,WORKDIR,READS,SAMPLESHEET,FTRUNC,RTRUNC,AMPLICONLENGTH,LENGTHVAR,FOR_PRIMER,REV_PRIMER)
))
```
Please check above and ensure that the variables were set appropriately.

***

## 1.1 Library Import {#S1.1}

```{r LibraryImport, message=FALSE, warning=FALSE}

library(foreach)
library(doParallel)
registerDoParallel(makeCluster(as.numeric(Sys.getenv("OMP_NUM_THREADS")))) #for parallel functions


library(dada2) #version 1.7.0
library(R.utils) 
library(tidyverse)
library(readxl)
library(stringr)
library(plotly)

library(ShortRead)
library(vegan)
library(DECIPHER)
library(ape)
library(ggtree)
library(phyloseq)


#setwd(WORKDIR)
dir.create("RDS")
dir.create("Filtered_reads")
dir.create("Trimmed_reads")
dir.create("Outputs")
print(sessionInfo())
```

***

## 1.2 Data Import {#S1.2}

```{r MetadataAndDepth}
metadata<-read_csv(SAMPLESHEET,skip = 19)

reads<-tibble(FASTQs=list.files(READS, include.dirs = FALSE, recursive = TRUE)) %>%
        mutate(Read=str_sub(FASTQs, -15) %>% str_sub(., 1, 2)) %>%
        mutate(Sample_ID=gsub("-[0-9]+/..+fastq\\.gz", "", FASTQs)) %>%
        spread(key="Read", value="FASTQs")

metadata<-metadata %>% left_join(reads)

if(nrow(metadata %>% filter(is.na(R1)))>0){
  print("The following samples were not returned by Illumina demultiplexing, check their barcodes and that fastq file exists:")
  metadata %>% filter(is.na(R1)) %>% MicrobeR::Nice.Table()
  metadata<-metadata %>% filter(!is.na(R1))
} else {
  print("All samples were returned from sequencing.")
}

```

***

## 1.3 Adapter Removal {#S1.3}

Now removing adapters allowing only 2 mismatches and no indels. This is running in parallel via foreach/doparallel so consider adjusting depending on your disk speeds and number of cores requested.
```{r AdapterTrim}
#Primer_For<-"GTGCCAGCMGCCGCGGTAA" #Forward primary PCR primer (515F)
#Primer_Rev<-"GGACTACHVGGGTWTCTAAT" #Reverse primary PCR primer (806R)

Primer_For<-FOR_PRIMER
Primer_Rev<-REV_PRIMER
MisMatches<-2 #allow 2 primer mismatches
if(!file.exists("RDS/adapter.RDS")){

ADAPTERTIME<-system.time({ReadStats<-
  foreach(Sample=metadata$Sample_ID) %dopar% {
	
    curjob<-subset(metadata, Sample_ID==Sample)
    if (nrow(curjob)>1){stop("Multiple samples with the same name detected!!! Please ensure no duplicate sample names, repicates should be processed separately and pooled later!!!")}
    
    Forward<-ShortRead::readFastq(paste0(READS,"/",curjob$R1))#assuming that per-sample fastq.gz files are sufficiently small<1Gb
    Reverse<-ShortRead::readFastq(paste0(READS,"/",curjob$R2))
      
    Forward_Fragment<-XVector::subseq(Forward@sread, start=1, end=nchar(Primer_For)) #grab first N bases to search for primer
    Reverse_Fragment<-XVector::subseq(Reverse@sread, start=1, end=nchar(Primer_Rev))
    
    Forward_Match<-Biostrings::vmatchPattern(Primer_For, Forward_Fragment, max.mismatch = 2, with.indels = FALSE, fixed=FALSE) #fixed=FALS for ambiguity codes
    Reverse_Match<-Biostrings::vmatchPattern(Primer_Rev, Reverse_Fragment, max.mismatch = 2, with.indels = FALSE, fixed=FALSE)
    Forward_Keep<-S4Vectors::elementNROWS(Forward_Match)==1
    Reverse_Keep<-S4Vectors::elementNROWS(Reverse_Match)==1
    Merger_Keep<-mapply(function(x,y) x==TRUE & y==TRUE, Forward_Keep, Reverse_Keep)
    
    Forward<-Forward[Merger_Keep]
    Reverse<-Reverse[Merger_Keep]
    
    Forward<-IRanges::narrow(Forward, start=nchar(Primer_For)+1)#trim off primer sequences
    Reverse<-IRanges::narrow(Reverse, start=nchar(Primer_Rev)+1)
    
    ShortRead::writeFastq(Forward, paste0("Trimmed_reads/",Sample,"_R1.fastq.gz"), compress=T)
    ShortRead::writeFastq(Reverse, paste0("Trimmed_reads/",Sample,"_R2.fastq.gz"), compress=T)
    
    return(list(Sample_ID=Sample, PreTrimming=length(Merger_Keep), PostTrimming=sum(Merger_Keep)))
  }
})[3]

ReadStatsDF<-as.data.frame(do.call(rbind, lapply(ReadStats, unlist))) %>%
          mutate(PreTrimming=as.integer(as.character(PreTrimming))) %>%
          mutate(PostTrimming=as.integer(as.character(PostTrimming))) %>%
          mutate(TrimmingLoss_Percent=round(100*(PostTrimming/PreTrimming),1))
          
metadata<-metadata %>% left_join(ReadStatsDF)

saveRDS(metadata, "RDS/adapter.RDS")

print(paste("Adapter removal finished in", round(ADAPTERTIME), "seconds with an average of", round(mean(metadata$TrimmingLoss_Percent),1),"% reads remaining after trimming."))
} else {
print("Adapters already removed...Skipping")
metadata<-readRDS("RDS/adapter.RDS")

}
```

***

## 1.4 Prefilter on reads {#S1.4}

Removing all samples with less than 5000 reads before going forward as these will typically represent failure to amplify or negative controls. If you had a very bad sequencing run, you could disable this step or reduce the cut off.

```{r Prefilter}
if(!file.exists("RDS/prefilter.RDS")){
print(paste("The following", nrow(metadata %>% filter(PostTrimming<=5000)),"samples had less than 5000 reads and have been discarded before moving forward"))
metadata %>% filter(PostTrimming<=5000) %>% MicrobeR::Nice.Table()
metadata<-metadata %>% filter(PostTrimming>5000)
saveRDS(metadata, "RDS/prefilter.RDS")
} else {
print("Prefilter already done.")
metadata<-readRDS("RDS/prefilter.RDS")
}

```

***

## 1.5 Read Quality {#S1.5}

Plotting 12 random samples quality profiles. Check these to make sure your trimming parameters were appropriate for your samples.

```{r QualityPlot}
Sample<-sample(metadata$Sample_ID, 12, replace = F)
plotQualityProfile(paste0("Trimmed_reads/",Sample,"_R1.fastq.gz")) + ggtitle("Forward Qualities")
plotQualityProfile(paste0("Trimmed_reads/",Sample,"_R2.fastq.gz")) + ggtitle("Reverse Qualities")
```

***

## 1.6 Filtering {#S1.6}

Filtering sequences allowing no Ns, max 2 expected errors and truncated to user setting.

```{r FilterAndTrim}
if(!file.exists("RDS/filtertrim.RDS")){
FILTERTIME<-system.time({
filt.sum<-filterAndTrim(fwd=paste0("Trimmed_reads/",metadata$Sample_ID,"_R1.fastq.gz"),
				    filt=paste0("Filtered_reads/",metadata$Sample_ID,"_R1.fastq.gz"), 
                        rev=paste0("Trimmed_reads/",metadata$Sample_ID,"_R2.fastq.gz"),
				    filt.rev=paste0("Filtered_reads/",metadata$Sample_ID,"_R2.fastq.gz"), 
                        truncLen=c(as.numeric(FTRUNC),as.numeric(RTRUNC)),
                        maxN=0,
                        maxEE=c(2,2),
                        truncQ=2,
                        rm.phix=TRUE,
                        compress=TRUE,
                        multithread=FALSE
                        )
})[3]


filt.sum<-filt.sum %>% as.data.frame %>%
  rownames_to_column("Sample_ID") %>%
  select(Sample_ID, PreFiltering=reads.in, PostFiltering=reads.out)

metadata<-metadata %>% left_join(filt.sum %>% mutate(Sample_ID=gsub("_R[12].fastq.gz","", Sample_ID)))
metadata<-metadata %>% mutate(FilteringLoss_Percent=round(100*(PostFiltering/PreTrimming),1))
print(paste("Quality filtering finished in", round(FILTERTIME), "seconds with an average of", round(mean(metadata$FilteringLoss_Percent),1),"% original reads remaining after quality filtering."))
saveRDS(metadata, "RDS/filtertrim.RDS")
} else {
	print("Filter and Trim already done.")
	metadata<-readRDS("RDS/filtertrim.RDS")
}


```

***

## 1.7 Error Profile Learning {#S1.7}

Learning error rate from 2 million reads.

```{r ErrorLearning}

if(!file.exists("RDS/R1.errprofile.RDS") & !file.exists("RDS/R2.errprofile.RDS")){
  LEARNTIME<-system.time({
    tmp<-capture.output(R1.errprofile<-learnErrors(paste0("Filtered_reads/",metadata$Sample_ID,"_R1.fastq.gz"), multithread=TRUE, nreads = 2e6, randomize=T))
      saveRDS(R1.errprofile,paste0("RDS","/R1.errprofile.RDS"))
  
    tmp<-capture.output(R2.errprofile<-learnErrors(paste0("Filtered_reads/",metadata$Sample_ID,"_R2.fastq.gz"), multithread=TRUE, nreads = 2e6, randomize=T))
      saveRDS(R2.errprofile,paste0("RDS","/R2.errprofile.RDS"))
  })[3]  
  
  print(paste("Error profile learning completed in", round(LEARNTIME), "seconds."))
} else {
  R1.errprofile<-readRDS("RDS/R1.errprofile.RDS")
  R2.errprofile<-readRDS("RDS/R2.errprofile.RDS")
  print("Error profiles loaded from disk.")
}


plotErrors(R1.errprofile, nominalQ=TRUE) + ggtitle("Forward Error Profile") + theme_bw()
plotErrors(R2.errprofile, nominalQ=TRUE) + ggtitle("Reverse Error Profile") + theme_bw()
```

The above error profiles should be an approximation of the expected. There is not a firm way to interpret this. Check [github.com/benjjneb/dada2/issues/](https://github.com/benjjneb/dada2/issues/) for guidance if these look very bad.

***

## 1.8 Dereplication {#S1.8}

```{r Dereplication}
if(!file.exists("RDS/derepFs.RDS") & !file.exists("RDS/derepRs.RDS")){
DEREPTIME<-system.time({
  derepFs <- derepFastq(paste0("Filtered_reads/",metadata$Sample_ID,"_R1.fastq.gz"), verbose=FALSE)
    names(derepFs)<-gsub("_R1.fastq.gz","", names(derepFs))
      saveRDS(derepFs, paste0("RDS","/derepFs.RDS"))


  derepRs <- derepFastq(paste0("Filtered_reads/",metadata$Sample_ID,"_R2.fastq.gz"), verbose=FALSE)
    names(derepRs)<-gsub("_R2.fastq.gz","", names(derepRs))
      saveRDS(derepRs, paste0("RDS","/derepRs.RDS"))
})[3]
print(paste("Dereplication completed in", round(DEREPTIME), "seconds."))
} else {
	derepFs<-readRDS("RDS/derepFs.RDS")
	derepRs<-readRDS("RDS/derepRs.RDS")
	print("Dereplicated reads loaded from file.")
}
```

***

## 1.9 Denoising {#S1.9}

```{r Denoising}
if(!file.exists("RDS/dadaFs.RDS") & !file.exists("RDS/dadaRs.RDS")){

DADATIME<-system.time({
  tmp<-capture.output(dadaFs <- dada(derepFs, err=R1.errprofile, multithread=TRUE))
  saveRDS(dadaFs, paste0("RDS","/dadaFs.RDS"))
  
  tmp<-capture.output(dadaRs <- dada(derepRs, err=R2.errprofile, multithread=TRUE))
  saveRDS(dadaRs, paste0("RDS","/dadaRs.RDS"))
})[3]

print(paste("Denoising completed in", round(DADATIME), "seconds."))

} else {
	dadaFs<-readRDS("RDS/dadaFs.RDS")
	dadaRs<-readRDS("RDS/dadaRs.RDS")
	print("Denoised reads loaded from file.")
	
}
```

***

## 1.10 Sequence Overlap {#S1.10}

```{r Overlap}
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=F)
saveRDS(mergers, "RDS/mergers.RDS")
print(paste("A total of", sum(mergers$accept),"reads successfully overlapped"))
```

***

## 1.11 Sequence Table and <i>In Silico</i> Size Selection {#S1.11}

```{r TableBuild}
ExpectedSize=as.numeric(AMPLICONLENGTH)
AllowVar=as.numeric(LENGTHVAR)

SVtable.unfilt<-makeSequenceTable(mergers)
saveRDS(SVtable.unfilt,"RDS/SVtable.unfilt.RDS")

print(paste("Before any filtering, there are", ncol(SVtable.unfilt),"SVs in", nrow(SVtable.unfilt), "samples with", sum(SVtable.unfilt), "reads."))

metadata<-metadata %>% left_join(data.frame(PostOverlap=rowSums(SVtable.unfilt)) %>% rownames_to_column("Sample_ID"))

SizeDist<-t(SVtable.unfilt) %>% as.data.frame %>%
  rownames_to_column("SVseq") %>%
  mutate(SVlength=nchar(SVseq)) %>%
  select(-SVseq) %>%
  group_by(SVlength) %>%
  summarize_all(sum) %>%
  mutate(AllSamples=rowSums(.[2:ncol(.)])) %>%
  select(SVlength, AllSamples)

saveRDS(SizeDist,"RDS/SizeDist.RDS")

ggplotly( 
  ggplot(SizeDist, aes(x=SVlength, y=AllSamples)) +
    geom_point() +
    theme_bw() +
    scale_y_log10() +
    ylab("Count across all samples") +
    xlab("SV length (bp)") +
    ggtitle("SV count by length")
)

SVtable.sizetrim<- SVtable.unfilt[,nchar(colnames(SVtable.unfilt)) %in% seq(round(ExpectedSize*(1-AllowVar)),round(ExpectedSize*(1+AllowVar)))]
metadata<-metadata %>% left_join(data.frame(PostSizeSelect=rowSums(SVtable.sizetrim)) %>% rownames_to_column("Sample_ID"))

print(paste("After in silico selection there are", ncol(SVtable.sizetrim),"SVs in", nrow(SVtable.sizetrim), "samples with", sum(SVtable.sizetrim), "reads."))
saveRDS(SVtable.sizetrim, "RDS/SVtable.sizetrim.RDS")
```


***

## 1.12 Chimera Removal {#S1.12}

Removing using pooled method.
```{r Chimeras}
SVtable<- removeBimeraDenovo(SVtable.sizetrim, method="pooled", multithread=TRUE, verbose=TRUE)
print(paste("Chimera removal removed", ncol(SVtable.sizetrim)-ncol(SVtable), "chimeras representing", round(100-100*(sum(SVtable)/sum(SVtable.sizetrim)),2), "% of all reads."))
print(paste("After chimera removal, there are", ncol(SVtable),"SVs in", nrow(SVtable), "samples in", sum(SVtable), "reads."))

metadata<-metadata %>% left_join(data.frame(FinalReads=rowSums(SVtable)) %>% rownames_to_column("Sample_ID"))
saveRDS(SVtable, "RDS/SVtable.final.RDS")
write_tsv(t(SVtable) %>% as.data.frame %>% rownames_to_column("Sample"), paste0("Outputs/",PROJECTNAME,".SVtable.tsv"))
```

Now will examine chimera rates on a per sample basis.

```{r}
metadata<-metadata %>% mutate(Percent_Chimera=round(100-100*(FinalReads/PostSizeSelect),2))

ggplotly(
  ggplot(metadata, aes(x=Percent_Chimera)) +
  geom_freqpoly(binwidth=0.5) +
  ylab("Number of Samples") +
  xlab("Percent Chimeric Reads") +
  ggtitle("Chimeric Reads per Sample") +
  theme_bw()
)
```

***

## 1.13 Read loss tracking {#S1.13}

The following plot shows the number of reads lost at any given step. Check to make sure there are no systematic heavey losses at any of the steps. If there are, parameters can be tuned to avoid this or problems with PCR and/or sequencing can be identified.

```{r ReadLossTracking}
metadata.final<-
  metadata %>%
    select(Sample_ID,
           Sample_Name,
           Sample_Plate,
           Sample_Well,
           I7_Index_ID,
           I7_Index_Seq=index,
           I5_Index_ID,
           I5_Index_Seq=index2,
           Sample_Project,
           Raw_reads=PreTrimming,
           PrimerTrimming=PostTrimming,
           QualityFilter=PostFiltering,
           Overlapping=PostOverlap,
           SizeSelection=PostSizeSelect,
           ChimeraRemoval=FinalReads,
           Percent_Chimera
           )

metadata.final<-metadata.final %>% mutate(TabledReads=ChimeraRemoval)

ggplotly(
  metadata.final %>% select(Sample_ID, 10:15) %>% mutate_if(is.numeric, function(x) 100*(x/metadata.final$Raw_reads)) %>%
    gather(-Sample_ID, key="Step", value="Reads") %>%
    mutate(Step=factor(Step, levels=c(colnames(metadata.final)[10:15]))) %>%
    ggplot(aes(x=Step, y=Reads, group=Sample_ID)) +
    geom_line(alpha=0.2) +
    theme_bw() +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(angle=45)) +
    ggtitle("Read loss by Sample")
)


```

***

## 1.14 Taxonomy Assignment {#S1.14}

Asigning taxonomy to the species level using SILVA 128 and Green Genes 13_8. These are being saved both as text and RDS files in Outputs and RDS respectively. RDS is better to load for work in R via `readRDS()` while TSV will be compatible with other approaches such as QIIME.

```{r TaxAssignment}
taxlevels<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

if(!file.exists("RDS/taxonomy_silva128.RDS") & !file.exists("RDS/taxonomy_gg_13_8")){
taxonomy <- assignTaxonomy(SVtable, "/turnbaugh/qb3share/shared_resources/databases/dada2_training_sets/silva_nr_v128_train_set.fa.gz", multithread=TRUE)
taxonomy <- addSpecies(taxonomy, "/turnbaugh/qb3share/shared_resources/databases/dada2_training_sets/silva_species_assignment_v128.fa.gz", allowMultiple=TRUE)
ggtaxonomy <- assignTaxonomy(SVtable, "/turnbaugh/qb3share/shared_resources/databases/dada2_training_sets/gg_13_8_train_set_97.fa.gz", multithread=TRUE)

colnames(taxonomy)<-taxlevels
colnames(ggtaxonomy)<-taxlevels

saveRDS(taxonomy, paste0("RDS","/taxonomy_silva128.RDS"))
saveRDS(ggtaxonomy, paste0("RDS","/taxonomy_gg_13_8.RDS"))
} else {
	print("Taxonomic assignment loaded from file.")
	taxonomy<-readRDS("RDS/taxonomy_silva128.RDS")
	ggtaxonomy<-readRDS("RDS/taxonomy_gg_13_8.RDS")
}
write_tsv(taxonomy %>% as.data.frame %>% rownames_to_column("SV"), paste0("Outputs/",PROJECTNAME,".taxonomy_silva128.tsv"))
write_tsv(ggtaxonomy %>% as.data.frame %>% rownames_to_column("SV"),paste0("Outputs/",PROJECTNAME,".taxonomy_gg_13_8.tsv"))


ggplotly(
data.frame(Assignment=apply(taxonomy, 2, function(x) {100*sum(!is.na(x))/length(x)})) %>%
  rownames_to_column("Level") %>%
  mutate(Level=factor(Level, levels=taxlevels)) %>%
  ggplot(aes(x=Level, y=Assignment)) + geom_point() +
  theme_bw() +
  xlab("Taxonomic Level") + 
  ylab("% SVs assigned") +
  ggtitle("SV Taxonomic Assignment by Level")
)
```

### Phylum Summary Plot

Use this to get a big picture look at your data. This plot will be useless unless you have very big differences in your data.

```{r,phylaplot}
physum<-t(SVtable) %>%
  as.data.frame %>%
  rownames_to_column("SV") %>%
  left_join(taxonomy %>% as.data.frame %>% rownames_to_column("SV") %>% select(SV, Phylum)) %>%
  select(-SV) %>%
  group_by(Phylum) %>%
  summarize_all(sum) %>%
  mutate_if(is.numeric, function(x) 100*(x/sum(x))) %>%
  gather(-Phylum, key="Sample", value="Abundance")

phyorder<-physum %>% select(-Sample) %>% group_by(Phylum) %>% summarize_all(mean) %>% arrange(Abundance)


sampleorder<-physum %>% spread(key="Sample", value="Abundance") %>% as.data.frame %>% mutate(Phylum=as.character(Phylum))
sampleorder$Phylum[is.na(sampleorder$Phylum)] <-"Unassigned"
sampleorder<-sampleorder %>% column_to_rownames("Phylum")

sampleorder<-hclust(dist(t(sampleorder))) 

  ggplotly(
  physum %>% mutate(Phylum=factor(Phylum, levels=phyorder$Phylum)) %>%
  mutate(Sample=factor(Sample, levels=sampleorder$labels[sampleorder$order])) %>%
  ggplot(aes(x=Sample, y=Abundance, fill=Phylum)) +
  geom_bar(stat="identity") +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  ylab("% Phylum Abundance")
)
```

***

## 1.15 Alpha Diversity Metrics {#S1.15}

Generating subsampled richness and non-subsampled Shannon diversity.
```{r,alphadiv}
richness<-data.frame(SV_Richness=specnumber(MicrobeR::Subsample.Table(t(SVtable)), MARGIN=2)) %>% rownames_to_column("Sample_ID")
shannon<-data.frame(Shannon_Diversity=diversity(SVtable, index="shannon"))%>% rownames_to_column("Sample_ID")

metadata.final<-metadata.final %>% left_join(richness) %>% left_join(shannon)

metadata.final %>% select(Sample_ID, Percent_Chimera, SV_Richness, Shannon_Diversity) %>%
  gather(-Sample_ID, -Percent_Chimera, key="Metric",value="Diversity") %>%
  ggplot(aes(x=Percent_Chimera, y=Diversity)) +
  geom_point(alpha=0.5, shape=20) +
  theme_bw() +
  xlab("% Chimeric Reads") +
  ylab("Alpha Diversity") +
  ggtitle("Relationship between Chimeric Reads and Alpha Diversity") +
  facet_wrap(~Metric, scales="free")


metadata.final %>% select(Sample_ID, Percent_Chimera, SV_Richness, Shannon_Diversity) %>%
  gather(-Sample_ID, -Percent_Chimera, key="Metric",value="Diversity") %>%
  ggplot(aes(x=Diversity)) +
  geom_freqpoly(bins=30) +
  theme_bw() +
  xlab("Alpha Diversity") +
  ylab("Number of Samples") +
  ggtitle("Alpha Diversity Distribution") +
  facet_wrap(~Metric, scales="free")

saveRDS(metadata.final, paste0("RDS/",PROJECTNAME,".metadata.tsv"))
write_tsv(metadata.final, paste0("Outputs/",PROJECTNAME,".metadata.tsv"))#write out final metadata that has alpha diversity metrics
```

## 1.16 Phylogenetic Tree Building {#S1.16}

Building as per [https://f1000research.com/articles/5-1492/v2](https://f1000research.com/articles/5-1492/v2) but using fasttreeMP with gtr instead. The parallel version will use the number of cores exported as OMP_NUM_THREADS. Tree is midpoint rooted.

```{r, phylotree}

if(!file.exists("RDS/tree.RDS")){
dir.create("alignments")
TREETIME<-system.time({
  lookup<-data.frame(Seq=colnames(SVtable), Name=paste0("SV_", 1:ncol(SVtable)))
  seqs<-lookup$Seq
  names(seqs)<-lookup$Name
  alignment <- DECIPHER::AlignSeqs(DNAStringSet(seqs), anchor=NA, verbose=F)
  names(alignment)<-lookup$Name #names lost for some reason
  writeXStringSet(alignment, paste0("Outputs/",PROJECTNAME,".SValigned.fasta")) #for some reason the DNAstring to DNAbin conversion was failing
 log<-system(paste("FastTreeMP -nt -gtr <",paste0("Outputs/",PROJECTNAME,".SValigned.fasta"), ">", paste0("Outputs/",PROJECTNAME,".tree")), intern=T)
})[3]

tree<-read.tree(paste0("Outputs/",PROJECTNAME,".tree"))
tree$tip.label<-lookup[match(tree$tip.label, lookup$Name),]$Seq #switch label back
tree<-phangorn::midpoint(tree)#midpoint root
write.tree(tree, paste0("Outputs/",PROJECTNAME,".tree"))
saveRDS(tree, "RDS/tree.RDS")

print(paste("Tree built in", round(TREETIME), "seconds."))
}

tree<-read.tree(paste0("Outputs/",PROJECTNAME,".tree"))

ggtree(tree) %<+%
  (as.data.frame(taxonomy[tree$tip.label, c("Phylum","Genus")]) %>% rownames_to_column("SV") ) +
  #geom_tiplab2(aes(color = Phylum, label=Genus),  size = 2) +
  geom_tippoint(aes(color = Phylum),  size = 2) +
  theme(legend.position = "right") +
  ggtitle("SV phylogenetic tree")

```

## 1.17 Beta Diversity Metrics {#S1.17}

```{r}
sub.SVtable<-vegan::rrarefy(SVtable, sample=min(rowSums(SVtable))) %>% apply(., 1, function(x) x/sum(x) )

UW.UniFrac<-phyloseq::UniFrac(phyloseq(otu_table(sub.SVtable, taxa_are_rows=T), phy_tree(tree)), weighted=F, parallel=T)
W.UniFrac<-phyloseq::UniFrac(phyloseq(otu_table(sub.SVtable, taxa_are_rows=T), phy_tree(tree)), weighted=T, parallel=T)
Bray<-vegan::vegdist(t(sub.SVtable), method="bray")

saveRDS(UW.UniFrac, paste0("RDS/",PROJECTNAME,".UnweightedUniFrac.RDS"))
saveRDS(W.UniFrac, paste0("RDS/",PROJECTNAME,".WeightedUniFrac.RDS"))
saveRDS(Bray, paste0("RDS/",PROJECTNAME,".BrayCurtis.RDS"))

write_tsv(UW.UniFrac %>% as.matrix %>% as.data.frame %>% rownames_to_column("Sample"), paste0("Outputs/",PROJECTNAME,".UnweightedUniFrac.tsv"))
write_tsv(W.UniFrac %>% as.matrix %>% as.data.frame %>% rownames_to_column("Sample"), paste0("Outputs/",PROJECTNAME,".WeightedUniFrac.tsv"))
write_tsv(Bray %>% as.matrix %>% as.data.frame %>% rownames_to_column("Sample"), paste0("Outputs/",PROJECTNAME,".BrayCurtis.tsv"))


distvectors<-rbind(
  ape::pcoa(UW.UniFrac)$vectors %>%
    as.data.frame %>% rownames_to_column("Sample_ID") %>%
    mutate(Metric="UW.Unifrac") %>% select(Sample_ID, Axis.1, Axis.2, Axis.3, Metric),
  ape::pcoa(W.UniFrac)$vectors %>%
    as.data.frame %>% rownames_to_column("Sample_ID") %>%
    mutate(Metric="W.Unifrac") %>% select(Sample_ID, Axis.1, Axis.2, Axis.3, Metric),
  ape::pcoa(Bray)$vectors %>% 
    as.data.frame %>% rownames_to_column("Sample_ID") %>%
    mutate(Metric="BrayCurtis") %>% select(Sample_ID, Axis.1, Axis.2, Axis.3, Metric)
)

ggplotly(
  ggplot(distvectors, aes(x=Axis.1, y=Axis.2, label=Sample_ID)) +
    geom_point() +
    facet_wrap(~Metric, scales="free") +
    theme_bw()
)

```

Distance matrices have been stored in RDS and tabular form.


# Standards and Controls {#S2}

## 2.1 Negative Controls and GF Mice {#S2.1}

Looking for negative controls and germ-free mice by the sample IDs containing the following: NTC, CON, GF. 

```{r dists}

control_samples<-metadata.final %>% filter(grepl("NTC", Sample_ID) | grepl("CON", Sample_ID) | grepl("GF", Sample_ID))

if(nrow(control_samples)==0){print("No control samples found")} else {
  tm<-metadata.final %>%
    mutate(Negative_Control=ifelse(Sample_ID %in% control_samples$Sample_ID, "Negative Control", "Sample"))
  ggplot(tm, aes(x=TabledReads)) +
    geom_freqpoly(binwidth=1000) +
    theme_classic() +
    xlab("Reads") +
    ylab("Number of samples") +
    geom_text(data=tm %>% filter(Negative_Control=="Negative Control"), aes(x=TabledReads, y=50, label=Sample_ID), color="darkred", position=position_jitter(width=0, height=20), hjust=0) +
    geom_vline(xintercept = (tm %>% filter(Negative_Control=="Negative Control"))$TabledReads, color="red", linetype="dashed")

ggplotly(
  metadata.final %>%
      mutate(Negative_Control=ifelse(Sample_ID %in% control_samples$Sample_ID, "Negative Control", "Sample")) %>%
      left_join(distvectors) %>%
      ggplot(aes(x=Axis.1, y=Axis.2, color=Negative_Control)) +
      geom_point(alpha=0.5) +
      theme_bw() +
      facet_wrap(~Metric, scales="free") +
      ggtitle("PCoA showing Negative Controls") +
      scale_color_manual(values = c("indianred","cornflowerblue"))
) 
  
  print("Closest matches to negative controls:")
  as.data.frame(as.matrix(Bray)) %>% rownames_to_column("Sample_ID") %>%
  filter(Sample_ID %in% control_samples$Sample_ID) %>%
  gather(-Sample_ID, key="Closest_Match", value="BrayCurtisDism") %>%
  group_by(Sample_ID) %>% arrange(BrayCurtisDism) %>%
  filter(Sample_ID != Closest_Match) %>%
  filter(BrayCurtisDism == min(BrayCurtisDism)) %>%
  MicrobeR::Nice.Table()
}

```


##2.2 Community Standards {#S2.2}

Looking for controls containing the following: EXTSTD (ZYMO Community Standard D6300), DNASTD (Zymo DNA Standard D6305/D6306), and HUMSTD (Previously sequenced human). Note: currently this is the expected copy number, not the lot controlled copy number. This will be adjusted in the future as there is definite variability between zymo batches. Note there is a particular sequence (`TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCGAAACTGGCAGGCTTGAGTCTTGTAGAGGGGGGGAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACAAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGG`) that appears to either be a sequencing error found across multiple batches or a true variant in the Salmonella genome that appears to be present at a ratio of 1:5. The actual genome from Zymo is highly fragmented with obvious assembly issues around the rRNA operon so this is possibly correct.

### Zymo Community Standard (Extraction)

```{r extstd}
control_samples<-metadata.final %>% filter(grepl("EXTSTD", Sample_ID))

if(nrow(control_samples)==0){print("No extraction standards found")} else{
  theoretical<-read_excel("/turnbaugh/qb3share/shared_resources/Zymo_Standards/ZymoLotInfo.xlsx", sheet="CopyNumber")
  theoretical<-theoretical %>% select(Name, `16S_Perc`, Sequence, Sequence2)  %>%
                gather(-Name, -`16S_Perc`, key="SequenceID", value="Seq") %>% filter(!is.na(Seq))
    
  #match independently as matchprobepair will not handle ambiguous bases
  theoretical$Forward_Start<-sapply(theoretical$Seq, function(x) matchPattern(FOR_PRIMER,DNAString(x), fixed=F)@ranges@start[1])
  theoretical$Forward_Stop<-sapply(theoretical$Seq, function(x) end(matchPattern(FOR_PRIMER,DNAString(x), fixed=F)@ranges[1]))
  theoretical$Reverse_Start<-sapply(theoretical$Seq, function(x) matchPattern(reverseComplement(DNAString(REV_PRIMER)),DNAString(x), fixed=F)@ranges@start[1])
  theoretical$Reverse_Stop<-sapply(theoretical$Seq, function(x) end(matchPattern(reverseComplement(DNAString(REV_PRIMER)),DNAString(x), fixed=F)@ranges[1]))
  
  theoretical$Amplicon<-mapply(function(fp,rp,seq){narrow(seq, start=fp+1, end=rp-1)},
                               theoretical$Forward_Stop, theoretical$Reverse_Start, theoretical$Seq)
  
  theoretical<-theoretical[!duplicated(theoretical$Amplicon),]
  
  controltab<-t(SVtable)[,control_samples$Sample_ID]
  controltab<-apply(controltab, 2, function(x) 100*(x/sum(x)))
  controltab<-controltab[rowSums(controltab)>0,] %>%
  as.data.frame %>%
    rownames_to_column("Amplicon") %>%
    left_join(theoretical %>% select(Name, `16S_Perc`, Amplicon) %>% dplyr::rename(Reference=`16S_Perc`)) %>%
    mutate(Name=replace(Name, is.na(Name), "Other")) %>%
    select(Name, Reference, get(control_samples$Sample_ID), Amplicon) %>%
    mutate(Name=factor(Name, levels=c("Lactobacillus_fermentum", 
                                      "Listeria_monocytogenes", 
                                      "Bacillus_subtilis", 
                                      "Staphylococcus_aureus", 
                                      "Salmonella_enterica", 
                                      "Enterococcus_faecalis", 
                                      "Escherichia_coli", 
                                      "Pseudomonas_aeruginosa",
                                      "Other"))) %>%
    arrange(Name)
  
controltab %>% left_join(taxonomy %>% as.data.frame %>% rownames_to_column("Amplicon")) %>% MicrobeR::Nice.Table()
controltab<-controltab %>% select(-Amplicon) %>% group_by(Name) %>% summarize_all(sum)

ggplotly(
 controltab%>%
  mutate(Name=factor(Name, levels=rev(levels(Name)))) %>%
  gather(-Name, key="Sample", value="Abundance") %>%
  ggplot(aes(x=Sample,y=Abundance,fill=Name)) +
  geom_bar(stat="identity") +
  theme_classic() +
   ggtitle("DNA Extraction Control Composition")
)


  print(paste("Average level of cross contamination in DNA extraction control samples is",
              controltab %>%
                select(-Reference) %>%
                gather(-Name, key="sample", value="abundance") %>%
                filter(Name=="Other") %>%
                summarize(mean(abundance)) %>% as.numeric(),
              "%"))
}



```

### Zymo DNA Standard (PCR)


```{r dnastd}
control_samples<-metadata.final %>% filter(grepl("DNASTD", Sample_ID))

if(nrow(control_samples)==0){print("No DNA standards found")} else{
  theoretical<-read_excel("/turnbaugh/qb3share/shared_resources/Zymo_Standards/ZymoLotInfo.xlsx", sheet="CopyNumber")
  theoretical<-theoretical %>% select(Name, `16S_Perc`, Sequence, Sequence2)  %>%
                gather(-Name, -`16S_Perc`, key="SequenceID", value="Seq") %>% filter(!is.na(Seq))
    
  #match independently as matchprobepair will not handle ambiguous bases
  theoretical$Forward_Start<-sapply(theoretical$Seq, function(x) matchPattern(FOR_PRIMER,DNAString(x), fixed=F)@ranges@start[1])
  theoretical$Forward_Stop<-sapply(theoretical$Seq, function(x) end(matchPattern(FOR_PRIMER,DNAString(x), fixed=F)@ranges[1]))
  theoretical$Reverse_Start<-sapply(theoretical$Seq, function(x) matchPattern(reverseComplement(DNAString(REV_PRIMER)),DNAString(x), fixed=F)@ranges@start[1])
  theoretical$Reverse_Stop<-sapply(theoretical$Seq, function(x) end(matchPattern(reverseComplement(DNAString(REV_PRIMER)),DNAString(x), fixed=F)@ranges[1]))
  
  theoretical$Amplicon<-mapply(function(fp,rp,seq){narrow(seq, start=fp+1, end=rp-1)},
                               theoretical$Forward_Stop, theoretical$Reverse_Start, theoretical$Seq)
  
  theoretical<-theoretical[!duplicated(theoretical$Amplicon),]
  
  controltab<-t(SVtable)[,control_samples$Sample_ID]
  controltab<-apply(controltab, 2, function(x) 100*(x/sum(x)))
  controltab<-controltab[rowSums(controltab)>0,] %>%
  as.data.frame %>%
    rownames_to_column("Amplicon") %>%
    left_join(theoretical %>% select(Name, `16S_Perc`, Amplicon) %>% dplyr::rename(Reference=`16S_Perc`)) %>%
    mutate(Name=replace(Name, is.na(Name), "Other")) %>%
    select(Name, Reference, get(control_samples$Sample_ID), Amplicon) %>%
    mutate(Name=factor(Name, levels=c("Lactobacillus_fermentum", 
                                      "Listeria_monocytogenes", 
                                      "Bacillus_subtilis", 
                                      "Staphylococcus_aureus", 
                                      "Salmonella_enterica", 
                                      "Enterococcus_faecalis", 
                                      "Escherichia_coli", 
                                      "Pseudomonas_aeruginosa",
                                      "Other"))) %>%
    arrange(Name)
  
controltab %>% left_join(taxonomy %>% as.data.frame %>% rownames_to_column("Amplicon")) %>% MicrobeR::Nice.Table()
controltab<-controltab %>% select(-Amplicon) %>% group_by(Name) %>% summarize_all(sum)

ggplotly(
 controltab%>%
  mutate(Name=factor(Name, levels=rev(levels(Name)))) %>%
  gather(-Name, key="Sample", value="Abundance") %>%
  ggplot(aes(x=Sample,y=Abundance,fill=Name)) +
  geom_bar(stat="identity") +
  theme_classic() +
   ggtitle("PCR Control Composition")
)

  print(paste("Average level of cross contamination in DNA extraction control samples is",
              controltab %>%
                select(-Reference) %>%
                gather(-Name, key="sample", value="abundance") %>%
                filter(Name=="Other") %>%
                summarize(mean(abundance)) %>% as.numeric(),
              "%"))
}

```

### Human DNA Standard (PCR)

This will only be valid for V4 amplicons and will be skipped if 515F and 806R not used.

```{r humstd}
control_samples<-metadata.final %>% filter(grepl("HUMSTD", Sample_ID))
#Primer_For<-"GTGCCAGCMGCCGCGGTAA" #Forward primary PCR primer (515F)
#Primer_Rev<-"GGACTACHVGGGTWTCTAAT" #Reverse primary PCR primer (806R)

if(nrow(control_samples)>0 & Primer_For=="GTGCCAGCMGCCGCGGTAA" & Primer_Rev=="GGACTACHVGGGTWTCTAAT"){
  hum3<-readRDS("/turnbaugh/qb3share/shared_resources/Zymo_Standards/HUMSTD.RDS")
  hum3$HUMSTD_reference<-100*(hum3$HUMSTD_reference/sum(hum3$HUMSTD_reference))
  hum3<-hum3 %>% filter(!is.na(HUMSTD_reference))
  
  merger<-t(SVtable)[,control_samples$Sample_ID] %>% data.frame %>% rownames_to_column("SV") %>%
    full_join(hum3)

  matr<-merger[,c("SV","HUMSTD_reference", control_samples$Sample_ID)] %>% data.frame %>% column_to_rownames("SV")
  matr[is.na(matr)]<-0
  matr<-matr[rowSums(matr)>0,]
  matr<-apply(matr, 2, function(x) 100*(x/sum(x)))
  matr<-matr %>% as.data.frame %>% rownames_to_column("SV")
  
  matr %>% gather(-SV, -HUMSTD_reference, key="Sample",value="Abundance") %>%
    ggplot(aes(x=HUMSTD_reference, y=Abundance)) +
    geom_point() +
    facet_wrap(~Sample) +
    geom_abline(color="red", linetype="dashed") +
    theme_bw() +
    xlab("% Abundance Reference") +
    ylab("% Abundance Sample")
}else{
	print("No DNA standards found or V4 primers not used")
}
```

```{r}
print(paste("Analysis finished at", date()))
```

***

Processing Complete!!! Find RDS outputs in RDS/ and TSV outputs in Outputs/. RDS files are ideal to read into R using `readRDS()`.

EOF

###
Rscript -e "rmarkdown::render('${PL_PROJECTNAME}.Rmd')"


if [ $CLEANONEXIT == "YES" ]; then
	rm -r Trimmed_reads
	rm -r Filtered_reads
fi
