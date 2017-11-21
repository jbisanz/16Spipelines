#Script to demultiplex inline barcodes on paired reads where structure is NNNN[8bp barcode]PrimerSeq
#If other length primers are used then positions at line 55-64 will need to be altered
#V1 J. Bisanz Nov 20
#Run time for 24million reads with 24 CPUs with 4GB availble per CPU via a cluster node(dual Xeon E5-2680v4) with 256GB total RAM is
#Limiting reagent appears to be reading/writing.

V4F="GTGCCAGCMGCCGCGGTAA"
V4R="GGACTACHVGGGTWTCTAAT"
ForwardReads<-"16Sreads/AEFA0_forward.fastq.gz"
ReverseReads<-"16Sreads/AEFA0_reverse.fastq.gz"
OutputFolder<-"demux/" #for sample_R1.fastq.gz and sample_R2.fastq.gz
metadata<-"metadata.txt" #must have a column called "SampleID" with name, "Lbarcode" and "Rbarcode" with nucleotide barcodes
NSTREAM=5e6 #work on NSTREAM reads at a time
NMISMATCH=3 #maximum number of primer mismatches
CORES=24 #number of CPUs to use
############################################
############################################
############################################
library(tidyverse)
library(ShortRead)
library(data.table)
library(doParallel) #using for parallel writting and compression
library(foreach)
dir.create(OutputFolder)

AmbiguousString <- function(x){ #create a version with ambiguous bases
  seq <- gsub("Y","\\[CT\\]",x)
  seq <- gsub("R","\\[AG\\]",x)
  seq <- gsub("W","\\[AT\\]",seq)
  seq <- gsub("S","\\[GC\\]",seq)
  seq <- gsub("K","\\[TG\\]",seq)
  seq <- gsub("M","\\[CA\\]",seq)
  seq <- gsub("D","\\[AGT\\]",seq)
  seq <- gsub("V","\\[ACG\\]",seq)
  seq <- gsub("H","\\[ACT\\]",seq)
  seq <- gsub("B","\\[CGT\\]",seq)
  seq <- gsub("N","\\[ATCG\\]",seq)
  return(seq)
}

V4F<-AmbiguousString(V4F)
V4R<-AmbiguousString(V4R)
#V4F<-DNAString(V4F)
#V4R<-DNAString(V4R)

metadata<-read_tsv(metadata) %>%
  mutate(key=paste0(toupper(Lbarcode),"_",toupper(Rbarcode)))

keys<-new.env(hash = T)
sink<-apply(metadata[,c("SampleID","key")], 1, function(x) assign(x['key'],x['SampleID'], keys)) #create a hash-type object where key is ForwardBC_ReverseBC and value is sample name.


message("Setting up cluster")
clust<-makeCluster(CORES) #set up a $CORES core cluster
registerDoParallel(clust)
sink<-clusterEvalQ(clust, library("ShortRead")) #import libraries into all clusters
sink<-clusterEvalQ(clust, library("data.table"))
clusterExport(clust, c("V4F","V4R", "NMISMATCH","OutputFolder"), envir = keys) #export variables

i=0
Npass=0
Nfail=0
Nnotfound=0
Ntotal=0
message("Starting demultiplex")
totaltime<-system.time({
  streamF<-FastqStreamer(ForwardReads, n=NSTREAM) 
  streamR<-FastqStreamer(ReverseReads, n=NSTREAM) 
  repeat{
    message(paste(date(),"Reading"))
    fqF<-yield(streamF)
    fqR<-yield(streamR)
    
    if(length(fqF)==0){break}
    message(paste(date(),"Parsing"))
    
    tbd<-data.table(
      Index=paste0(narrow(fqF@sread,5,12),"_", narrow(fqR@sread,5,12)),
      Freadid=as.character(fqF@id),
      Rreadid=as.character(fqR@id),
      
      Fprimer=as.character(narrow(fqF@sread,13, 31)),
      Fseq=as.character(narrow(fqF@sread,32)),
      Fqual=as.character(narrow(fqF@quality,32)@quality),
      
      Rprimer=as.character(narrow(fqR@sread,13, 32)),
      Rseq=as.character(narrow(fqR@sread,33)),
      Rqual=as.character(narrow(fqR@quality,33)@quality)
    )
    
    rm(fqF) #clear up RAM
    rm(fqR)
    
    message(paste(date(),"Assigning IDs"))
    
    tbd$SampleID<-as.character(sapply(tbd$Index, function(x) keys[[x]]))
    #tbd$SampleID<-unlist(parSapply(clust, split(tbd$Index, ceiling(seq_along(tbd$Index)/round(NSTREAM/CORES))), function(x) keys[[x]]  ))
    
    #This is very slow, better to use base agrepl
    #   tbd$ValidFp<-parSapply(clust ,tbd$Fprimer, function(x) isMatchingAt(V4F, DNAString(x), at=1, max.mismatch=NMISMATCH, fixed="subject", with.indels=F))
    #   tbd$ValidRp<-parSapply(clust, tbd$Rprimer, function(x) isMatchingAt(V4R, DNAString(x), at=1, max.mismatch=NMISMATCH, fixed="subject", with.indels=F))
    #fixed=F allows for IUPAC matching of ambiguous bases but only in the primer when fixed="subject", ie Ns in the read are mismatches
    #using parsapply with 12 nodes decreases processing time from by order of magnitude as compared to regular sapply.
    #This took approximately 45  minutesfor 24 million reads using 24 cores.
    
    #tbd$ValidFp<-agrepl(V4F, tbd$Fprimer, max.distance=list(cost=1,substitutions=NMISMATCH, deletions=0, insertions=0), fixed=F)
    #tbd$ValidRp<-agrepl(V4R, tbd$Rprimer, max.distance=list(cost=1,substitutions=NMISMATCH, deletions=0, insertions=0), fixed=F)
    
    message(paste(date(),"Checking Primers"))
    #Break into smaller jobs to better to use parallel threads
    tbd$ValidFp<-unlist(parSapply(clust, split(tbd$Fprimer, ceiling(seq_along(tbd$Fprimer)/round(NSTREAM/CORES))), function(x) agrepl(V4F, x, max.distance=list(cost=1,substitutions=NMISMATCH, deletions=0, insertions=0), fixed=F)))
    tbd$ValidRp<-unlist(parSapply(clust, split(tbd$Rprimer, ceiling(seq_along(tbd$Rprimer)/round(NSTREAM/CORES))), function(x) agrepl(V4R, x, max.distance=list(cost=1,substitutions=NMISMATCH, deletions=0, insertions=0), fixed=F)))
    
    message(paste(date(),"Writing FASTQs"))
    
    tmp<-nrow(tbd)
    Ntotal=Ntotal+tmp
    Nnotfound=sum(tbd$SampleID=="NULL")
    tbd<-subset(tbd, ValidFp==T & ValidRp==T & SampleID!="NULL")
    Npass=Npass+nrow(tbd)
    Nfail=Nfail+(tmp-nrow(tbd))
    
    
    tbd$NewFreadid<-paste(tbd$SampleID, tbd$Index, tbd$Freadid)
    tbd$NewRreadid<-paste(tbd$SampleID, tbd$Index, tbd$Rreadid)
    
    demux<-split(tbd, tbd$SampleID)
    
    #sink<-foreach(sample=names(demux)) %dopar%{ #parallel writting seems to slow down.
    for(sample in names(demux)){
      writeFastq(ShortReadQ(
        sread=DNAStringSet(demux[[sample]]$Fseq), 
        quality=BStringSet(demux[[sample]]$Fqual),
        id=BStringSet(demux[[sample]]$NewFreadid)
      ), 
      paste0(OutputFolder,sample,"_R1.fastq"), 
      compress=F, 
      mode="a")
      
      writeFastq(ShortReadQ(
        sread=DNAStringSet(demux[[sample]]$Rseq), 
        quality=BStringSet(demux[[sample]]$Rqual),
        id=BStringSet(demux[[sample]]$NewRreadid)
      ), 
      paste0(OutputFolder,sample,"_R2.fastq"), 
      compress=F, 
      mode="a")
      
    }
    
    i<-i+1
    message(paste(date(),"-------->Processed", i*NSTREAM, "reads."))
  }
  close(streamF)
  close(streamR)
  
  message("Compressing reads.")
  
  sink<-foreach(file=list.files(OutputFolder)) %dopar%{
    system(paste0("gzip ", OutputFolder,"/",file)) #assuming system has gzip
  }
})[3]


stopCluster(clust)
print(paste("Finished in ", lubridate::seconds_to_period(totaltime)))
print(paste("Total reads processed:", Ntotal))
print(paste("Total reads passing filter:", Npass))
print(paste("Total reads failing filter:", Nfail))
print(paste("Total indices failing to map to sample:", Nnotfound))
print(paste("Percentage reads passing:", round(100*(Npass/Ntotal),2)))

