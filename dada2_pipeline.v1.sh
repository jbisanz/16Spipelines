#!/bin/bash
#$ -S /bin/bash
#$ -o dada2.outputlog.txt
#$ -e dada2.errorlog.txt
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=6G
#$ -l arch=linux-x64
#$ -l netapp=20G,scratch=1G
#$ -l h_rt=96:0:0
#$ -pe smp 24
date
hostname

#########################################################################
#JB Sep 24/2016
# Version 1
# Generates a table of ISUs/RSVs and corresponding taxonomies from multiple DBs
# Using QIIME 1.9.1, on UCSF QB3 cluster under OGS/GE 2011.11p1 with dada2_1.1.5 
# To do, remove reads mapping to PhiX which appear to show up as singleton ISUs (usearch or vsearch?)
# Will strip out these ISUs during the analysis with usearch

WORKING_DIR=/scrapp2/human_sep25_2016
FORWARD_READ=$WORKING_DIR/Reiner_Human_S0_L001_R1_001.fastq.gz
REVERSE_READ=$WORKING_DIR/Reiner_Human_S0_L001_R2_001.fastq.gz
INDEX_READ=$WORKING_DIR/Reiner_Human_S0_L001_I1_001.fastq.gz
MAPPING_FILE=$WORKING_DIR/mapping_file.txt

ERROR_LEARN=60 #learn error profile from this many randomly picked samples
SEED=182 #randomization seed
FTrim=10 # 5' trim on forward read
RTrim=10 # 5' trim on reverse read
ForTrunc=240 # forward read truncation length
RevTrunc=120 # reverse read truncation length


#########################################################################
source /netapp/home/jbisanz/qiime_for_cluster_1.9.1/settings/sourceme.sh
#########################################################################


#########################################################################
#demultiplexing using QIIME
if [ ! -d demultiplexed ]; then 

mkdir demultiplexed
mkdir demultiplexed/F
mkdir demultiplexed/R

echo "$(date)     Demultiplexing forward"
split_libraries_fastq.py -m $MAPPING_FILE --rev_comp_mapping_barcodes --barcode_type golay_12 -r 999 -n 999 -q 0 -p 0.0001 -i $FORWARD_READ -b $INDEX_READ --store_demultiplexed_fastq -o demultiplexed/
echo "$(date)     Splitting forward to individual files"
split_sequence_file_on_sample_ids.py --file_type fastq -i demultiplexed/seqs.fastq -o demultiplexed/F
	rm demultiplexed/histograms.txt
	rm demultiplexed/seqs.fna
	rm demultiplexed/seqs.fastq
	mv demultiplexed/split_library_log.txt demultiplexed/forward_split_library_log.txt
	
echo "$(date)     Demultiplexing reverse"
split_libraries_fastq.py -m $MAPPING_FILE --rev_comp_mapping_barcodes --barcode_type golay_12 -r 999 -n 999 -q 0 -p 0.0001 -i $REVERSE_READ -b $INDEX_READ --store_demultiplexed_fastq -o demultiplexed/
echo "$(date)     Splitting reverse to individual files"
split_sequence_file_on_sample_ids.py --file_type fastq -i demultiplexed/seqs.fastq -o demultiplexed/R
	rm demultiplexed/histograms.txt
	rm demultiplexed/seqs.fna
	rm demultiplexed/seqs.fastq
	mv demultiplexed/split_library_log.txt demultiplexed/reverse_split_library_log.txt


for f in demultiplexed/F/*
do
	mv $f $(echo $f | sed 's/\.fastq/_F.fastq/')
done

for r in demultiplexed/R/*
do
	mv $r $(echo $r | sed 's/\.fastq/_R.fastq/')
done

else
	echo "Expected that reads have already been demultiplexed based on directory structure"
fi
#########################################################################
echo "$(date)	Handing over to R"
export R_LIBS="/netapp/home/jbisanz/R/x86_64-pc-linux-gnu-library/3.2/"
export OMP_NUM_THREADS=$NSLOTS
module load MRO

#########################################################################
#########################################################################
#########################################################################
#########################################################################

#Rcode below using shell variables

/netopt/MRO/R/lib64/R/bin/R --no-save <<EOF
require(dada2)
forwards<-list.files("demultiplexed/F/")
reverses<-list.files("demultiplexed/R/")
sample.names <- gsub("_..+","",forwards) #this is ok because QIIME does not allow for _s in sample names
forwards<-file.path("demultiplexed/F/", forwards)
reverses<-file.path("demultiplexed/R/", reverses)


dir.create("results")
dir.create("results/read_quality")

pdf(paste("results/read_quality/","Forward_Quality.pdf",sep=""))
plotQualityProfile(forwards[1])
dev.off()

pdf(paste("results/read_quality/","Reverse_Quality.pdf",sep=""))
plotQualityProfile(reverses[1])
dev.off()

filtpath="filteredreads"
dir.create(filtpath)
filtFs <- file.path(filtpath, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filtpath, paste0(sample.names, "_R_filt.fastq.gz"))

for(i in seq_along(forwards)) {
  fastqPairedFilter(c(forwards[i], reverses[i]), c(filtFs[i], filtRs[i]), trimLeft=c($FTrim, $RTrim), truncLen=c($ForTrunc,$RevTrunc), maxN=0, maxEE=2, truncQ=2, compress=TRUE, verbose=TRUE)
}

##############################
#learn errors on subset 
set.seed($SEED) # Blink182/Scarface inspired seed
filts.learn <- sample(c(filtFs, filtRs), $ERROR_LEARN) # Pick subset samples (>100k reads/sample) to learn from
message(date(),      "Learning error rates from: ", filts.learn)
drp.learn <- derepFastq(filts.learn)
dd.learn <- dada(drp.learn, err=NULL, selfConsist=TRUE, multithread=TRUE)
err <- dd.learn[[1]]$err_out
rm(drp.learn);rm(dd.learn)
#############################

message(date(), "     Dereplicating")
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

message(date(), "     DADAing forward...")
dadaFs <- dada(derepFs, err=err, selfConsist = TRUE, multithread=TRUE)
message(date(), "     DADAing reverse...")
dadaRs <- dada(derepRs, err=err, selfConsist = TRUE, multithread=TRUE)

pdf("results/ErrorProfiles.pdf")
  plotErrors(dadaFs[[1]], nominalQ=TRUE)
  plotErrors(dadaRs[[1]], nominalQ=TRUE)
dev.off()

message(date(), "     Merging...")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

message(date(), "     Making table..")
seqtab <- makeSequenceTable(mergers)
message("Dimensions of seqtab: ", dim(seqtab))
write.table(seqtab, file="results/ISUtable.withchimeras.txt", sep='\t', quote=F, col.names=NA)

message(date(), "     Removing chimeras...")
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)
message("Fraction surviving: ", sum(seqtab.nochim)/sum(seqtab))
write.table(seqtab.nochim, file="results/ISUtable.txt", sep='\t', quote=F, col.names=NA)



message(date(), "     Assigning Taxonomy to species level with SILVA...")
taxa <- assignTaxonomy(seqtab.nochim, "/netapp/home/jbisanz/dada2_training_sets/silva_nr_v123_train_set.fa.gz")
taxa <- addSpecies(taxa, "/netapp/home/jbisanz/dada2_training_sets/silva_species_assignment_v123.fa.gz", verbose=TRUE, allowMultiple=TRUE)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
write.table(taxa, file="results/ISUtaxonomy.SILVA.txt", sep='\t', quote=F, col.names=NA)
##
message(date(), "     Assigning Taxonomy to species level with RDP...")
taxa <- assignTaxonomy(seqtab.nochim, "/netapp/home/jbisanz/dada2_training_sets/rdp_train_set_14.fa.gz")
taxa <- addSpecies(taxa, "/netapp/home/jbisanz/dada2_training_sets/rdp_species_assignment_14.fa.gz", verbose=TRUE, allowMultiple=TRUE)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
write.table(taxa, file="results/ISUtaxonomy.RDP.txt", sep='\t', quote=F, col.names=NA)
##
message(date(), "     Assigning Taxonomy to genus level with Green Genes...")
taxa <- assignTaxonomy(seqtab.nochim, "/netapp/home/jbisanz/dada2_training_sets/gg_13_8_train_set_97.fa.gz")
write.table(taxa, file="results/ISUtaxonomy.GG.txt", sep='\t', quote=F, col.names=NA)

message(date(), "    Complete")

EOF

