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
#$ -pe smp 12

echo "##################################################################################"
echo "$(date)	Running DADA2 Pipeline v1.4 on node $(hostname)"


#########################################################################
#JB Oct31/2016
# Version 1.5
# Generates a table of sequence variants and corresponding taxonomies from multiple DBs
# Using QIIME 1.9.1, on UCSF QB3 cluster under OGS/GE 2011.11p1 with dada2_1.1.5 
# Using code derived from http://benjjneb.github.io/dada2/bigdata.html
#
# Changes
#	1.4: now removing phix during dereplication and saving intermediate files for later use if required. 
#	1.5: saving computationally intensive intermediate files and not regenerating if rerun
#


WORKING_DIR=/scrapp2/INSERTHERE
FORWARD_READ=$WORKING_DIR/INSERTHERE
REVERSE_READ=$WORKING_DIR/INSERTHERE
INDEX_READ=$WORKING_DIR/INSERTHERE
MAPPING_FILE=$WORKING_DIR/INSERTHERE

ERROR_LEARN=60 #learn error profile from this many randomly picked samples, aim for ~2million reads
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

echo "$(date)	Beginning demultiplexing"

if [ ! -d demultiplexed ]; then 

mkdir demultiplexed
mkdir demultiplexed/F
mkdir demultiplexed/R

echo "$(date)     Demultiplexing $FORWARD_READ"
split_libraries_fastq.py -m $MAPPING_FILE --rev_comp_mapping_barcodes --barcode_type golay_12 -r 999 -n 999 -q 0 -p 0.0001 -i $FORWARD_READ -b $INDEX_READ --store_demultiplexed_fastq -o demultiplexed/
echo "$(date)     Splitting forward to individual files"
split_sequence_file_on_sample_ids.py --file_type fastq -i demultiplexed/seqs.fastq -o demultiplexed/F
	rm demultiplexed/histograms.txt
	rm demultiplexed/seqs.fna
	rm demultiplexed/seqs.fastq
	mv demultiplexed/split_library_log.txt demultiplexed/forward_split_library_log.txt
	
echo "$(date)     Demultiplexing $REVERSE_READ"
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

echo "$(date)     Summary of demultiplexing:"
head -n 12  $WORKING_DIR/demultiplexed/forward_split_library_log.txt


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


dir.create("dada2_results")
dir.create("dada2_intermediates")

pdf(paste("dada2_results/","Forward_Quality.pdf",sep=""))
	plotQualityProfile(forwards[1])
dev.off()

pdf(paste("dada2_results/","Reverse_Quality.pdf",sep=""))
	plotQualityProfile(reverses[1])
dev.off()

filtpath="filteredreads"
filtFs <- file.path(filtpath, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filtpath, paste0(sample.names, "_R_filt.fastq.gz"))


if (!file.exists(filtpath)){
dir.create(filtpath)
	for(i in seq_along(forwards)) {
  		fastqPairedFilter(c(forwards[i], reverses[i]), c(filtFs[i], filtRs[i]), trimLeft=c($FTrim, $RTrim), truncLen=c($ForTrunc,$RevTrunc), rm.phix=c(TRUE,TRUE), maxN=0, maxEE=2, truncQ=2, compress=TRUE, verbose=TRUE)
	}
}else { message("------->", date(),"	Fastq files already filered. Reading in from ", filtpath)}

##############################
#learn errors on subset 
set.seed($SEED)

filts.learn <- sample(filtFs, $ERROR_LEARN)
message("------->", date(), "	Learning Forward error rates from: ", filts.learn)
	
	if (!file.exists("dada2_intermediates/dada2_forward_error.rds")){
		drp.learn <- derepFastq(filts.learn)
		dd.learn <- dada(drp.learn, err=NULL, selfConsist=TRUE, multithread=TRUE)
		err_F <- dd.learn[[1]]$err_out
		rm(drp.learn);rm(dd.learn)
		saveRDS(err_F, "dada2_intermediates/dada2_forward_error.rds")
	} else {
		message("------->", date(),      "	Forward error already determined, locating from dada2_intermediates/dada2_forward_error.rds")
		err_F<-readRDS("dada2_intermediates/dada2_forward_error.rds")
	}


filts.learn <- sample(filtRs, $ERROR_LEARN)
message("------->", date(),  "	Learning Reverse error rates from: ", filts.learn)

	if (!file.exists("dada2_intermediates/dada2_reverse_error.rds")){

		drp.learn <- derepFastq(filts.learn)
		dd.learn <- dada(drp.learn, err=NULL, selfConsist=TRUE, multithread=TRUE)
		err_R <- dd.learn[[1]]$err_out
		rm(drp.learn);rm(dd.learn)
		saveRDS(err_R, "dada2_intermediates/dada2_reverse_error.rds")
	} else {
		message("------->", date(), "	Reverse error already determined, locating from dada2_intermediates/dada2_reverse_error.rds")
		err_R<-readRDS("dada2_intermediates/dada2_reverse_error.rds")
	}
	
	
#############################

message("------->", date(), "     Dereplicating")
if (!file.exists("dada2_intermediates/dada2_derepRs.rds")){
	derepFs <- derepFastq(filtFs, verbose=TRUE)
	derepRs <- derepFastq(filtRs, verbose=TRUE)
	names(derepFs) <- sample.names
	names(derepRs) <- sample.names
		saveRDS(derepFs, "dada2_intermediates/dada2_derepFs.rds")
		saveRDS(derepFs, "dada2_intermediates/dada2_derepRs.rds")
 }else{
 	message("------->", date(),"	reading in deprelicated reads from dada2_intermediates/dada2_derep[FR]s.rds")
	derepFs<-readRDS("dada2_intermediates/dada2_derepFs.rds")
	derepRs<-readRDS("dada2_intermediates/dada2_derepRs.rds")
	}



message("------->", date(), "     DADAing forward...")
if (!file.exists("dada2_intermediates/dada2_Forwards.rds")){
	dadaFs <- dada(derepFs, err=err_F, selfConsist = TRUE, multithread=TRUE)
	saveRDS(dadaFs, "dada2_intermediates/dada2_Forwards.rds")
	}else{
	message(date(),"	Forward already denoised. Reading from dada2_intermediates/dada2_Forwards.rds")
		dadaFs<-readRDS("dada2_intermediates/dada2_Forwards.rds")
	}
	
message("------->", date(), "     DADAing reverse...")
if (!file.exists("dada2_intermediates/dada2_Reverses.rds")){
	dadaRs <- dada(derepRs, err=err_R, selfConsist = TRUE, multithread=TRUE)
	saveRDS(dadaRs, "dada2_intermediates/dada2_Reverses.rds")
	}else{
	message(date(),"	Reverse already denoised. Reading from dada2_intermediates/dada2_Reverses.rds")
		dadaRs<-readRDS("dada2_intermediates/dada2_Reverses.rds")
	}
	
	
pdf("dada2_results/ErrorProfiles.pdf")
  plotErrors(dadaFs[[1]], nominalQ=TRUE)
  plotErrors(dadaRs[[1]], nominalQ=TRUE)
dev.off()

message("------->", date(), "     Merging...")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

message("------->", date(), "     Making table..")
	seqtab <- makeSequenceTable(mergers)
message("Dimensions of seqtab: ", dim(seqtab))

message("------->", date(), "     Removing chimeras...")
	seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)
message("Fraction surviving: ", sum(seqtab.nochim)/sum(seqtab))

message("------->", date(), "     Assigning Taxonomy to species level with SILVA...")
	taxa <- assignTaxonomy(seqtab.nochim, "/netapp/home/jbisanz/dada2_training_sets/silva_nr_v123_train_set.fa.gz")
	taxa <- addSpecies(taxa, "/netapp/home/jbisanz/dada2_training_sets/silva_species_assignment_v123.fa.gz", verbose=TRUE, allowMultiple=TRUE)
	colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
write.table(taxa, file="dada2_results/ISUtaxonomy.SILVA.txt", sep='\t', quote=F, col.names=NA)

##
message("------->", date(), "     Assigning Taxonomy to species level with RDP...")
	taxa <- assignTaxonomy(seqtab.nochim, "/netapp/home/jbisanz/dada2_training_sets/rdp_train_set_14.fa.gz")
	taxa <- addSpecies(taxa, "/netapp/home/jbisanz/dada2_training_sets/rdp_species_assignment_14.fa.gz", verbose=TRUE, allowMultiple=TRUE)
	colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
write.table(taxa, file="dada2_results/ISUtaxonomy.RDP.txt", sep='\t', quote=F, col.names=NA)

##
message("------->", date(), "     Assigning Taxonomy to genus level with Green Genes...")
	taxa <- assignTaxonomy(seqtab.nochim, "/netapp/home/jbisanz/dada2_training_sets/gg_13_8_train_set_97.fa.gz")
write.table(taxa, file="dada2_results/ISUtaxonomy.GG.txt", sep='\t', quote=F, col.names=NA)

message("------->", date(), "    Complete")

EOF

