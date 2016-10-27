#!/bin/bash
#$ -S /bin/bash
#$ -o qiime.outputlog.txt
#$ -e qiime.errorlog.txt
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=64G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=100G
#$ -l h_rt=336:0:0

#Using QIIME 1.9.1, on UCSF QB3 cluster under OGS/GE 2011.11p1


#################################################################
#JB Oct25/2016
# Version 1.3
# Will merge the forward and reverse reads  from illumina sequencers, demultiplex and pick OTUs. I only did a base install on the cluster so after completion move files to biggut.ucsf.edu, and download biom files and trees and use macqiime to analyze on own computer
# All work will be done in what ever is described in WORKING_DIR so make sure your reads are here
# Expecting V4 reads with Caparaso golay primers
#
#
#Change log:
# Oct25, 2016- Minor cosmetic changes.
# Sept20, 2016- Adding Chimera removal with usearch61, and changing overlap to vsearch which requires separate demultiplexing
# Mar3, 2016- Reduced the number of jobs being submitted as most weren't clearing the que before the others were finished, reduced request for netapp as this shouldn't be used for this analysis
#
# July6, 2016- Added the SILVA DB as an option, also changed the method of counting reads as it was taking obscenely long, also added a default taxa summary plot and pcoa as standard outputs
#			- Now also putting in a minimum abundance filter, an OTU needs at least 10 counts in at least 2 samples to be retained, this can be modified below
#			- Also trimming back the seed sequences and tree to only otus passing filter
#			- Be aware that SILVA is using the 90% majority taxonomy assignments and not the 100% consensus which could result in a slight over calling of taxa ids, you could change this by creating your own version of settings/qiime_config_silva.txt and adjusting line 125
#			- Reduced number of jobs to 14 such that it will fill our allocated number of lab.q slots (14 total)
#			- Total run time for 16million 251x151 reads timed at 1.5h with SILVA, closer to 1h with GG.
#
#################################################################

echo "##################################################################################"
echo "$(date)	Running QIIME Pipeline v1.3 on node $(hostname)"

#INDICATE YOUR FORWARD/REVERSE/INDEX READS AND MAPPING FILE
#PREPARE YOUR MAPPING FILE LOCALLY WITH MACQIIME (validate_mapping_file.py)
#PUT YOUR INFORMATION HERE

WORKING_DIR=/scrapp2/INSERTHERE
FORWARD_READ=$WORKING_DIR/INSERTHERE
REVERSE_READ=$WORKING_DIR/INSERTHERE
INDEX_READ=$WORKING_DIR/INSERTHERE
MAPPING_FILE=$WORKING_DIR/INSERTHERE

DBCHOICE="GG" #Use "SILVA" to use SILVA123 or "GG" for Green Genes May 2013

MINREADS=10 #For an OTU to be retained, needs to be at counted at least 10x
MINSAMPLES=2  #and present in at least 2 samples


################################################################
#DON'T CHANGE WHAT IS BELOW HERE


source /netapp/home/jbisanz/.bash_profile
source /netapp/home/jbisanz/qiime_for_cluster_1.9.1/settings/sourceme.sh 

if [ "$DBCHOICE" = "SILVA" ]
then
	echo "Using SILVA 123"
	export QIIME_CONFIG_FP=/netapp/home/jbisanz/qiime_for_cluster_1.9.1/settings/qiime_config_silva.txt
else
	echo "Using Green Genes"
fi

echo "$(date)	Validating Mapping File"
	validate_mapping_file.py -m $MAPPING_FILE -o mapping_validation
echo "...done"

echo "$(date)     Demultiplexing $FORWARD_READ"
if [ ! -e $WORKING_DIR/splitlib/forward_seqs.fastq ]; then
split_libraries_fastq.py -m $MAPPING_FILE --rev_comp_mapping_barcodes --barcode_type golay_12 -r 999 -n 999 -q 0 -p 0.0001 -i $FORWARD_READ -b $INDEX_READ -o $WORKING_DIR/splitlib --store_demultiplexed_fastq
	mv $WORKING_DIR/splitlib/seqs.fastq $WORKING_DIR/splitlib/forward_seqs.fastq
	rm $WORKING_DIR/splitlib/seqs.fna
	mv $WORKING_DIR/splitlib/split_library_log.txt $WORKING_DIR/splitlib/forward_split_library_log.txt
else
	echo "...already done, skipping"
fi

echo "$(date)     Demultiplexing $REVERSE_READ"
if [ ! -e $WORKING_DIR/splitlib/reverse_seqs.fastq ]; then
split_libraries_fastq.py -m $MAPPING_FILE --rev_comp_mapping_barcodes --barcode_type golay_12 -r 999 -n 999 -q 0 -p 0.0001 -i $REVERSE_READ -b $INDEX_READ -o $WORKING_DIR/splitlib --store_demultiplexed_fastq
	mv $WORKING_DIR/splitlib/seqs.fastq $WORKING_DIR/splitlib/reverse_seqs.fastq
	rm $WORKING_DIR/splitlib/seqs.fna
	mv $WORKING_DIR/splitlib/split_library_log.txt $WORKING_DIR/splitlib/reverse_split_library_log.txt
else
        echo "...already done, skipping"
fi

echo "$(date)     Summary of demultiplexing:"
head -n 12  $WORKING_DIR/splitlib/forward_split_library_log.txt

echo "$(date)     Overlapping"
if [ ! -e $WORKING_DIR/splitlib/seqs.fna ]; then
	vsearch --fastq_mergepairs $WORKING_DIR/splitlib/forward_seqs.fastq  --fastaout $WORKING_DIR/splitlib/seqs.fna --fastqout $WORKING_DIR/splitlib/seqs.fastq --fastq_maxmergelen 260  --fastq_minmergelen 240 --fastq_maxns 0 --reverse $WORKING_DIR/splitlib/reverse_seqs.fastq
else
        echo "...already done, skipping"
fi

echo "$(date)	Removing Chimeras with USEARCH 6.1"
if [ ! -e $WORKING_DIR/splitlib/seqs.nochim.fna ]; then
	identify_chimeric_seqs.py -i $WORKING_DIR/splitlib/seqs.fna -o $WORKING_DIR/splitlib/ -m usearch61 -r /netapp/home/jbisanz/qiime_for_cluster_1.9.1/dependencies/gg_13_8_otus/rep_set/97_otus.fasta
	filter_fasta.py -f $WORKING_DIR/splitlib/seqs.fna -o $WORKING_DIR/splitlib/seqs.nochim.fna -s $WORKING_DIR/splitlib/chimeras.txt -n	
else
        echo "...already done, skipping"
fi


echo "$(date) Summary of demultiplexing, overlapping and chimera removal:"
declare -i SEQS_INPUT
declare -i SEQS_F
declare -i SEQS_R
declare -i SEQS_OVERLAP
declare -i SEQS_CHIMERA

	SEQS_INPUT=$(zcat $INDEX_READ | wc -l)/4
	SEQS_F=$(cat $WORKING_DIR/splitlib/forward_seqs.fastq | wc -l)/4
	SEQS_R=$(cat $WORKING_DIR/splitlib/reverse_seqs.fastq | wc -l)/4
	SEQS_OVERLAP=$(cat $WORKING_DIR/splitlib/seqs.fastq | wc -l)/4
	SEQS_CHIMERA=$(cat $WORKING_DIR/splitlib/seqs.nochim.fna | wc -l)/2

echo "$SEQS_INPUT	#Input Sequences"
echo "$SEQS_F	#Demultiplexed Forward Sequences"
echo "$SEQS_R	#Demultiplexed Reverse Sequences"
echo "$SEQS_OVERLAP	#Overlapped sequences"
echo "$SEQS_CHIMERA	#Post-chimera filtering sequences"


echo "$(date)	Starting Open Reference Picking"
if [ "$DBCHOICE" = "GG" ]
then
	echo "Using the Green Genes lane mask"
	pick_open_reference_otus.py -o $WORKING_DIR/otus/ -i $WORKING_DIR/splitlib/seqs.nochim.fna -p /netapp/home/jbisanz/qiime_for_cluster_1.9.1/settings/otu_picking_prams.txt  -aO 14
else
	echo "Using the entropy (0.10) and gap fraction (0.8) instead of lane mask for filter_alignments.py as reccomended in Silva_123_notes.txt"
	pick_open_reference_otus.py -o $WORKING_DIR/otus/ -i $WORKING_DIR/splitlib/seqs.fna -p /netapp/home/jbisanz/qiime_for_cluster_1.9.1/settings/otu_picking_prams_silva.txt  -aO 14
fi
echo "...done picking"


echo "$(date)	Cleaning up"
mkdir paralleljob_logs
mv **.o[0-9]* paralleljob_logs


echo "Doing default trimming as specified above: At least $MINREADS Reads in at least $MINSAMPLES Samples. Originals have been retained in otus/"
	filter_otus_from_otu_table.py -n $MINREADS -s $MINSAMPLES -i otus/otu_table_mc2_w_tax_no_pynast_failures.biom -o master_otu_table.biom
	filter_fasta.py -f otus/rep_set.fna -b master_otu_table.biom -o otu_seed_seqs.fna
	filter_tree.py -i otus/rep_set.tre -f otu_seed_seqs.fna -o otu_seed_seqs.tree


echo "################Some output statistics for your reference"
NEW=$(grep -c '>New' otu_seed_seqs.fna)
#GG=$(grep -c '>[0-9]' otu_seed_seqs.fna)
biom summarize-table -i master_otu_table.biom -o summary.master_otu_table.biom.txt
SAMPSUM=$(head -n 11 summary.master_otu_table.biom.txt)
echo "################OTU summary:"
count_seqs.py -i otu_seed_seqs.fna


echo "For your reference: there are $NEW de novo OTUs..."
echo "################Biom summary:"
echo "$SAMPSUM"

echo "Doing preliminary Taxa Summaries, find results in Taxa_Summaries, NOTE: Plots will not be generated"
summarize_taxa.py -i master_otu_table.biom -L 2,3,4,5,6,7 -a -o Taxa_Summaries

echo "Doing preliminary PCoA subsampled to 10000, find results in PCoAs"
echo "beta_diversity:metrics  bray_curtis,unweighted_unifrac,weighted_unifrac" >>bdiv_prams.txt
beta_diversity_through_plots.py -i master_otu_table.biom -m $MAPPING_FILE -o PCoAs -e 10000 -t otu_seed_seqs.tree -p bdiv_prams.txt


echo "#####################################################################"
echo "Pipeline complete at: $(date)"
echo "Move folder back to biggut for futher analysis! You have ~ 2 weeks till it will be automatically erased!"
echo "Try the following command while ssh'd into biggut: scp -r yourusername@pass1.compbio.ucsf.edu:/scrapp/yourfolder /where ever you want it"


