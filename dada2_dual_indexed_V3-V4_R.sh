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
#$ -pe smp 16

echo "##################################################################################"
echo "$(date)	Running DADA2 Pipeline v1.0 for demultiplexed sequences from Basespace on node $(hostname)"
echo "$(date)	Assuming overlapping V3-V4 reads with primer sequences still embedded"

#Enter your settings below which will be exported to the environment for import to R
#Pipeline has been updated to version 1.5.2 of dada2 installed from github

export PROJECTNAME="test"
export WORKDIR="~/testing/16S_pipeline/"             #"~/Research/MHGEP/MGHEP_16S/" #working directory
export READS="reads/"                               #name of read folder from Basespace folder
export SAMPLESHEET="sample_sheet_trunc.csv"          #illumina sample sheet from sequencer, the result of MakeSampleSheet.R
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


Rscript dada2_dual_indexed_V3-V4_R.R


