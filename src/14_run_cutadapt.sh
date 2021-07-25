#!/bin/bash

# This script performs demultiplexing and primer/adapter removal with cutadapt

cd ~/Documents/lab/implant_microbiota_study

mkdir cutadapt_reads

source /home/ccarr/anaconda3/etc/profile.d/conda.sh
conda activate cutadapt_env

# Demultiplex
cutadapt \
    -e 0.15 --no-indels -m 1 --discard-untrimmed -j 4 \
    -g file:data/fwd_barcodes.fasta \
    -G file:data/rev_barcodes.fasta \
    -o cutadapt_reads/{name1}-{name2}.1.fastq.gz -p cutadapt_reads/{name1}-{name2}.2.fastq.gz \
    Burton-SMV_S1_L001_R1_001.fastq Burton-SMV_S1_L001_R2_001.fastq

cd cutadapt_reads

mkdir trimmed_reads

# Trim
for i in *.1.fastq.gz 
do 
	sample=$(echo ${i} | sed "s/\.1\.fastq\.gz//")
	cutadapt \
		-e 0.1 -m 1 --discard-untrimmed -j 4 \
		-a GTGCCAGCMGCCGCGGTAA...ATTAGAWACCCBDGTAGTCC \
		-A GGACTACHVGGGTWTCTAAT...TTACCGCGGCKGCTGGCAC \
		-o trimmed_reads/${sample}.1.fastq.gz \
		-p trimmed_reads/${sample}.2.fastq.gz \
		${sample}.1.fastq.gz ${sample}.2.fastq.gz
done
