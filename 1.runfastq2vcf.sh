#!/bin/bash

## HOW TO USE:

# bash ./runfastq2vcfvcf.sh MyReferenceGenome.fa

# 2>&1 means "copy standard errors to standard output", and the > redirects the standard output (now containing error if there are any) to the file fastq2vcf_log.txt


# Check if the user provided a reference genome file
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <reference_genome_file>"
    exit 1
fi

# Assign the first command-line argument to a variable
REF_GENOME=$1

# Check if the reference genome file exists
if [ ! -f "$REF_GENOME" ]; then
    echo "Error: Reference genome file '$REF_GENOME' does not exist."
    exit 1
fi

# Log information
echo "Using reference genome file: $REF_GENOME"

# 2>&1 means "copy standard errors to standard output", and the > redirects the standard output (now containing error if there are any) to the file fastq2vcf_log.txt

ls *_1.fastq.gz | sed 's/_1.fastq.gz//' > fastq2vcfsamples.txt

cat fastq2vcfsamples.txt | parallel -j 10 \
"fastq2vcf.py all --read1 {}_1.fastq.gz --read2 {}_2.fastq.gz \
--ref $REF_GENOME \
--prefix {}" > fastq2vcf_log.txt 2>&1
