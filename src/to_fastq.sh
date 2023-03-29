#!/bin/bash
# to_fastq.sh
#
# This script produces FASTQ files from the reads of a BAM file. 

while getopts "i:" opt; do
  case $opt in
    i) input="$OPTARG" ;;
    *) echo "Invalid option -$OPTARG" >&2 ;;
  esac
done

if [ -z "$input" ]; then
  echo "Usage: to_fastq.sh -i <input.bam>"
  exit 1
fi

# Remove everything after the first dot (path/to/file.bam -> path/to/file)
path=$(echo ${input} | sed 's/\..*//')

# Obtaining the FASTQ files from the BAM file
samtools sort -n ${input} -o ${path}.nsorted.bam
samtools fastq ${path}.nsorted.bam -1 ${path}_R1.fastq -2 ${path}_R2.fastq -n  # -0 /dev/null -s /dev/null  # Discard unpaired reads 
rm ${path}.nsorted.bam