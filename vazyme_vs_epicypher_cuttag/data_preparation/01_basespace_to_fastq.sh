#!/bin/bash

# ============================================
# STEP 01: BASESPACE TO FASTQ
# PURPOSE:
#   Copy and merge gzipped FASTQ files from BaseSpace
#   for MM217 samples into a local project directory.
# INPUT:
#   BaseSpace sample folders containing R1/R2 FASTQ chunks.
# OUTPUT:
#   Merged FASTQ files in the project fastq directory.
# NOTES:
#   This step is kept commented because it depends on the
#   local storage layout and on BaseSpace being mounted.
# ============================================

# Mount BaseSpace if needed
# bs-mount BaseSpace

# Create project directories
# mkdir -p /mnt/nas-safu03/analysis/vazyme_chip_0525
# mkdir -p /mnt/nas-safu03/analysis/vazyme_chip_0525/fastq

# Copy and merge Vazyme sample R1 reads
# for i in /mnt/nas-safu03/analysis/Basespace/BSSH/Projects/New-Cut-Tag-09-24/Samples/*-MM217-* ; do
#     sample_name=${i##*/}
#     echo $sample_name
#     cat $i/Files/*R1*.gz > /mnt/nas-safu03/analysis/vazyme_chip_0525/fastq/${sample_name}_R1.fastq.gz
# done &

# Copy and merge Vazyme sample R2 reads
# for i in /mnt/nas-safu03/analysis/Basespace/BSSH/Projects/New-Cut-Tag-09-24/Samples/*-MM217-* ; do
#     sample_name=${i##*/}
#     echo $sample_name
#     cat $i/Files/*R2*.gz > /mnt/nas-safu03/analysis/vazyme_chip_0525/fastq/${sample_name}_R2.fastq.gz
# done &
