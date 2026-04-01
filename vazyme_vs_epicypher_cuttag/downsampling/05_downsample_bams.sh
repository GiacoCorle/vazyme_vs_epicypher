#!/bin/bash

# ============================================
# STEP 05: BAM DOWNSAMPLING
# PURPOSE:
#   Generate downsampled BAM files at 10% to 90% of the
#   original read depth to evaluate the effect of read count
#   on downstream peak calling.
# INPUT:
#   Mark-duplicated sorted BAM files from the Nextflow run.
# OUTPUT:
#   Downsampled BAM files in a dedicated output directory.
# REQUIREMENTS:
#   samtools.
# NOTES:
#   The script first collates each BAM so that paired reads
#   remain adjacent before random downsampling.
# ============================================

# Move to the directory containing markdup BAM files
# cd /mnt/nas-safu03/analysis/vazyme_chip_0525/nextflow_analysis080126/02_alignment/bowtie2/target/markdup

# Define output directory for downsampled BAMs
# out="/mnt/nas-safu03/analysis/vazyme_chip_0525/nextflow_analysis080126/02_alignment/bowtie2/target/downsampled"
# mkdir -p "$out"

# Loop over BAM files
# for f in *sorted.bam; do
#     echo "Collating $f"
#
#     # Create a simplified sample name
#     base=$(echo "$f" | sed -E 's/_vazyme//; s/\.target.*//')
#
#     # Collate BAM to keep paired reads together
#     namesorted="${f%.bam}.namesorted.bam"
#     samtools collate -@ 60 -o "$namesorted" "$f"
#
#     # Generate downsampled BAMs at 10% to 90%
#     for i in .1 .2 .3 .4 .5 .6 .7 .8 .9; do
#         frac=$(echo $i | sed 's/\.//')
#         outbam="$out/${base}.ds${frac}0.bam"
#
#         echo "Downsampling $f to $i -> $outbam"
#         samtools view -@ 60 -s 42$i -b "$namesorted" > "$outbam"
#
#         # Optional QC summary
#         # samtools flagstat -@ 60 "$outbam" > "$outbam.flagstat"
#     done
#
#     # Remove temporary collated BAM
#     rm "$namesorted"
# done
