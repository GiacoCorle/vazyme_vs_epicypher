#!/bin/bash

# ============================================
# STEP 07: BLACKLIST FILTERING
# PURPOSE:
#   Remove peaks overlapping ENCODE blacklist regions from
#   the MACS2 narrowPeak output files.
# INPUT:
#   MACS2 narrowPeak files.
# OUTPUT:
#   Filtered BED files prefixed with blacklist_.
# REQUIREMENTS:
#   bedtools and hg38 blacklist BED file.
# ============================================

out_peak="/mnt/nas-safu03/analysis/vazyme_chip_0525/nextflow_analysis080126/03_peak_calling/macs2/downsampled"

for filenames in $out_peak/*.narrowPeak; do
    filename2=${filenames##*/}
    filename3=${filename2%.*}
    bedtools intersect -v -a $filenames -b /data/genomes/hg38_UCSC/hg38-blacklist.v3.bed > $out_peak/blacklist_$filename3.bed
done
