#!/bin/bash

# ============================================
# STEP 06: MACS2 PEAK CALLING
# PURPOSE:
#   Call peaks on the downsampled BAM files and generate
#   coverage tracks for visualization.
# INPUT:
#   Downsampled BAM files.
# OUTPUT:
#   MACS2 peak calls, pileup bedGraph files, and BigWig files.
# REQUIREMENTS:
#   Conda environment with MACS2, bedGraphToBigWig.
# NOTES:
#   The peak-calling loop is kept commented to preserve the
#   original execution history without changing behavior.
# ============================================

source ~/miniconda3/bin/activate
conda activate macs2

bams="/mnt/nas-safu03/analysis/vazyme_chip_0525/nextflow_analysis080126/02_alignment/bowtie2/target/downsampled/*.bam"
out_peak="/mnt/nas-safu03/analysis/vazyme_chip_0525/nextflow_analysis080126/03_peak_calling/macs2/downsampled"
tempdir="/data/tmpdir/"

# Loop over downsampled BAM files
# count_1=0
# for filenames in $bams ; do
#     filename2=${filenames##*/}
#     filename3=${filename2%.*}
#     let "count_1++"
#     echo $count_1
#     echo $filenames
#     echo $filename3
#
#     macs2 callpeak \
#         --nomodel \
#         --shift -75 \
#         --extsize 150 \
#         --keep-dup all \
#         -q 0.01 \
#         --gsize 2.7E+9 \
#         --format BAMPE \
#         --treatment $filenames \
#         --tempdir $tempdir \
#         --outdir $out_peak \
#         -n $filename3
#
#     sort -k1,1 -k2,2n $out_peak/"$filename3"_treat_pileup.bdg > $out_peak/"$filename3"_treat_pileup_S.bdg
#     bedGraphToBigWig $out_peak/"$filename3"_treat_pileup_S.bdg /data/genomes/hg38_UCSC/hg38.chrom.sizes $out_peak/"$filename3"_treat_pileup_S.bw
#     echo $count_1 "done"
#     cd $out_peak
# done
