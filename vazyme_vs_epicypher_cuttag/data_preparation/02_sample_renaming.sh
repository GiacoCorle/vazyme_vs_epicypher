#!/bin/bash

# ============================================
# STEP 02: SAMPLE RENAMING
# PURPOSE:
#   Harmonize file names for Epicypher and Vazyme samples.
# INPUT:
#   FASTQ files in the epicypher subdirectory.
# OUTPUT:
#   Renamed files with consistent kit labels.
# NOTES:
#   This step is kept commented because it is intended to
#   be run only after manual inspection of the input files.
# ============================================

# Create Epicypher subdirectory and move into it
# mkdir -p /mnt/nas-safu03/analysis/vazyme_chip_0525/fastq/epicypher
# cd /mnt/nas-safu03/analysis/vazyme_chip_0525/fastq/epicypher

# Replace ChIP with epicypher in file names
# for f in *ChIP*; do
#     mv -n -- "$f" "${f//ChIP/epicypher}"
# done

# Replace CutTag with vazyme in file names
# for f in *CutTag*; do
#     mv -n -- "$f" "${f//CutTag/vazyme}"
# done
