#!/bin/bash

# ============================================
# STEP 03: RUN NF-CORE/CUTANDRUN
# PURPOSE:
#   Run the main nf-core/cutandrun workflow on the samples.
# INPUT:
#   Sample manifest CSV and hg38 reference resources.
# OUTPUT:
#   Nextflow analysis directory with alignment and peak results.
# REQUIREMENTS:
#   Nextflow 24.10.0, Docker, nf-core/cutandrun v3.2.1.
# ============================================

# Install or activate the required Nextflow version
# nextflow self-update 24.10.0 || true
# nextflow -v || true

# Fix Windows line endings in the manifest file if needed
# sed -i 's/\r$//' ~/scripts/manifest_vazume_epicypher_cuttag_2026.csv

# Run nf-core/cutandrun with the requested version
# NXF_VER=24.10.0 nextflow run nf-core/cutandrun \
#     --input ~/scripts/manifest_vazume_epicypher_cuttag_2026.csv \
#     --outdir /mnt/nas-safu03/analysis/vazyme_chip_0525/nextflow_analysis080126 \
#     --genome hg38 \
#     --blacklist /mnt/nas-safu03/analysis/vazyme_chip_0525/hg38-blacklist.v2.bed \
#     -r 3.2.1 \
#     --use_control 'false' \
#     -profile docker \
#     --normalisation_mode CPM \
#     --peakcaller macs2
