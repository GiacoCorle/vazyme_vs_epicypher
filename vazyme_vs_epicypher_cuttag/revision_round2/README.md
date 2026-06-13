# Revision Round 2 — Vazyme CUT&Tag vs EpiCypher ChIP-seq (H3K27ac)

This folder contains the full analysis pipeline for **Revision Round 2** of the Vazyme vs EpiCypher comparison manuscript. The pipeline assesses peak quality, method agreement, and biological enrichment for H3K27ac signal across two patient cohorts (**MM217**).

---

## Scripts

| File | Description |
|---|---|
| `vazyme_chip_1225_revision_round2.sh` | Main analysis pipeline (steps 1–10, MM217 cohort) |
| `plot_abs_distance.py` | Python helper: grouped barplot of GREAT region–gene associations by absolute distance to TSS |

---

## Pipeline steps

### [1] UpSet plots at fold-enrichment thresholds
Filters MACS2 narrowPeak files at signalValue > 3×, 4×, 5×, 10× and plots UpSet diagrams (ComplexUpset/ggplot2) showing intersection sizes across the four samples (Vazyme R1/R2, EpiCypher R1/R2).

### [2] Signal value rank plot
Percentile-ranked line plot of column 7 (fold enrichment over background) across all four samples. Used to identify a reasonable peak quality threshold.

### [3] Peak retention across -log10(q) thresholds
For peaks passing FC > 3, counts how many survive progressively stricter q-value cutoffs (−log10(q) = 3–10). Reveals that EpiCypher_R1 in MM217 loses peaks disproportionately at higher thresholds, indicating background enrichment.

### [4] UpSet at combined threshold FC > 3 & −log10(q) > 8
Applies both filters simultaneously to bring peak counts into balance across all four samples.

### [5] Merge replicates + UpSet + BED export
Union-merges each method's replicates into a single peak set, plots a 2-set UpSet (Vazyme merged vs EpiCypher merged), and exports three BED files:
- `peaks_shared.bed`
- `peaks_vazyme_only.bed`
- `peaks_epicypher_only.bed`

### [6] deepTools heatmap
Runs `computeMatrix reference-point` (±2 kb, 20 bp bins) and `plotHeatmap` on the three BED subsets using all four bigWig tracks. Colour scale: blue → yellow → red (0–2.5 CPM).

### [7] rGREAT ontology enrichment (GO:BP, MSigDB Hallmarks, MSigDB C5)
Runs local GREAT analysis on each BED subset. Saves top-10 bar plots and `.rds` objects for interactive `shinyReport()`.

### [8] rGREAT myeloma-specific enrichment (MSigDB C2)
Tests 59 published myeloma gene sets (Boylan, Chng, Corre, Davies, IRF4/Shaffer, Zhan, etc.) from MSigDB C2 curated collection against each peak subset.

### [9] Myeloma enrichment dot plot
Reads the GREAT enrichment CSV and produces a publication-quality dot plot (size = fold enrichment, colour = −log10(adj. p-value), x = peak subset, y = gene set). Filtered to gene sets with fold enrichment ≥ 1.8× in at least one group.

### [10] Combined MM217 vs MM196 dot plot
Merges enrichment tables from both cohorts into a single faceted figure for cross-patient comparison, retaining gene sets with fold enrichment ≥ 1.8× in at least one sample across either cohort.

---

## Dependencies

### R packages
```r
# Bioconductor
BiocManager::install(c("GenomicRanges", "rtracklayer", "rGREAT",
                       "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"))

# CRAN
install.packages(c("ComplexUpset", "ggplot2", "dplyr", "tidyr",
                   "scales", "ggrepel", "msigdbr"))
```

### Python (step 9 helper script)
```bash
pip install pandas numpy matplotlib
```

### Command-line tools
- [deepTools](https://deeptools.readthedocs.io) ≥ 3.5 (`computeMatrix`, `plotHeatmap`)
- [bedtools](https://bedtools.readthedocs.io) ≥ 2.30

---

## Input files required

The script expects the following files (paths are hardcoded and should be updated for your system):

| File | Description |
|---|---|
| `h3k27ac_vazyme_R1/R2.macs2_peaks.narrowPeak` | Vazyme CUT&Tag MACS2 peak calls |
| `h3k27_epicypher_R1/R2.macs2_peaks.narrowPeak` | EpiCypher ChIP-seq MACS2 peak calls |
| `h3k27ac_vazyme_R1/R2.bigWig` | Vazyme CPM-normalised signal tracks |
| `h3k27_epicypher_R1/R2.bigWig` | EpiCypher CPM-normalised signal tracks |
| `GREAT_H3K27ac_Myeloma_enrichment.csv` | GREAT MSigDB C2 myeloma enrichment table (MM217) |
| `GREAT_H3K27ac_MM196_C2_enrichment.csv` | GREAT MSigDB C2 myeloma enrichment table (MM196) |
| `great_abs_distance.csv` | Region–gene association distance distribution |

---

## Output files

All outputs are written to `revision_round_2/`:

| File | Step | Description |
|---|---|---|
| `upset_fc{3,4,5,10}.png` | 1 | UpSet plots at FC thresholds |
| `signal_rank_rev2.png` | 2 | Signal value rank plot |
| `peak_retention_qval_fc3_rev2.png` | 3 | Peak retention line plot |
| `upset_fc3_q8.png` | 4 | UpSet at FC>3 & q>8 |
| `upset_merged_replicates_fc3_q8.png` | 5 | Merged replicate UpSet |
| `peaks_{shared,vazyme_only,epicypher_only}.bed` | 5 | Peak subset BED files |
| `heatmap_fc3_q8.png` | 6 | deepTools heatmap |
| `great_top10_{GO_BP,msigdb_H,msigdb_C5}_{subset}.png` | 7 | GREAT top-10 plots |
| `great_res_{subset}_{ontology}.rds` | 7–8 | rGREAT objects for shinyReport |
| `great_top10_myeloma_C2_{subset}.png` | 8 | Myeloma gene set top-10 |
| `H3K27ac_myeloma_enrichment.png/.pdf` | 9 | Myeloma enrichment dot plot |
| `MM196_MM217_myeloma_enrichment_combined.png/.pdf` | 10 | Cross-cohort comparison dot plot |
| `abs_distance_combined.png` | 9 | TSS distance distribution barplot |

---

## Citation

If you use this pipeline, please cite the associated manuscript (in preparation).

---

*Pipeline developed by Giacomo Corleone. Generated with [Claude Code](https://claude.ai/claude-code).*
