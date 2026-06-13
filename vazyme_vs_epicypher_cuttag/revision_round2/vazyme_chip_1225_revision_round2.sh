#!/bin/bash

# =============================================================================
# vazyme_chip_1225_revision_round2.sh
# Revision Round 2 analyses — Vazyme vs EpiCypher CUT&Tag / ChIP-seq
# =============================================================================

# --------------------------------------------------------------------------
# [1] UpSet plot of significant peaks across the four samples
# --------------------------------------------------------------------------

mkdir -p /Users/giacomocorleone/Desktop/Projects/vazyme_vs_epicypher/revision_round_2

Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ComplexUpset)
  library(ggplot2)
})

peaks_dir <- "/Users/giacomocorleone/Desktop/Projects/vazyme_vs_epicypher/peaks_hg38"
out_dir   <- "/Users/giacomocorleone/Desktop/Projects/vazyme_vs_epicypher/revision_round_2"

extra_cols <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
peaks_data <- list(
  Vazyme_R1    = import(file.path(peaks_dir, "h3k27ac_vazyme_R1.macs2_peaks.narrowPeak"),  format = "BED", extraCols = extra_cols),
  Vazyme_R2    = import(file.path(peaks_dir, "h3k27ac_vazyme_R2.macs2_peaks.narrowPeak"),  format = "BED", extraCols = extra_cols),
  EpiCypher_R1 = import(file.path(peaks_dir, "h3k27_epicypher_R1.macs2_peaks.narrowPeak"), format = "BED", extraCols = extra_cols),
  EpiCypher_R2 = import(file.path(peaks_dir, "h3k27_epicypher_R2.macs2_peaks.narrowPeak"), format = "BED", extraCols = extra_cols)
)

sample_names <- c("Vazyme_R1", "Vazyme_R2", "EpiCypher_R1", "EpiCypher_R2")

bar_colors <- c(
  Vazyme_R1    = "#1b9e77",
  Vazyme_R2    = "#1b9e77",
  EpiCypher_R1 = "royalblue3",
  EpiCypher_R2 = "royalblue3"
)

make_mat <- function(gr_list) {
  combined  <- suppressWarnings(Reduce(c, unname(gr_list)))
  all_peaks <- reduce(combined)
  as.data.frame(sapply(gr_list, function(gr)
    suppressWarnings(countOverlaps(all_peaks, gr)) > 0
  ))
}

thresholds <- list(
  list(fc = 3,  label = "Signal > 3x",  file = "upset_fc3.png"),
  list(fc = 4,  label = "Signal > 4x",  file = "upset_fc4.png"),
  list(fc = 5,  label = "Signal > 5x",  file = "upset_fc5.png"),
  list(fc = 10, label = "Signal > 10x", file = "upset_fc10.png")
)

for (th in thresholds) {
  filtered <- lapply(peaks_data, function(gr) gr[gr$signalValue > th$fc])
  mat      <- make_mat(filtered)

  p <- upset(
    mat,
    intersect      = sample_names,
    sort_sets      = FALSE,
    sort_intersections_by = "cardinality",
    base_annotations = list(
      "Intersection size" = intersection_size(
        counts = TRUE,
        text   = list(size = 3),
        bar_number_threshold = 1
      ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        theme(axis.title.y = element_text(face = "bold"))
    ),
    set_sizes = (
      upset_set_size(
        filter_intersections = FALSE,
        geom = geom_bar(aes(fill = group), width = 0.6)
      ) +
        scale_fill_manual(values = bar_colors, guide = "none") +
        theme(axis.title.x = element_text(face = "bold"))
    ),
    matrix = (
      intersection_matrix(
        geom    = geom_point(size = 3),
        segment = geom_segment(linewidth = 0.8)
      ) +
        scale_color_manual(
          values = c("TRUE" = "grey25", "FALSE" = "grey85"),
          guide  = "none"
        )
    ),
    themes = upset_modify_themes(list(
      "overall_sizes"        = theme(text = element_text(size = 12)),
      "intersections_matrix" = theme(text = element_text(size = 12))
    ))
  ) +
    labs(title = th$label) +
    theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))

  out_file <- file.path(out_dir, th$file)
  ggsave(out_file, plot = p, width = 8, height = 5.5, units = "in", dpi = 300, bg = "white")
  cat("Saved:", out_file, "\n")
}
EOF

# --------------------------------------------------------------------------
# [2] Signal value rank plot (col 7 — fold enrichment over background)
#     All 4 samples plotted as percentile-ranked line plot to assess the
#     distribution of enrichment and identify a reasonable peak threshold.
# --------------------------------------------------------------------------

Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(rtracklayer)
  library(ggplot2)
  library(dplyr)
})

peaks_dir  <- "/Users/giacomocorleone/Desktop/Projects/vazyme_vs_epicypher/peaks_hg38"
out_dir    <- "/Users/giacomocorleone/Desktop/Projects/vazyme_vs_epicypher/revision_round_2"
extra_cols <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")

files <- list(
  Vazyme_R1    = file.path(peaks_dir, "h3k27ac_vazyme_R1.macs2_peaks.narrowPeak"),
  Vazyme_R2    = file.path(peaks_dir, "h3k27ac_vazyme_R2.macs2_peaks.narrowPeak"),
  EpiCypher_R1 = file.path(peaks_dir, "h3k27_epicypher_R1.macs2_peaks.narrowPeak"),
  EpiCypher_R2 = file.path(peaks_dir, "h3k27_epicypher_R2.macs2_peaks.narrowPeak")
)

df <- bind_rows(lapply(names(files), function(s) {
  vals <- sort(import(files[[s]], format = "BED", extraCols = extra_cols)$signalValue)
  data.frame(
    sample    = s,
    pct       = 100 * seq_along(vals) / length(vals),
    signal    = vals,
    method    = ifelse(grepl("Vazyme", s), "Vazyme", "EpiCypher"),
    replicate = ifelse(grepl("R1", s), "R1", "R2")
  )
}))

p <- ggplot(df, aes(x = pct, y = signal, color = method, linetype = replicate)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("Vazyme" = "#1b9e77", "EpiCypher" = "royalblue3")) +
  scale_x_continuous(breaks = seq(0, 100, 20), labels = paste0(seq(0, 100, 20), "%")) +
  labs(
    x = "Peak rank (percentile)", y = "Signal value (fold enrichment over background)",
    color = "Method", linetype = "Replicate"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"), axis.title = element_text(face = "bold"),
    legend.position = "right", panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(out_dir, "signal_rank_rev2.png"), plot = p,
       width = 7, height = 4.5, units = "in", dpi = 300, bg = "white")
cat("Saved: signal_rank_rev2.png\n")
EOF

# --------------------------------------------------------------------------
# [3] Peak retention table and line plot across -log10(q) thresholds
#     For peaks passing FC > 3, count how many survive each -log10(q) cutoff
#     from 3 to 10. EpiCypher_R1 drops disproportionately, confirming that
#     many of its peaks are weakly supported background calls.
#     Convergence point (~q > 9–10) indicates the threshold where all samples
#     reach comparable peak numbers.
# --------------------------------------------------------------------------

Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(rtracklayer)
  library(ggplot2)
  library(dplyr)
})

peaks_dir  <- "/Users/giacomocorleone/Desktop/Projects/vazyme_vs_epicypher/peaks_hg38"
out_dir    <- "/Users/giacomocorleone/Desktop/Projects/vazyme_vs_epicypher/revision_round_2"
extra_cols <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")

files <- list(
  Vazyme_R1    = file.path(peaks_dir, "h3k27ac_vazyme_R1.macs2_peaks.narrowPeak"),
  Vazyme_R2    = file.path(peaks_dir, "h3k27ac_vazyme_R2.macs2_peaks.narrowPeak"),
  EpiCypher_R1 = file.path(peaks_dir, "h3k27_epicypher_R1.macs2_peaks.narrowPeak"),
  EpiCypher_R2 = file.path(peaks_dir, "h3k27_epicypher_R2.macs2_peaks.narrowPeak")
)

thresholds <- 3:10

# Print summary table
rows <- lapply(names(files), function(s) {
  gr <- import(files[[s]], format = "BED", extraCols = extra_cols)
  gr <- gr[gr$signalValue > 3]
  n  <- length(gr)
  counts <- sapply(thresholds, function(t) sum(gr$qValue > t))
  pcts   <- round(100 * counts / n, 1)
  row    <- as.data.frame(t(c(Peaks_FC_gt3 = n, counts)))
  colnames(row) <- c("Peaks_FC_gt3", paste0("q>", thresholds))
  pct_row <- as.data.frame(t(pcts))
  colnames(pct_row) <- paste0("pct_q>", thresholds)
  cbind(Sample = s, row, pct_row)
})
df_table <- bind_rows(rows)
cat("\nCounts:\n"); print(df_table[, c("Sample", "Peaks_FC_gt3", paste0("q>", thresholds))], row.names = FALSE)
cat("\nPercentages:\n"); print(df_table[, c("Sample", paste0("pct_q>", thresholds))], row.names = FALSE)

# Line plot of peak counts across thresholds
df_plot <- bind_rows(lapply(names(files), function(s) {
  gr <- import(files[[s]], format = "BED", extraCols = extra_cols)
  gr <- gr[gr$signalValue > 3]
  data.frame(
    sample    = s,
    threshold = thresholds,
    count     = sapply(thresholds, function(t) sum(gr$qValue > t)),
    method    = ifelse(grepl("Vazyme", s), "Vazyme", "EpiCypher"),
    replicate = ifelse(grepl("R1", s), "R1", "R2")
  )
}))

p <- ggplot(df_plot, aes(x = threshold, y = count, color = method, linetype = replicate)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Vazyme" = "#1b9e77", "EpiCypher" = "royalblue3")) +
  scale_linetype_manual(values = c("R1" = "solid", "R2" = "dashed")) +
  scale_x_continuous(breaks = thresholds, labels = paste0(">", thresholds)) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = expression(-log[10]~"(q-value) threshold"), y = "Number of peaks",
    color = "Method", linetype = "Replicate",
    title = "Peak retention across -log10(q) thresholds (FC > 3)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"), axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    legend.position = "right", panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(out_dir, "peak_retention_qval_fc3_rev2.png"), plot = p,
       width = 7, height = 4.5, units = "in", dpi = 300, bg = "white")
cat("Saved: peak_retention_qval_fc3_rev2.png\n")
EOF

# --------------------------------------------------------------------------
# [4] UpSet plot at combined threshold: FC > 3 & -log10(q) > 8
#     At this combined threshold the four samples reach comparable peak counts
#     (~35k–52k), removing the excess weakly-enriched peaks in EpiCypher_R1.
# --------------------------------------------------------------------------

Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ComplexUpset)
  library(ggplot2)
})

peaks_dir  <- "/Users/giacomocorleone/Desktop/Projects/vazyme_vs_epicypher/peaks_hg38"
out_dir    <- "/Users/giacomocorleone/Desktop/Projects/vazyme_vs_epicypher/revision_round_2"
extra_cols <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")

peaks_data <- list(
  Vazyme_R1    = import(file.path(peaks_dir, "h3k27ac_vazyme_R1.macs2_peaks.narrowPeak"),  format = "BED", extraCols = extra_cols),
  Vazyme_R2    = import(file.path(peaks_dir, "h3k27ac_vazyme_R2.macs2_peaks.narrowPeak"),  format = "BED", extraCols = extra_cols),
  EpiCypher_R1 = import(file.path(peaks_dir, "h3k27_epicypher_R1.macs2_peaks.narrowPeak"), format = "BED", extraCols = extra_cols),
  EpiCypher_R2 = import(file.path(peaks_dir, "h3k27_epicypher_R2.macs2_peaks.narrowPeak"), format = "BED", extraCols = extra_cols)
)

sample_names <- c("Vazyme_R1", "Vazyme_R2", "EpiCypher_R1", "EpiCypher_R2")
bar_colors   <- c(Vazyme_R1 = "#1b9e77", Vazyme_R2 = "#1b9e77", EpiCypher_R1 = "royalblue3", EpiCypher_R2 = "royalblue3")

filtered  <- lapply(peaks_data, function(gr) gr[gr$signalValue > 3 & gr$qValue > 8])
combined  <- suppressWarnings(Reduce(c, unname(filtered)))
all_peaks <- reduce(combined)
mat <- as.data.frame(sapply(filtered, function(gr)
  suppressWarnings(countOverlaps(all_peaks, gr)) > 0
))

p <- upset(
  mat, intersect = sample_names, sort_sets = FALSE, sort_intersections_by = "cardinality",
  base_annotations = list(
    "Intersection size" = intersection_size(counts = TRUE, text = list(size = 3), bar_number_threshold = 1) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      theme(axis.title.y = element_text(face = "bold"))
  ),
  set_sizes = (
    upset_set_size(filter_intersections = FALSE, geom = geom_bar(aes(fill = group), width = 0.6)) +
      scale_fill_manual(values = bar_colors, guide = "none") +
      theme(axis.title.x = element_text(face = "bold"))
  ),
  matrix = (
    intersection_matrix(geom = geom_point(size = 3), segment = geom_segment(linewidth = 0.8)) +
      scale_color_manual(values = c("TRUE" = "grey25", "FALSE" = "grey85"), guide = "none")
  ),
  themes = upset_modify_themes(list(
    "overall_sizes" = theme(text = element_text(size = 12)),
    "intersections_matrix" = theme(text = element_text(size = 12))
  ))
) +
  labs(title = "FC > 3  &  -log10(q) > 8") +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))

ggsave(file.path(out_dir, "upset_fc3_q8.png"), plot = p,
       width = 8, height = 5.5, units = "in", dpi = 300, bg = "white")
cat("Saved: upset_fc3_q8.png\n")
EOF

# --------------------------------------------------------------------------
# [5] Merge replicates per method and UpSet at FC > 3 & -log10(q) > 8
#     Vazyme R1+R2 and EpiCypher R1+R2 are each collapsed into a single
#     non-redundant peak set (union + reduce). The UpSet shows shared vs.
#     method-specific peaks between the two merged sets.
#     After filtering, both methods yield ~60k peaks, confirming balance.
#     The 3 resulting BED files (shared, Vazyme only, EpiCypher only) are
#     saved for downstream deepTools analysis.
# --------------------------------------------------------------------------

Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ComplexUpset)
  library(ggplot2)
})

peaks_dir  <- "/Users/giacomocorleone/Desktop/Projects/vazyme_vs_epicypher/peaks_hg38"
out_dir    <- "/Users/giacomocorleone/Desktop/Projects/vazyme_vs_epicypher/revision_round_2"
extra_cols <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")

load_filter <- function(f) {
  gr <- import(f, format = "BED", extraCols = extra_cols)
  gr[gr$signalValue > 3 & gr$qValue > 8]
}

vazyme    <- suppressWarnings(reduce(c(
  load_filter(file.path(peaks_dir, "h3k27ac_vazyme_R1.macs2_peaks.narrowPeak")),
  load_filter(file.path(peaks_dir, "h3k27ac_vazyme_R2.macs2_peaks.narrowPeak"))
)))
epicypher <- suppressWarnings(reduce(c(
  load_filter(file.path(peaks_dir, "h3k27_epicypher_R1.macs2_peaks.narrowPeak")),
  load_filter(file.path(peaks_dir, "h3k27_epicypher_R2.macs2_peaks.narrowPeak"))
)))

all_peaks <- suppressWarnings(reduce(c(vazyme, epicypher)))
mat <- data.frame(
  Vazyme    = suppressWarnings(countOverlaps(all_peaks, vazyme))    > 0,
  EpiCypher = suppressWarnings(countOverlaps(all_peaks, epicypher)) > 0
)

bar_colors <- c(Vazyme = "#1b9e77", EpiCypher = "royalblue3")

p <- upset(
  mat, intersect = c("Vazyme", "EpiCypher"), sort_sets = FALSE,
  sort_intersections_by = "cardinality",
  base_annotations = list(
    "Intersection size" = intersection_size(counts = TRUE, text = list(size = 4), bar_number_threshold = 1) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      theme(axis.title.y = element_text(face = "bold"))
  ),
  set_sizes = (
    upset_set_size(filter_intersections = FALSE, geom = geom_bar(aes(fill = group), width = 0.5)) +
      scale_fill_manual(values = bar_colors, guide = "none") +
      theme(axis.title.x = element_text(face = "bold"))
  ),
  matrix = (
    intersection_matrix(geom = geom_point(size = 4), segment = geom_segment(linewidth = 1)) +
      scale_color_manual(values = c("TRUE" = "grey25", "FALSE" = "grey85"), guide = "none")
  ),
  themes = upset_modify_themes(list(
    "overall_sizes" = theme(text = element_text(size = 13)),
    "intersections_matrix" = theme(text = element_text(size = 13))
  ))
) +
  labs(title = "Merged replicates — FC > 3  &  -log10(q) > 8") +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))

ggsave(file.path(out_dir, "upset_merged_replicates_fc3_q8.png"), plot = p,
       width = 7, height = 5, units = "in", dpi = 300, bg = "white")
cat("Saved: upset_merged_replicates_fc3_q8.png\n")

# Export 3 BED subsets
shared         <- suppressWarnings(subsetByOverlaps(vazyme,    epicypher))
vazyme_only    <- suppressWarnings(subsetByOverlaps(vazyme,    epicypher, invert = TRUE))
epicypher_only <- suppressWarnings(subsetByOverlaps(epicypher, vazyme,    invert = TRUE))

cat("Shared:", length(shared), "| Vazyme only:", length(vazyme_only),
    "| EpiCypher only:", length(epicypher_only), "\n")

export(shared,         file.path(out_dir, "peaks_shared.bed"),         format = "BED")
export(vazyme_only,    file.path(out_dir, "peaks_vazyme_only.bed"),    format = "BED")
export(epicypher_only, file.path(out_dir, "peaks_epicypher_only.bed"), format = "BED")
cat("BED files saved: peaks_shared.bed, peaks_vazyme_only.bed, peaks_epicypher_only.bed\n")
EOF

# --------------------------------------------------------------------------
# [6] deepTools heatmap on the 3 peak subsets (shared / Vazyme only /
#     EpiCypher only) using all 4 bigWig files.
#     Signal computed ±2 kb around peak centres (bin 20 bp).
#     Colour scale: blue (0) → yellow (1.25) → red (2.5).
# --------------------------------------------------------------------------

source ~/miniconda3/bin/activate

out_dir="/Users/giacomocorleone/Desktop/Projects/vazyme_vs_epicypher/revision_round_2"
bw_dir="/Users/giacomocorleone/Desktop/Projects/vazyme_vs_epicypher/peaks_hg38"

computeMatrix reference-point \
  --referencePoint center \
  --scoreFileName \
    "$bw_dir/h3k27ac_vazyme_R1.bigWig" \
    "$bw_dir/h3k27ac_vazyme_R2.bigWig" \
    "$bw_dir/h3k27_epicypher_R1.bigWig" \
    "$bw_dir/h3k27_epicypher_R2.bigWig" \
  --regionsFileName \
    "$out_dir/peaks_shared.bed" \
    "$out_dir/peaks_vazyme_only.bed" \
    "$out_dir/peaks_epicypher_only.bed" \
  --beforeRegionStartLength 2000 \
  --afterRegionStartLength  2000 \
  --binSize 20 \
  --missingDataAsZero \
  --numberOfProcessors 8 \
  --outFileName "$out_dir/matrix_heatmap.gz"

plotHeatmap \
  --matrixFile "$out_dir/matrix_heatmap.gz" \
  --outFileName "$out_dir/heatmap_fc3_q8.png" \
  --samplesLabel "Vazyme_R1" "Vazyme_R2" "EpiCypher_R1" "EpiCypher_R2" \
  --regionsLabel "Shared" "Vazyme only" "EpiCypher only" \
  --colorList "blue,yellow,red" \
  --zMin 0 --zMax 2.5 \
  --missingDataColor white \
  --heatmapHeight 15 \
  --heatmapWidth 6 \
  --dpi 300 \
  --plotTitle "H3K27ac signal — FC>3 & -log10(q)>8" \
  --sortRegions descend \
  --sortUsing mean

echo "Saved: heatmap_fc3_q8.png"

# --------------------------------------------------------------------------
# [7] rGREAT ontology enrichment on the 3 peak subsets
#     For each of: shared, Vazyme-only, EpiCypher-only peaks (FC>3 & q>8)
#     - plotRegionGeneAssociations: distance to TSS, genomic context
#     - Top 10 bar plots for GO:BP, msigdb:H (Hallmarks), msigdb:C5
#     - Results saved as .rds for interactive shinyReport
#
#     Output files per subset (e.g. "shared"):
#       great_region_gene_shared.png
#       great_top10_GO_BP_shared.png
#       great_top10_msigdb_H_shared.png
#       great_top10_msigdb_C5_shared.png
#       great_res_shared.rds
#
#     To launch interactive Shiny report:
#       shinyReport(readRDS("great_res_shared.rds")[["GO:BP"]])
# --------------------------------------------------------------------------

Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(rGREAT)
  library(rtracklayer)
  library(ggplot2)
  library(dplyr)
})

out_dir <- "/Users/giacomocorleone/Desktop/Projects/vazyme_vs_epicypher/revision_round_2"

bed_files <- list(
  shared         = "peaks_shared.bed",
  vazyme_only    = "peaks_vazyme_only.bed",
  epicypher_only = "peaks_epicypher_only.bed"
)

bar_fill <- list(
  shared         = "#7570b3",
  vazyme_only    = "#1b9e77",
  epicypher_only = "royalblue3"
)

ontologies <- c("GO:BP", "msigdb:H", "msigdb:C5")

for (name in names(bed_files)) {
  cat("\n========== Processing:", name, "==========\n")
  gr <- suppressWarnings(import(file.path(out_dir, bed_files[[name]]), format = "BED"))

  res_list <- list()
  for (gs in ontologies) {
    cat("  GREAT:", gs, "\n")
    res_list[[gs]] <- suppressWarnings(
      great(gr, gene_sets = gs,
            tss_source = "TxDb.Hsapiens.UCSC.hg38.knownGene",
            cores = 4)
    )
  }

  # 1. Region-gene association plot
  p_assoc <- plotRegionGeneAssociations(res_list[["GO:BP"]])
  ggsave(file.path(out_dir, paste0("great_region_gene_", name, ".png")),
         plot = p_assoc, width = 10, height = 4, units = "in", dpi = 300, bg = "white")
  cat("  Saved: great_region_gene_", name, ".png\n", sep = "")

  # 2. Top 10 bar plots per ontology
  for (gs in ontologies) {
    tb    <- getEnrichmentTable(res_list[[gs]])
    label <- if ("description" %in% colnames(tb)) "description" else "id"
    top10 <- tb %>%
      arrange(p_adjust) %>%
      slice_head(n = 10) %>%
      mutate(term = .data[[label]],
             term = factor(term, levels = rev(term)))

    p <- ggplot(top10, aes(x = -log10(p_adjust + 1e-300), y = term)) +
      geom_bar(stat = "identity", fill = bar_fill[[name]], width = 0.7) +
      labs(
        title = paste0(gs, "  |  ", gsub("_", " ", name), " — top 10"),
        x     = expression(-log[10]~"(adj. p-value)"),
        y     = NULL
      ) +
      theme_classic(base_size = 12) +
      theme(
        axis.text.y  = element_text(color = "black", size = 10),
        axis.title.x = element_text(face = "bold"),
        plot.title   = element_text(face = "bold", size = 12, hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
      )

    gs_safe  <- gsub(":", "_", gs)
    out_file <- file.path(out_dir, paste0("great_top10_", gs_safe, "_", name, ".png"))
    ggsave(out_file, plot = p, width = 9, height = 5, units = "in", dpi = 300, bg = "white")
    cat("  Saved:", basename(out_file), "\n")
  }

  # 3. Save each GreatObject as a separate .rds (shinyReport requires a single
  #    GreatObject, not a named list — one file per ontology per subset)
  for (gs in ontologies) {
    gs_safe  <- gsub("[^a-zA-Z0-9]", "_", gs)
    out_file <- file.path(out_dir, paste0("great_res_", name, "_", gs_safe, ".rds"))
    saveRDS(res_list[[gs]], out_file)
    cat("  Saved:", basename(out_file), "\n")
  }
}

cat("\nAll done. Launch Shiny reports interactively in R:\n")
for (name in c("shared", "vazyme_only", "epicypher_only")) {
  for (gs in ontologies) {
    gs_safe <- gsub("[^a-zA-Z0-9]", "_", gs)
    cat("  shinyReport(readRDS('",
        file.path(out_dir, paste0("great_res_", name, "_", gs_safe, ".rds")),
        "'))\n", sep = "")
  }
}
EOF

# --------------------------------------------------------------------------
# [8] rGREAT — myeloma-specific gene sets (MSigDB C2)
#     59 myeloma gene sets extracted from MSigDB C2 curated collection using
#     msigdbr. These cover published transcriptional signatures from multiple
#     myeloma studies (Boylan, Chng, Corre, Davies, Munshi, Zhan, etc.).
#     Enrichment is tested for each of the 3 peak subsets to assess whether
#     the chromatin landscape reflects myeloma-relevant transcriptional programs.
#
#     COMPATIBILITY NOTE: msigdbr >= 10.0 renamed the Entrez ID column from
#     entrez_gene → ncbi_gene, and deprecated the category= argument in favour
#     of collection=. Both changes are applied here. rGREAT requires Entrez IDs
#     because its TSS extension uses TxDb (Entrez Gene ID type).
#
#     Outputs per subset:
#       great_top10_myeloma_C2_<subset>.png  — top 10 enriched myeloma sets
#       great_res_<subset>_myeloma_C2.rds    — GreatObject for shinyReport
#
#     shinyReport: shinyReport(readRDS("great_res_shared_myeloma_C2.rds"))
# --------------------------------------------------------------------------

Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(rGREAT)
  library(rtracklayer)
  library(msigdbr)
  library(ggplot2)
  library(dplyr)
})

out_dir <- "/Users/giacomocorleone/Desktop/Projects/vazyme_vs_epicypher/revision_round_2"

bed_files <- list(
  shared         = "peaks_shared.bed",
  vazyme_only    = "peaks_vazyme_only.bed",
  epicypher_only = "peaks_epicypher_only.bed"
)

bar_fill <- list(
  shared         = "#7570b3",
  vazyme_only    = "#1b9e77",
  epicypher_only = "royalblue3"
)

# Fetch C2 curated sets and filter for myeloma signatures.
# msigdbr >= 10.0: use collection= (not category=) and ncbi_gene (not entrez_gene)
c2_sets      <- msigdbr(species = "Homo sapiens", collection = "C2")
myeloma_df   <- c2_sets[grepl("MYELOMA", c2_sets$gs_name, ignore.case = TRUE), ]
myeloma_list <- split(as.character(myeloma_df$ncbi_gene), myeloma_df$gs_name)
myeloma_list <- lapply(myeloma_list, function(x) unique(x[!is.na(x) & x != "NA"]))
cat("Myeloma gene sets loaded:", length(myeloma_list), "\n\n")

for (name in names(bed_files)) {
  cat("Processing:", name, "\n")
  gr <- suppressWarnings(import(file.path(out_dir, bed_files[[name]]), format = "BED"))

  res <- suppressWarnings(
    great(gr, gene_sets = myeloma_list,
          tss_source = "TxDb.Hsapiens.UCSC.hg38.knownGene",
          cores = 4)
  )

  # Top 10 bar plot
  tb    <- getEnrichmentTable(res)
  label <- if ("description" %in% colnames(tb)) "description" else "id"
  top10 <- tb %>%
    arrange(p_adjust) %>%
    slice_head(n = 10) %>%
    mutate(term = gsub("_", " ", .data[[label]]),
           term = factor(term, levels = rev(term)))

  p <- ggplot(top10, aes(x = -log10(p_adjust + 1e-300), y = term)) +
    geom_bar(stat = "identity", fill = bar_fill[[name]], width = 0.7) +
    labs(
      title = paste0("MSigDB C2 Myeloma  |  ", gsub("_", " ", name), " — top 10"),
      x     = expression(-log[10]~"(adj. p-value)"),
      y     = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y  = element_text(color = "black", size = 9),
      axis.title.x = element_text(face = "bold"),
      plot.title   = element_text(face = "bold", size = 12, hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )

  out_plot <- file.path(out_dir, paste0("great_top10_myeloma_C2_", name, ".png"))
  ggsave(out_plot, plot = p, width = 10, height = 5, units = "in", dpi = 300, bg = "white")
  cat("  Saved:", basename(out_plot), "\n")

  # Save individual GreatObject for shinyReport
  out_rds <- file.path(out_dir, paste0("great_res_", name, "_myeloma_C2.rds"))
  saveRDS(res, out_rds)
  cat("  Saved:", basename(out_rds), "\n")
}

cat("\nLaunch Shiny reports:\n")
for (name in names(bed_files)) {
  cat("  shinyReport(readRDS('",
      file.path(out_dir, paste0("great_res_", name, "_myeloma_C2.rds")),
      "'))\n", sep = "")
}
EOF


# ═══════════════════════════════════════════════════════════════════════════════
# H3K27ac ChIP-seq — GREAT enrichment dot plot
# Multiple myeloma oncogenic gene sets · MSigDB C6
#
# Input:  /Users/giacomocorleone/Downloads/GREAT_H3K27ac_Myeloma_enrichment.csv  (same directory)
# Output: H3K27ac_myeloma_enrichment.pdf / .png
#
# Requires: tidyverse, scales
#   install.packages(c("tidyverse", "scales"))
# ═══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(scales)

# ── 1. Load data ─────────────────────────────────────────────────────────────
df <- read.csv("/Users/giacomocorleone/Downloads/GREAT_H3K27ac_Myeloma_enrichment.csv", stringsAsFactors = FALSE)

# ── 2. Short display labels for Y axis ───────────────────────────────────────
label_map <- c(
  "MOREAUX_MULTIPLE_MYELOMA_BY_TACI_DN"                    = "MM by TACI DN (Moreaux)",
  "BOYLAN_MULTIPLE_MYELOMA_PCA3_UP"                         = "MM PCA3 UP (Boylan)",
  "ZHAN_MULTIPLE_MYELOMA_MF_DN"                             = "MM MF subtype DN (Zhan)",
  "BOYLAN_MULTIPLE_MYELOMA_C_D_DN"                          = "MM C+D cluster DN (Boylan)",
  "SHAFFER_IRF4_TARGETS_IN_MYELOMA_VS_MATURE_B_LYMPHOCYTE"  = "IRF4 targets: MM vs B cells",
  "BOYLAN_MULTIPLE_MYELOMA_C_D_UP"                          = "MM C+D cluster UP (Boylan)",
  "CORRE_MULTIPLE_MYELOMA_DN"                               = "MM DN (Corre)",
  "ZHAN_MULTIPLE_MYELOMA_HP_DN"                             = "MM HP subtype DN (Zhan)",
  "ZHAN_MULTIPLE_MYELOMA_CD1_UP"                            = "MM CD1 UP (Zhan)",
  "BOYLAN_MULTIPLE_MYELOMA_C_DN"                            = "MM C subtype DN (Boylan)",
  "ZHAN_MULTIPLE_MYELOMA_CD1_VS_CD2_DN"                     = "MM CD1 vs CD2 DN (Zhan)",
  "ZHAN_MULTIPLE_MYELOMA_CD2_UP"                            = "MM CD2 UP (Zhan)",
  "ZHAN_MULTIPLE_MYELOMA_CD1_AND_CD2_UP"                    = "MM CD1+CD2 UP (Zhan)",
  "BOYLAN_MULTIPLE_MYELOMA_D_UP"                            = "MM subtype D UP (Boylan)",
  "ZHAN_MULTIPLE_MYELOMA_PR_DN"                             = "MM PR subtype DN (Zhan)",
  "BOYLAN_MULTIPLE_MYELOMA_C_CLUSTER_DN"                    = "MM C cluster DN (Boylan)",
  "BOYLAN_MULTIPLE_MYELOMA_D_DN"                            = "MM subtype D DN (Boylan)",
  "ZHAN_MULTIPLE_MYELOMA_LB_UP"                             = "MM LB subtype UP (Zhan)",
  "ZHAN_MULTIPLE_MYELOMA_UP"                                = "Multiple myeloma UP (Zhan)",
  "ZHAN_MULTIPLE_MYELOMA_CD1_DN"                            = "MM CD1 DN (Zhan)",
  "ZHAN_MULTIPLE_MYELOMA_CD1_AND_CD2_DN"                    = "MM CD1+CD2 DN (Zhan)",
  "SHAFFER_IRF4_MULTIPLE_MYELOMA_PROGRAM"                   = "IRF4 myeloma program"
)

df <- df %>%
  mutate(
    Gene_Set_Short = label_map[Gene_Set_Name],
    log10_p = -log10(Hyper_Adj_Pvalue),
    Peak_Set = factor(
      Peak_Set,
      levels = c("VAZ_only", "Shared_peaks", "EPI_only"),
      labels = c("VAZ only", "Shared peaks", "EPI only")
    )
  )

# ── 3. Order gene sets: most significant at top, least at bottom ──────────────
term_order <- df %>%
  group_by(Gene_Set_Short) %>%
  summarise(best_logp = max(log10_p), .groups = "drop") %>%
  arrange(best_logp) %>%          # ascending → most significant plots at top
  pull(Gene_Set_Short)

df$Gene_Set_Short <- factor(df$Gene_Set_Short, levels = term_order)

# ── 4. Plot ───────────────────────────────────────────────────────────────────
# Per-condition colours for X axis labels
x_label_colours <- c(
  "VAZ only"     = "#b04800",
  "Shared peaks" = "#444444",
  "EPI only"     = "#1050aa"
)

p <- ggplot(df, aes(
    x    = Peak_Set,
    y    = Gene_Set_Short,
    size = Binom_FoldEnrichment,
    fill = log10_p
  )) +

  # Alternating row background (mimics the figure stripes)
  geom_tile(
    aes(x = 2, width = Inf, fill = NULL),   # dummy tile spanning all columns
    data = df %>%
      distinct(Gene_Set_Short) %>%
      mutate(row = as.integer(Gene_Set_Short)) %>%
      filter(row %% 2 == 0),
    fill = "#f4f7fc", colour = NA, height = 1, inherit.aes = FALSE
  ) +

  geom_point(shape = 21, colour = "white", stroke = 0.4) +

  # Colour: RColorBrewer "Blues" (same palette used in the JS figure)
  scale_fill_distiller(
    palette   = "Blues",
    direction = 1,
    limits    = c(1.5, 8.5),
    oob       = scales::squish,
    name      = expression(-log[10](adj.~italic(p))),
    guide     = guide_colorbar(
      barwidth       = 0.7,
      barheight      = 7,
      title.position = "top",
      title.hjust    = 0.5,
      ticks.colour   = "grey40"
    )
  ) +

  # Size: binomial fold enrichment
  scale_size_continuous(
    name   = "Fold enrichment",
    range  = c(2, 10),
    breaks = c(1.3, 1.7, 2.1, 2.5),
    labels = c("1.3\u00d7", "1.7\u00d7", "2.1\u00d7", "2.5\u00d7"),
    guide  = guide_legend(
      override.aes   = list(fill = "#9ecae1", colour = "#3182bd", stroke = 0.6),
      title.position = "top",
      title.hjust    = 0.5
    )
  ) +

  scale_x_discrete(position = "top") +

  labs(
    title    = "H3K27ac ChIP-seq \u2014 multiple myeloma gene set enrichment",
    subtitle = "GREAT \u00b7 MSigDB Oncogenic Signatures \u00b7 hypergeometric adj. p < 0.05",
    x = NULL,
    y = NULL,
    caption = "Sorted by best hypergeometric adj. p \u00b7 top 10 terms per peak set \u00b7 MM = multiple myeloma"
  ) +

  theme_minimal(base_family = "Arial", base_size = 11) +
  theme(
    # Titles
    plot.title    = element_text(size = 12, face = "bold",
                                 hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 9, colour = "grey50",
                                 hjust = 0.5, margin = margin(b = 10)),
    plot.caption  = element_text(size = 7.5, colour = "grey65",
                                 hjust = 0, margin = margin(t = 8)),

    # Axes
    axis.text.x = element_text(
      size   = 11,
      face   = "bold",
      colour = unname(x_label_colours[levels(df$Peak_Set)])
    ),
    axis.text.y = element_text(size = 10, colour = "#111111"),

    # Grid
    panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.4,
                                       linetype = "dashed"),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(colour = "grey75", fill = NA,
                                      linewidth = 0.6),

    # Legend
    legend.position  = "right",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 9),
    legend.box       = "vertical",
    legend.spacing.y = unit(0.5, "cm"),

    # Background & margins
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin     = margin(12, 12, 12, 12)
  )

# ── 5. Export ─────────────────────────────────────────────────────────────────
# PDF (vector, best for journal submission) — requires cairo
ggsave("H3K27ac_myeloma_enrichment.pdf", p,
       width = 9, height = 8, units = "in",
       device = cairo_pdf)

# PNG (300 dpi raster, for presentations / preprints)
ggsave("H3K27ac_myeloma_enrichment.png", p,
       width = 9, height = 8, units = "in",
       dpi = 300, bg = "white")

message("Saved: H3K27ac_myeloma_enrichment.pdf  /  H3K27ac_myeloma_enrichment.png")

# --------------------------------------------------------------------------
# [9] Absolute distance to TSS grouped barplot (Python / matplotlib)
#     For each of the three peak subsets (shared, vazyme_only, epi_only),
#     region-gene associations from GREAT are binned by absolute distance
#     to the nearest TSS: 0 kb (overlap), 0–5 kb, 5–50 kb, 50–500 kb, >500 kb.
#     Bars show the fraction of associations in each bin; raw counts are
#     annotated above each bar.
#
#     Input:  revision_round_2/great_abs_distance.csv
#               columns: group, bin, bin_order, count, fraction
#     Output: revision_round_2/abs_distance_combined.png
#
#     Python dependencies: pandas, numpy, matplotlib
#       pip3 install pandas numpy matplotlib
# --------------------------------------------------------------------------

out_dir="/Users/giacomocorleone/Desktop/Projects/vazyme_vs_epicypher/revision_round_2"

python3 /Users/giacomocorleone/Desktop/Scripts/plot_abs_distance.py \
  "$out_dir/great_abs_distance.csv" \
  "$out_dir/abs_distance_combined.png"

# --------------------------------------------------------------------------
# [10] Combined MM217 vs MM196 myeloma enrichment dot plot
#      Merges the GREAT MSigDB C2 enrichment tables from both cohorts into
#      a single figure with MM217 and MM196 as side-by-side facets, sharing
#      the same y-axis (gene sets). Only gene sets with fold enrichment >= 1.8x
#      in at least one sample across both cohorts are shown.
#
#      This comparison directly reveals:
#        - Conserved myeloma signatures: enriched in both patients
#        - Patient-specific enrichments: present in only one cohort
#        - Method agreement: whether shared vs method-specific peaks drive
#          the same biological programmes across patients
#
#      Retained gene sets (FC >= 1.8x, union of both cohorts):
#        IRF4 targets: MM vs B cells  — both cohorts
#        IRF4 myeloma program         — both cohorts
#        Multiple myeloma UP (Zhan)   — both cohorts
#        MM CD1 DN (Zhan)             — MM217 Shared
#        MM hyperdiploid UP (Chng)    — MM196 Shared
#
#      Inputs:
#        revision_round_2/GREAT_H3K27ac_Myeloma_enrichment.csv      (MM217)
#        revision_round_2_mm196/GREAT_H3K27ac_MM196_C2_enrichment.csv (MM196)
#      Output:
#        revision_round_2_mm196/MM196_MM217_myeloma_enrichment_combined.png
#        revision_round_2_mm196/MM196_MM217_myeloma_enrichment_combined.pdf
# --------------------------------------------------------------------------

Rscript - <<'EOF'
suppressPackageStartupMessages({ library(tidyverse); library(scales) })

mm217_path <- "/Users/giacomocorleone/Desktop/Projects/vazyme_vs_epicypher/revision_round_2/GREAT_H3K27ac_Myeloma_enrichment.csv"
mm196_path <- "/Users/giacomocorleone/Desktop/Projects/vazyme_vs_epicypher/revision_round_2_mm196/GREAT_H3K27ac_MM196_C2_enrichment.csv"
out_dir    <- "/Users/giacomocorleone/Desktop/Projects/vazyme_vs_epicypher/revision_round_2_mm196"

label_map <- c(
  "BOYLAN_MULTIPLE_MYELOMA_C_CLUSTER_DN"                   = "MM C cluster DN (Boylan)",
  "BOYLAN_MULTIPLE_MYELOMA_C_D_DN"                         = "MM C+D cluster DN (Boylan)",
  "BOYLAN_MULTIPLE_MYELOMA_C_D_UP"                         = "MM C+D cluster UP (Boylan)",
  "BOYLAN_MULTIPLE_MYELOMA_C_DN"                           = "MM C subtype DN (Boylan)",
  "BOYLAN_MULTIPLE_MYELOMA_D_DN"                           = "MM subtype D DN (Boylan)",
  "BOYLAN_MULTIPLE_MYELOMA_D_UP"                           = "MM subtype D UP (Boylan)",
  "BOYLAN_MULTIPLE_MYELOMA_PCA3_UP"                        = "MM PCA3 UP (Boylan)",
  "CHNG_MULTIPLE_MYELOMA_HYPERPLOID_UP"                    = "MM hyperdiploid UP (Chng)",
  "CORRE_MULTIPLE_MYELOMA_DN"                              = "MM DN (Corre)",
  "MAGRANGEAS_MULTIPLE_MYELOMA_IGLL_VS_IGLK_UP"            = "MM IgLL vs IgLK UP (Magrangeas)",
  "MOREAUX_MULTIPLE_MYELOMA_BY_TACI_DN"                    = "MM by TACI DN (Moreaux)",
  "MOREAUX_MULTIPLE_MYELOMA_BY_TACI_UP"                    = "MM by TACI UP (Moreaux)",
  "SHAFFER_IRF4_MULTIPLE_MYELOMA_PROGRAM"                  = "IRF4 myeloma program",
  "SHAFFER_IRF4_TARGETS_IN_MYELOMA_VS_MATURE_B_LYMPHOCYTE" = "IRF4 targets: MM vs B cells",
  "ZHAN_MULTIPLE_MYELOMA_CD1_AND_CD2_DN"                   = "MM CD1+CD2 DN (Zhan)",
  "ZHAN_MULTIPLE_MYELOMA_CD1_AND_CD2_UP"                   = "MM CD1+CD2 UP (Zhan)",
  "ZHAN_MULTIPLE_MYELOMA_CD1_DN"                           = "MM CD1 DN (Zhan)",
  "ZHAN_MULTIPLE_MYELOMA_CD1_UP"                           = "MM CD1 UP (Zhan)",
  "ZHAN_MULTIPLE_MYELOMA_CD1_VS_CD2_DN"                    = "MM CD1 vs CD2 DN (Zhan)",
  "ZHAN_MULTIPLE_MYELOMA_CD2_UP"                           = "MM CD2 UP (Zhan)",
  "ZHAN_MULTIPLE_MYELOMA_HP_DN"                            = "MM HP subtype DN (Zhan)",
  "ZHAN_MULTIPLE_MYELOMA_LB_DN"                            = "MM LB subtype DN (Zhan)",
  "ZHAN_MULTIPLE_MYELOMA_LB_UP"                            = "MM LB subtype UP (Zhan)",
  "ZHAN_MULTIPLE_MYELOMA_MF_DN"                            = "MM MF subtype DN (Zhan)",
  "ZHAN_MULTIPLE_MYELOMA_MF_UP"                            = "MM MF subtype UP (Zhan)",
  "ZHAN_MULTIPLE_MYELOMA_PR_DN"                            = "MM PR subtype DN (Zhan)",
  "ZHAN_MULTIPLE_MYELOMA_UP"                               = "Multiple myeloma UP (Zhan)"
)

df_all <- bind_rows(
  read.csv(mm217_path, stringsAsFactors=FALSE) %>% mutate(Cohort="MM217"),
  read.csv(mm196_path, stringsAsFactors=FALSE) %>% mutate(Cohort="MM196")
) %>%
  mutate(
    Gene_Set_Short = coalesce(label_map[Gene_Set_Name], Gene_Set_Name),
    log10_p  = -log10(Hyper_Adj_Pvalue),
    Peak_Set = factor(Peak_Set,
                      levels = c("VAZ_only","Shared_peaks","EPI_only"),
                      labels = c("Vazyme only","Shared peaks","EpiCypher only")),
    Cohort   = factor(Cohort, levels=c("MM217","MM196"))
  )

# Keep gene sets with FC >= 1.8 in at least one row across both cohorts
keep_sets <- df_all %>%
  group_by(Gene_Set_Name) %>%
  summarise(max_fc = max(Binom_FoldEnrichment), .groups="drop") %>%
  filter(max_fc >= 1.8) %>%
  pull(Gene_Set_Name)

df_filt <- df_all %>% filter(Gene_Set_Name %in% keep_sets)
cat("Gene sets retained (FC >= 1.8x):", length(unique(df_filt$Gene_Set_Name)), "
")

term_order <- df_filt %>%
  group_by(Gene_Set_Short) %>%
  summarise(best_logp = max(log10_p), .groups="drop") %>%
  arrange(best_logp) %>%
  pull(Gene_Set_Short)

df_filt$Gene_Set_Short <- factor(df_filt$Gene_Set_Short, levels=term_order)

x_label_colours <- c(
  "Vazyme only"    = "#b04800",
  "Shared peaks"   = "#444444",
  "EpiCypher only" = "#08306b"
)

p <- ggplot(df_filt, aes(x=Peak_Set, y=Gene_Set_Short,
                          size=Binom_FoldEnrichment, fill=log10_p)) +
  geom_point(shape=21, colour="white", stroke=0.4) +
  facet_grid(. ~ Cohort, switch="x") +
  scale_fill_gradient(
    low="#6baed6", high="#08306b",
    name  = expression(-log[10](adj.~italic(p))),
    guide = guide_colorbar(barwidth=0.7, barheight=7,
                           title.position="top", title.hjust=0.5,
                           ticks.colour="grey40")
  ) +
  scale_size_continuous(
    name   = "Fold enrichment",
    range  = c(1.5, 9),
    breaks = c(1.5, 2.0, 2.5),
    labels = c("1.5x","2.0x","2.5x"),
    guide  = guide_legend(
      override.aes   = list(fill="#2171b5", colour="white", stroke=0.4),
      title.position = "top", title.hjust=0.5
    )
  ) +
  scale_x_discrete(position="top") +
  labs(
    title    = "H3K27ac ChIP-seq — multiple myeloma gene set enrichment",
    subtitle = "GREAT · MSigDB C2 · MM217 vs MM196 · fold enrichment >= 1.8x in >= 1 sample · hyper adj. p < 0.05",
    x=NULL, y=NULL,
    caption  = "Sorted by best hypergeometric adj. p across both cohorts · MM = multiple myeloma"
  ) +
  theme_minimal(base_size=11) +
  theme(
    plot.title         = element_text(size=12, face="bold", hjust=0.5, margin=margin(b=4)),
    plot.subtitle      = element_text(size=9, colour="grey50", hjust=0.5, margin=margin(b=10)),
    plot.caption       = element_text(size=7.5, colour="grey65", hjust=0, margin=margin(t=8)),
    axis.text.x        = element_text(size=10, face="bold",
                                       colour=unname(x_label_colours[levels(df_filt$Peak_Set)])),
    axis.text.y        = element_text(size=9.5, colour="#111111"),
    strip.text         = element_text(size=12, face="bold"),
    strip.placement    = "outside",
    panel.spacing      = unit(1.2, "lines"),
    panel.grid.major.x = element_line(colour="grey85", linewidth=0.4, linetype="dashed"),
    panel.grid.major.y = element_line(colour="grey90", linewidth=0.4),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(colour="grey75", fill=NA, linewidth=0.6),
    legend.position    = "right",
    legend.title       = element_text(size=9, face="bold"),
    legend.text        = element_text(size=9),
    legend.box         = "vertical",
    legend.spacing.y   = unit(0.5, "cm"),
    plot.background    = element_rect(fill="white", colour=NA),
    plot.margin        = margin(12, 12, 12, 12)
  )

ggsave(file.path(out_dir, "MM196_MM217_myeloma_enrichment_combined.pdf"), p,
       width=14, height=7, units="in", device=cairo_pdf)
ggsave(file.path(out_dir, "MM196_MM217_myeloma_enrichment_combined.png"), p,
       width=14, height=7, units="in", dpi=300, bg="white")
message("Saved: MM196_MM217_myeloma_enrichment_combined.pdf / .png")
EOF
