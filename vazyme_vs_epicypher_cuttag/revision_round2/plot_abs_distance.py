#!/usr/bin/env python3
"""
Plot GREAT region-gene associations binned by absolute distance to TSS
for three sets (shared, vazyme_only, epi_only) on one grouped barplot.

Usage:
    python plot_abs_distance.py great_abs_distance.csv out.png
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ---- inputs ----
csv_path = sys.argv[1] if len(sys.argv) > 1 else "great_abs_distance.csv"
out_path = sys.argv[2] if len(sys.argv) > 2 else "abs_distance_combined.png"

df = pd.read_csv(csv_path)

# Fixed ordering for bins and groups
bin_levels   = ["0", "(0, 5]", "(5, 50]", "(50, 500]", "> 500"]
group_levels = ["shared", "vazyme_only", "epi_only"]
colors       = {"shared": "#1f4e8c", "vazyme_only": "#3a86ff", "epi_only": "#8ecae6"}

# Pivot to a bin x group matrix of fractions
mat = (df.pivot(index="bin", columns="group", values="fraction")
         .reindex(index=bin_levels, columns=group_levels))

# Same for raw counts (used for the labels above each bar)
cnt = (df.pivot(index="bin", columns="group", values="count")
         .reindex(index=bin_levels, columns=group_levels))

# ---- plot ----
x     = np.arange(len(bin_levels))
width = 0.26

fig, ax = plt.subplots(figsize=(9, 5.5))

for i, g in enumerate(group_levels):
    offset = (i - 1) * width
    bars = ax.bar(x + offset, mat[g].values, width,
                  label=g, color=colors[g], edgecolor="black", linewidth=0.4)
    # annotate raw counts above bars
    for rect, c in zip(bars, cnt[g].values):
        ax.text(rect.get_x() + rect.get_width() / 2,
                rect.get_height() + 0.008,
                f"{int(c)}", ha="center", va="bottom", fontsize=7, rotation=90)

ax.set_xticks(x)
ax.set_xticklabels(bin_levels)
ax.set_xlabel("Absolute distance to TSS (kb)")
ax.set_ylabel("Region-gene associations (fraction)")
ax.set_title("Binned by absolute distance to TSS")
ax.set_ylim(0, max(mat.max()) * 1.18)
ax.legend(title="Set", frameon=False)
ax.spines[["top", "right"]].set_visible(False)

fig.tight_layout()
fig.savefig(out_path, dpi=200)
print(f"Saved: {out_path}")
