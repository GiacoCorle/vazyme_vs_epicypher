# ============================================
# STEP 08: PEAK COUNT VS READ DEPTH PLOT
# PURPOSE:
#   Plot the number of detected peaks as a function of
#   sequencing depth for Vazyme and EpiCypher samples.
# INPUT:
#   Peak counts summarized manually from filtered peak files.
# OUTPUT:
#   Publication-ready PNG figure.
# REQUIREMENTS:
#   R packages: tidyr, ggplot2, dplyr.
# ============================================

df <- data.frame(
  percentage = seq(10, 100, by = 10),
  vazyme_R1 = c(10375, 21581, 31925, 38371, 47251, 50998, 59007, 61192, 68695, 69842),
  vazyme_R2 = c(8292, 21598, 30788, 41978, 47463, 57531, 60483, 69925, 71430, 79934),
  epycypher_R1 = c(5225, 18969, 33706, 53417, 66468, 87476, 97535, 118710, 125736, 146361),
  epycypher_R2 = c(3346, 12983, 22633, 32694, 41098, 48275, 57693, 63793, 69264, 78113)
)

df$percentage <- seq(10, 100, by = 10)

library(tidyr)
library(ggplot2)
library(dplyr)

# Reshape data into long format for plotting
peak_counts_long <- pivot_longer(
  df,
  cols = c(vazyme_R1, vazyme_R2, epycypher_R1, epycypher_R2),
  names_to = "sample",
  values_to = "peaks"
)

# Annotate method and replicate based on sample name
peak_counts_long <- peak_counts_long %>%
  mutate(
    method = ifelse(grepl("vazyme", sample), "Vazyme", "EpiCypher"),
    replicate = ifelse(grepl("R1", sample), "R1", "R2")
  )

p <- ggplot(peak_counts_long, aes(x = percentage, y = log10(peaks), color = method, linetype = replicate)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Vazyme" = "#1b9e77", "EpiCypher" = "royalblue3")) +
  labs(
    x = "Percentage of Reads (%)",
    y = expression(Log[10]~"(Number of Peaks)"),
    color = "Method",
    linetype = "Replicate"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(10, 10, 10, 10)
  )

# Save figure
setwd("/Users/giacomocorleone/Desktop/Projects/vazyme_vs_epicypher/paper/")
ggsave(
  "peaks_rev1_rep1.png",
  plot = p,
  width = 6,
  height = 4,
  units = "in",
  dpi = 300,
  bg = "white"
)
