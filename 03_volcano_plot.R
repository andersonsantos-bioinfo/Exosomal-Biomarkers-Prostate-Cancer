# =========================================================
# 03_volcano_plot.R
# Volcano Plot - TCGA PRAD
# Exosomal Biomarkers in Prostate Cancer
# =========================================================

# ---------------------------------------------------------
# Define project directory
# ---------------------------------------------------------

dir_path <- "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer"

# ---------------------------------------------------------
# Create figures directory
# ---------------------------------------------------------

dir.create(
  file.path(dir_path, "figures"),
  showWarnings = FALSE
)

# =========================================================
# Load libraries
# =========================================================

library(ggplot2)
library(dplyr)

# =========================================================
# Load DESeq2 results
# =========================================================

res <- read.table(
  file.path(
    dir_path,
    "results",
    "DESeq2_FullResults_Tumor_vs_Normal.txt"
  ),
  header = TRUE,
  sep = "\t"
)

# Inspect dimensions
dim(res)

# View first rows
head(res)

# =========================================================
# Remove NA adjusted p-values
# =========================================================

res <- subset(
  res,
  !is.na(padj)
)

# Check dimensions after filtering
dim(res)

# =========================================================
# Create significance categories
# =========================================================

res$Significance <- "Not Significant"

# Upregulated genes
res$Significance[
  res$log2FoldChange >= 1 &
    res$padj < 0.05
] <- "Upregulated"

# Downregulated genes
res$Significance[
  res$log2FoldChange <= -1 &
    res$padj < 0.05
] <- "Downregulated"

# Check category distribution
table(res$Significance)

# =========================================================
# Create volcano plot
# =========================================================

volcano_plot <- ggplot(
  res,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = Significance
  )
) +
  
  geom_point(
    alpha = 0.7,
    size = 1.5
  ) +
  
  scale_color_manual(
    values = c(
      "Upregulated" = "#B22222",
      "Downregulated" = "#1F4E79",
      "Not Significant" = "gray70"
    )
  ) +
  
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "black"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "black"
  ) +
  
  labs(
    title = "Differential Gene Expression in TCGA-PRAD",
    subtitle = "Tumor versus Normal prostate tissue",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    
    legend.title = element_blank()
  )

# =========================================================
# Show volcano plot
# =========================================================

print(volcano_plot)

# =========================================================
# Save PDF
# =========================================================

ggsave(
  filename = file.path(
    dir_path,
    "figures",
    "VolcanoPlot_TCGA_PRAD.pdf"
  ),
  plot = volcano_plot,
  width = 8,
  height = 6
)

# =========================================================
# Save PNG
# =========================================================

ggsave(
  filename = file.path(
    dir_path,
    "figures",
    "VolcanoPlot_TCGA_PRAD.png"
  ),
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# =========================================================
# Session information
# =========================================================

sessionInfo()