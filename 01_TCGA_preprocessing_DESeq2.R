# =========================================================
# TCGA-PRAD preprocessing and differential expression
# Exosomal Biomarkers in Prostate Cancer
# =========================================================

# ---------------------------------------------------------
# Load libraries
# ---------------------------------------------------------

library(SummarizedExperiment)
library(DESeq2)

# ---------------------------------------------------------
# Define project directory
# ---------------------------------------------------------

dir_path <- "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer"

# ---------------------------------------------------------
# Create output directories
# ---------------------------------------------------------

dir.create(
  file.path(dir_path, "results"),
  showWarnings = FALSE
)

dir.create(
  file.path(dir_path, "scripts"),
  showWarnings = FALSE
)

# ---------------------------------------------------------
# Load SummarizedExperiment object
# ---------------------------------------------------------

se <- readRDS(
  file.path(
    dir_path,
    "data_processed",
    "TCGA_PRAD_STARcounts_se.rds"
  )
)

# ---------------------------------------------------------
# Inspect object
# ---------------------------------------------------------

print(se)

# ---------------------------------------------------------
# Extract raw count matrix
# ---------------------------------------------------------

counts_matrix <- assay(se)

# Matrix dimensions
dim(counts_matrix)

# Expected:
# 60660 genes x 553 samples

# ---------------------------------------------------------
# Extract metadata
# ---------------------------------------------------------

metadata <- as.data.frame(colData(se))

# Add sample IDs
metadata$Sample <- rownames(metadata)

# ---------------------------------------------------------
# Define sample groups
# ---------------------------------------------------------

sample_type <- metadata$shortLetterCode

group <- ifelse(
  sample_type == "NT",
  "Normal",
  "Tumor"
)

group <- factor(
  group,
  levels = c("Normal", "Tumor")
)

# Add group to metadata
metadata$Group <- group

# Check distribution
table(metadata$Group)

# Expected:
# Normal = 52
# Tumor = 501

# ---------------------------------------------------------
# Filter low-expression genes
# ---------------------------------------------------------

keep <- rowSums(counts_matrix >= 10) >= 5

counts_final <- counts_matrix[keep, ]

# Dimensions after filtering
dim(counts_final)

# ---------------------------------------------------------
# Build DESeq2 dataset
# ---------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_final),
  colData = metadata,
  design = ~ Group
)

# ---------------------------------------------------------
# Run differential expression analysis
# ---------------------------------------------------------

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("Group", "Tumor", "Normal")
)

# ---------------------------------------------------------
# Order by adjusted p-value
# ---------------------------------------------------------

res <- res[order(res$padj), ]

# ---------------------------------------------------------
# Convert to data frame
# ---------------------------------------------------------

res_df <- as.data.frame(res)

# Add Ensembl IDs
res_df$GeneID <- rownames(res_df)

# ---------------------------------------------------------
# Remove NA adjusted p-values
# ---------------------------------------------------------

res_df <- subset(
  res_df,
  !is.na(padj)
)

# ---------------------------------------------------------
# Define significant DEGs
# ---------------------------------------------------------

deg_standard <- subset(
  res_df,
  abs(log2FoldChange) >= 1 &
    padj < 0.05
)

deg_strict <- subset(
  res_df,
  abs(log2FoldChange) >= 2 &
    padj < 0.01
)

# ---------------------------------------------------------
# Save full DESeq2 results
# ---------------------------------------------------------

write.table(
  res_df,
  file = file.path(
    dir_path,
    "results",
    "DESeq2_FullResults_Tumor_vs_Normal.txt"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ---------------------------------------------------------
# Save standard DEG list
# ---------------------------------------------------------

write.table(
  deg_standard,
  file = file.path(
    dir_path,
    "results",
    "DEGs_Standard_log2FC1_padj005.txt"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ---------------------------------------------------------
# Save strict DEG list
# ---------------------------------------------------------

write.table(
  deg_strict,
  file = file.path(
    dir_path,
    "results",
    "DEGs_Strict_log2FC2_padj001.txt"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ---------------------------------------------------------
# Save filtered count matrix
# ---------------------------------------------------------

write.table(
  counts_final,
  file = file.path(
    dir_path,
    "results",
    "PRAD_TCGA_Filtered_Counts.txt"
  ),
  sep = "\t",
  quote = FALSE
)

# ---------------------------------------------------------
# Session information
# ---------------------------------------------------------

sessionInfo()