# =========================================================
# 04_enrichment_analysis.R
# Functional Enrichment Analysis - TCGA PRAD
# Exosomal Biomarkers in Prostate Cancer
# =========================================================

# =========================================================
# Load libraries
# =========================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# =========================================================
# Define project directory
# =========================================================

dir_path <- "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer"

# =========================================================
# Create enrichment directory
# =========================================================

dir.create(
  file.path(dir_path, "enrichment"),
  showWarnings = FALSE
)

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

# =========================================================
# Remove NA adjusted p-values
# =========================================================

res <- subset(
  res,
  !is.na(padj)
)

# =========================================================
# Select upregulated genes
# =========================================================

up_genes <- subset(
  res,
  log2FoldChange >= 1 &
    padj < 0.05
)

# =========================================================
# Extract Gene IDs
# =========================================================

gene_list <- up_genes$GeneID

# Remove Ensembl version numbers
gene_list <- sub("\\..*", "", gene_list)

# =========================================================
# Convert Ensembl IDs to ENTREZ IDs
# =========================================================

gene_conversion <- bitr(
  gene_list,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Extract ENTREZ IDs
entrez_genes <- unique(gene_conversion$ENTREZID)

# =========================================================
# GO Biological Process enrichment
# =========================================================

ego_bp <- enrichGO(
  gene = entrez_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# =========================================================
# KEGG enrichment
# =========================================================

ekegg <- enrichKEGG(
  gene = entrez_genes,
  organism = "hsa",
  pvalueCutoff = 0.05
)

# =========================================================
# Reactome enrichment
# =========================================================

library(ReactomePA)

ereact <- enrichPathway(
  gene = entrez_genes,
  organism = "human",
  pvalueCutoff = 0.05,
  readable = TRUE
)

# =========================================================
# Save enrichment tables
# =========================================================

write.table(
  as.data.frame(ego_bp),
  file = file.path(
    dir_path,
    "enrichment",
    "GO_Enrichment_Up_genes.txt"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  as.data.frame(ekegg),
  file = file.path(
    dir_path,
    "enrichment",
    "KEGG_Enrichment_Up_genes.txt"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  as.data.frame(ereact),
  file = file.path(
    dir_path,
    "enrichment",
    "Reactome_Enrichment_Up_genes.txt"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# =========================================================
# Create GO dotplot
# =========================================================

go_plot <- dotplot(
  ego_bp,
  showCategory = 15
) +
  ggtitle("GO Biological Process Enrichment")

# Show plot
print(go_plot)

# =========================================================
# Save GO enrichment plot
# =========================================================

ggsave(
  filename = file.path(
    dir_path,
    "enrichment",
    "GO_Enrichment_Dotplot.pdf"
  ),
  plot = go_plot,
  width = 10,
  height = 8
)

# =========================================================
# Create KEGG dotplot
# =========================================================

kegg_plot <- dotplot(
  ekegg,
  showCategory = 15
) +
  ggtitle("KEGG Pathway Enrichment")

# Show plot
print(kegg_plot)

# =========================================================
# Save KEGG plot
# =========================================================

ggsave(
  filename = file.path(
    dir_path,
    "enrichment",
    "KEGG_Enrichment_Dotplot.pdf"
  ),
  plot = kegg_plot,
  width = 10,
  height = 8
)

# =========================================================
# Create Reactome dotplot
# =========================================================

reactome_plot <- dotplot(
  ereact,
  showCategory = 15
) +
  ggtitle("Reactome Pathway Enrichment")

# Show plot
print(reactome_plot)

# =========================================================
# Save Reactome plot
# =========================================================

ggsave(
  filename = file.path(
    dir_path,
    "enrichment",
    "Reactome_Enrichment_Dotplot.pdf"
  ),
  plot = reactome_plot,
  width = 10,
  height = 8
)

# =========================================================
# Session information
# =========================================================

sessionInfo()