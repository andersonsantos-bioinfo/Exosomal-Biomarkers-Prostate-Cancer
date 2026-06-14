# Load DEG results
res <- read.table(
  "results/DESeq2_FullResults_Tumor_vs_Normal.txt",
  header = TRUE,
  sep = "\t"
)

# Up-regulated genes
up_genes <- subset(
  res,
  log2FoldChange > 1 &
    padj < 0.05
)

# Down-regulated genes
down_genes <- subset(
  res,
  log2FoldChange < -1 &
    padj < 0.05
)

# Save outputs
write.table(
  up_genes,
  "results/Upregulated_Genes.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  down_genes,
  "results/Downregulated_Genes.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
