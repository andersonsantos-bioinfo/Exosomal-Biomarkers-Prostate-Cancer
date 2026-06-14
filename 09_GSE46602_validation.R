# =========================================================
# 09_GSE46602_annotation_normalization.R
# GSE46602 - Annotation and Normalization
# External Validation Cohort
# =========================================================

# =========================================================
# Load required packages
# =========================================================

library(GEOquery)
library(limma)
library(dplyr)
library(stringr)

# =========================================================
# Define project directory
# =========================================================

dir_path <- "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer"

# =========================================================
# Create folders
# =========================================================

dir.create(
  file.path(dir_path, "validation"),
  showWarnings = FALSE
)

dir.create(
  file.path(dir_path, "validation", "GSE46602"),
  showWarnings = FALSE
)

# =========================================================
# Load normalized expression matrix
# =========================================================

expr <- read.table(
  file.path(
    dir_path,
    "validation",
    "GSE46602",
    "GSE46602_expression_normalized.txt"
  ),
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

# =========================================================
# Load annotation table (robust GEO import)
# =========================================================

annot <- read.table(
  file.path(
    dir_path,
    "validation",
    "GSE46602",
    "GSE46602_annotation.txt"
  ),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fill = TRUE,
  quote = "",
  comment.char = ""
)

# =========================================================
# Inspect annotation columns
# =========================================================

print(colnames(annot))

# =========================================================
# Detect gene symbol column automatically
# =========================================================

possible_cols <- c(
  "Gene Symbol",
  "Gene symbol",
  "GENE_SYMBOL",
  "Symbol",
  "gene_assignment",
  "GeneSymbol"
)

coluna_gene <- possible_cols[
  possible_cols %in% colnames(annot)
][1]

# =========================================================
# Stop if no annotation column is found
# =========================================================

if(is.na(coluna_gene)){
  
  stop(
    paste(
      "Gene symbol column not found.\n",
      "Check annotation columns manually using:\n",
      "colnames(annot)"
    )
  )
  
}

cat("Gene symbol column detected:\n")
print(coluna_gene)

# =========================================================
# Create annotation dataframe
# =========================================================

annot_df <- data.frame(
  Probe = rownames(annot),
  GeneSymbol = annot[[coluna_gene]],
  stringsAsFactors = FALSE
)

# =========================================================
# Remove empty annotations
# =========================================================

annot_df <- annot_df %>%
  filter(
    !is.na(GeneSymbol),
    GeneSymbol != "",
    GeneSymbol != "---"
  )

# =========================================================
# Remove multiple mappings
# Example: "GENE1 /// GENE2"
# =========================================================

annot_df$GeneSymbol <- str_split_fixed(
  annot_df$GeneSymbol,
  " /// ",
  2
)[,1]

# =========================================================
# Remove spaces
# =========================================================

annot_df$GeneSymbol <- trimws(
  annot_df$GeneSymbol
)

# =========================================================
# Match probes with expression matrix
# =========================================================

common_probes <- intersect(
  rownames(expr),
  annot_df$Probe
)

cat("Common probes:\n")
print(length(common_probes))

expr <- expr[common_probes, ]

annot_df <- annot_df %>%
  filter(
    Probe %in% common_probes
  )

# =========================================================
# Convert expression matrix to dataframe
# =========================================================

expr_df <- data.frame(
  Probe = rownames(expr),
  expr,
  stringsAsFactors = FALSE
)

# =========================================================
# Merge annotation + expression
# =========================================================

expr_annot <- merge(
  annot_df,
  expr_df,
  by = "Probe"
)

# =========================================================
# Calculate average expression
# =========================================================

expr_annot$MeanExpr <- rowMeans(
  expr_annot[, -(1:2)],
  na.rm = TRUE
)

# =========================================================
# Keep highest expressed probe per gene
# =========================================================

expr_annot <- expr_annot %>%
  arrange(
    desc(MeanExpr)
  )

expr_final <- expr_annot %>%
  distinct(
    GeneSymbol,
    .keep_all = TRUE
  )

# =========================================================
# Remove helper columns
# =========================================================

expr_final <- expr_final %>%
  dplyr::select(
    -Probe,
    -MeanExpr
  )

# =========================================================
# Set gene symbols as rownames
# =========================================================

rownames(expr_final) <- expr_final$GeneSymbol

expr_final <- expr_final[, -1]

# =========================================================
# Final dimensions
# =========================================================

cat("Final annotated matrix dimensions:\n")
print(dim(expr_final))

# =========================================================
# Define multigene panel
# =========================================================

genes_painel <- c(
  "ETV1",
  "HOXC6",
  "SLC45A2",
  "ZIC2",
  "SCHLAP1",
  "LINC01475"
)

# =========================================================
# Check genes present in platform
# =========================================================

genes_existentes <- genes_painel[
  genes_painel %in% rownames(expr_final)
]

cat("Genes identified in GSE46602:\n")
print(genes_existentes)

# =========================================================
# Check missing genes
# =========================================================

genes_ausentes <- setdiff(
  genes_painel,
  genes_existentes
)

cat("Genes absent from platform:\n")
print(genes_ausentes)

# =========================================================
# Extract multigene panel matrix
# =========================================================

expr_panel <- expr_final[
  genes_existentes,
]

# =========================================================
# Transpose matrix
# =========================================================

expr_panel_t <- as.data.frame(
  t(expr_panel)
)

# =========================================================
# Add sample IDs
# =========================================================

expr_panel_t$Sample <- rownames(
  expr_panel_t
)

# =========================================================
# Reorder columns
# =========================================================

expr_panel_t <- expr_panel_t %>%
  dplyr::select(
    Sample,
    everything()
  )

# =========================================================
# Save annotated matrix
# =========================================================

write.table(
  expr_final,
  file = file.path(
    dir_path,
    "validation",
    "GSE46602",
    "GSE46602_Annotated_Matrix.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

# =========================================================
# Save panel matrix
# =========================================================

write.table(
  expr_panel_t,
  file = file.path(
    dir_path,
    "validation",
    "GSE46602",
    "GSE46602_Panel_Expression.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================================================
# Save identified genes
# =========================================================

write.table(
  data.frame(
    Identified_Genes = genes_existentes
  ),
  file = file.path(
    dir_path,
    "validation",
    "GSE46602",
    "Genes_Identified_GSE46602.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================================================
# Save missing genes
# =========================================================

write.table(
  data.frame(
    Missing_Genes = genes_ausentes
  ),
  file = file.path(
    dir_path,
    "validation",
    "GSE46602",
    "Genes_Missing_GSE46602.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================================================
# Save R objects
# =========================================================

saveRDS(
  expr_final,
  file = file.path(
    dir_path,
    "validation",
    "GSE46602",
    "GSE46602_Annotated_Matrix.rds"
  )
)

saveRDS(
  expr_panel_t,
  file = file.path(
    dir_path,
    "validation",
    "GSE46602",
    "GSE46602_Panel_Expression.rds"
  )
)

# =========================================================
# Save session info
# =========================================================

writeLines(
  capture.output(sessionInfo()),
  file.path(
    dir_path,
    "session_info",
    "sessionInfo_GSE46602_annotation.txt"
  )
)

# =========================================================
# Final message
# =========================================================

cat("\n")
cat("=========================================\n")
cat("GSE46602 annotation completed successfully\n")
cat("=========================================\n")