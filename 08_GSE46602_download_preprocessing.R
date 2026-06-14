# =========================================================
# 08_GSE46602_download_preprocessing.R
# External Validation Cohort - GSE46602
# Download and Preprocessing
# =========================================================

# =========================================================
# Load required packages
# =========================================================

library(GEOquery)
library(limma)
library(dplyr)

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
# Set working directory
# =========================================================

setwd(
  file.path(
    dir_path,
    "validation",
    "GSE46602"
  )
)

# =========================================================
# Download GEO dataset
# =========================================================

gse <- getGEO(
  "GSE46602",
  GSEMatrix = TRUE,
  AnnotGPL = TRUE
)

# =========================================================
# Select expression object
# =========================================================

if(length(gse) > 1){
  idx <- grep("GPL", attr(gse, "names"))
  expr_set <- gse[[idx[1]]]
} else {
  expr_set <- gse[[1]]
}

# =========================================================
# Extract expression matrix
# =========================================================

expr <- exprs(expr_set)

dim(expr)

head(expr[,1:5])

# =========================================================
# Extract phenotype metadata
# =========================================================

pheno <- pData(expr_set)

dim(pheno)

head(pheno)

# =========================================================
# Extract feature annotation
# =========================================================

annot <- fData(expr_set)

dim(annot)

head(annot)

# =========================================================
# Check expression distribution
# =========================================================

summary(expr)

qx <- as.numeric(
  quantile(
    expr,
    c(0, 0.25, 0.5, 0.75, 0.99, 1.0),
    na.rm = TRUE
  )
)

qx

# =========================================================
# Log2 transformation if necessary
# =========================================================

LogC <- (
  qx[5] > 100 ||
    (qx[6] - qx[1] > 50 && qx[2] > 0)
)

if(LogC){
  
  expr[expr <= 0] <- NA
  
  expr <- log2(expr)
  
  cat("Log2 transformation applied\n")
  
} else {
  
  cat("Expression matrix already log2 transformed\n")
}

# =========================================================
# Remove probes with many missing values
# =========================================================

keep <- rowSums(is.na(expr)) < ncol(expr) * 0.5

expr <- expr[keep, ]

dim(expr)

# =========================================================
# Replace remaining missing values
# =========================================================

expr <- avereps(expr)

expr[is.na(expr)] <- median(
  expr,
  na.rm = TRUE
)

# =========================================================
# Quantile normalization
# =========================================================

expr_norm <- normalizeBetweenArrays(
  expr,
  method = "quantile"
)

# =========================================================
# Verify normalized matrix
# =========================================================

summary(expr_norm)

dim(expr_norm)

# =========================================================
# Save processed matrix
# =========================================================

write.table(
  expr_norm,
  file = file.path(
    dir_path,
    "validation",
    "GSE46602",
    "GSE46602_expression_normalized.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

# =========================================================
# Save phenotype data
# =========================================================

write.table(
  pheno,
  file = file.path(
    dir_path,
    "validation",
    "GSE46602",
    "GSE46602_metadata.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================================================
# Save annotation table
# =========================================================

write.table(
  annot,
  file = file.path(
    dir_path,
    "validation",
    "GSE46602",
    "GSE46602_annotation.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================================================
# Save R objects
# =========================================================

saveRDS(
  expr_norm,
  file = file.path(
    dir_path,
    "validation",
    "GSE46602",
    "GSE46602_expression_normalized.rds"
  )
)

saveRDS(
  pheno,
  file = file.path(
    dir_path,
    "validation",
    "GSE46602",
    "GSE46602_metadata.rds"
  )
)

# =========================================================
# Boxplot before normalization
# =========================================================

pdf(
  file.path(
    dir_path,
    "validation",
    "GSE46602",
    "Boxplot_Before_Normalization.pdf"
  ),
  width = 12,
  height = 6
)

boxplot(
  expr,
  outline = FALSE,
  las = 2,
  main = "GSE46602 - Before Normalization",
  col = "lightgray"
)

dev.off()

# =========================================================
# Boxplot after normalization
# =========================================================

pdf(
  file.path(
    dir_path,
    "validation",
    "GSE46602",
    "Boxplot_After_Normalization.pdf"
  ),
  width = 12,
  height = 6
)

boxplot(
  expr_norm,
  outline = FALSE,
  las = 2,
  main = "GSE46602 - After Quantile Normalization",
  col = "lightblue"
)

dev.off()

# =========================================================
# Session information
# =========================================================

sessionInfo()