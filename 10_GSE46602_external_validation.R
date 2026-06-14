# =========================================================
# 11_GSE46602_external_validation.R
# External Validation of Multigene Panel
# GSE46602 Cohort
# =========================================================

# =========================================================
# Load required packages
# =========================================================

library(GEOquery)
library(Biobase)
library(dplyr)
library(pROC)
library(caret)
library(ggplot2)

# =========================================================
# Define project directory
# =========================================================

dir_path <- "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer"

# =========================================================
# Create output folders
# =========================================================

dir.create(
  file.path(dir_path, "results"),
  showWarnings = FALSE
)

dir.create(
  file.path(dir_path, "figures"),
  showWarnings = FALSE
)

# =========================================================
# Load panel expression matrix
# =========================================================

expr_panel <- read.table(
  file.path(
    dir_path,
    "validation",
    "GSE46602",
    "GSE46602_Panel_Expression.txt"
  ),
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

cat("\nExpression matrix loaded:\n")
print(dim(expr_panel))

cat("\nColumns:\n")
print(colnames(expr_panel))

# =========================================================
# Load GEO metadata
# =========================================================

gse <- getGEO(
  "GSE46602",
  GSEMatrix = TRUE
)

gse <- gse[[1]]

# =========================================================
# Extract phenotype data
# =========================================================

pheno <- pData(gse)

cat("\nPhenotype columns:\n")
print(colnames(pheno))

# =========================================================
# Inspect source names
# =========================================================

if("source_name_ch1" %in% colnames(pheno)){
  
  cat("\nSource names:\n")
  print(table(pheno$source_name_ch1))
  
}

# =========================================================
# Create condition labels
# =========================================================

condition <- ifelse(
  grepl(
    "tumor|cancer|primary",
    pheno$source_name_ch1,
    ignore.case = TRUE
  ),
  "Tumor",
  "Normal"
)

# =========================================================
# Create sample information dataframe
# =========================================================

sample_info <- data.frame(
  Sample = rownames(pheno),
  Condition = condition,
  stringsAsFactors = FALSE
)

# =========================================================
# Merge expression + phenotype
# =========================================================

expr_panel <- merge(
  expr_panel,
  sample_info,
  by = "Sample"
)

cat("\nMerged dataset dimensions:\n")
print(dim(expr_panel))

# =========================================================
# Verify class distribution
# =========================================================

cat("\nCondition distribution:\n")
print(table(expr_panel$Condition))

# =========================================================
# Convert condition to factor
# =========================================================

expr_panel$Condition <- factor(
  expr_panel$Condition,
  levels = c("Normal", "Tumor")
)

# =========================================================
# Define TCGA-derived coefficients
# Replace with your final coefficients if needed
# =========================================================

coef_4genes <- c(
  ETV1 = 0.85,
  HOXC6 = 1.10,
  SLC45A2 = 0.72,
  ZIC2 = 0.66
)

intercept <- -2.15

# =========================================================
# Check genes present
# =========================================================

genes_modelo <- names(coef_4genes)

genes_presentes <- genes_modelo[
  genes_modelo %in% colnames(expr_panel)
]

cat("\nGenes present in validation matrix:\n")
print(genes_presentes)

# =========================================================
# Stop if no genes found
# =========================================================

if(length(genes_presentes) == 0){
  
  stop("No model genes found in validation matrix.")
  
}

# =========================================================
# Keep only available coefficients
# =========================================================

coef_4genes <- coef_4genes[
  genes_presentes
]

# =========================================================
# Calculate risk score
# =========================================================

risk_score <- intercept

for(gene in genes_presentes){
  
  risk_score <- risk_score +
    (
      expr_panel[[gene]] *
        coef_4genes[gene]
    )
  
}

# =========================================================
# Add risk score
# =========================================================

expr_panel$Risk_Score <- risk_score

# =========================================================
# Calculate predicted probability
# =========================================================

expr_panel$Predicted_Probability <- 1 / (
  1 + exp(-risk_score)
)

# =========================================================
# Verify probabilities
# =========================================================

cat("\nProbability summary:\n")
print(summary(expr_panel$Predicted_Probability))

# =========================================================
# ROC analysis
# =========================================================

roc_obj <- roc(
  response = expr_panel$Condition,
  predictor = expr_panel$Predicted_Probability,
  levels = c("Normal", "Tumor"),
  direction = "<"
)

# =========================================================
# Calculate AUC
# =========================================================

auc_value <- auc(roc_obj)

ci_auc <- ci.auc(roc_obj)

cat("\n=====================================\n")
cat("External Validation - GSE46602\n")
cat("=====================================\n")

cat("\nAUC:\n")
print(auc_value)

cat("\n95% CI:\n")
print(ci_auc)

# =========================================================
# Determine optimal threshold
# =========================================================

best_coords <- coords(
  roc_obj,
  "best",
  ret = c(
    "threshold",
    "sensitivity",
    "specificity"
  ),
  best.method = "youden"
)

cat("\nOptimal coordinates:\n")
print(best_coords)

# =========================================================
# Extract numeric threshold
# =========================================================

threshold <- as.numeric(
  best_coords["threshold"]
)

cat("\nOptimal threshold:\n")
print(threshold)

# =========================================================
# Generate predicted classes
# =========================================================

expr_panel$Predicted_Class <- ifelse(
  expr_panel$Predicted_Probability >= threshold,
  "Tumor",
  "Normal"
)

# =========================================================
# Convert to factor
# =========================================================

expr_panel$Predicted_Class <- factor(
  expr_panel$Predicted_Class,
  levels = c("Normal", "Tumor")
)

# =========================================================
# Verify predictions
# =========================================================

cat("\nPredicted class distribution:\n")
print(table(expr_panel$Predicted_Class))

# =========================================================
# Confusion matrix
# =========================================================

cm <- confusionMatrix(
  expr_panel$Predicted_Class,
  expr_panel$Condition,
  positive = "Tumor"
)

cat("\nConfusion matrix:\n")
print(cm)

# =========================================================
# Save validation results
# =========================================================

write.table(
  expr_panel,
  file = file.path(
    dir_path,
    "results",
    "GSE46602_Validation_Results.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================================================
# Save ROC coordinates
# =========================================================

roc_df <- data.frame(
  Sensitivity = roc_obj$sensitivities,
  Specificity = roc_obj$specificities,
  Threshold = roc_obj$thresholds
)

write.table(
  roc_df,
  file = file.path(
    dir_path,
    "results",
    "GSE46602_ROC_coordinates.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================================================
# Save ROC object
# =========================================================

saveRDS(
  roc_obj,
  file = file.path(
    dir_path,
    "results",
    "GSE46602_ROC_Object.rds"
  )
)

# =========================================================
# ROC Plot - PDF
# =========================================================

pdf(
  file.path(
    dir_path,
    "figures",
    "GSE46602_ROC_Curve.pdf"
  ),
  width = 6,
  height = 6
)

plot(
  roc_obj,
  col = "#1f77b4",
  lwd = 3,
  legacy.axes = TRUE,
  main = "External Validation - GSE46602",
  print.auc = TRUE,
  print.auc.cex = 1.2
)

abline(
  a = 0,
  b = 1,
  lty = 2,
  col = "gray60"
)

dev.off()

# =========================================================
# ROC Plot - PNG
# =========================================================

png(
  file.path(
    dir_path,
    "figures",
    "GSE46602_ROC_Curve.png"
  ),
  width = 2000,
  height = 2000,
  res = 300
)

plot(
  roc_obj,
  col = "#1f77b4",
  lwd = 3,
  legacy.axes = TRUE,
  main = "External Validation - GSE46602",
  print.auc = TRUE,
  print.auc.cex = 1.2
)

abline(
  a = 0,
  b = 1,
  lty = 2,
  col = "gray60"
)

dev.off()

# =========================================================
# Save session info
# =========================================================

dir.create(
  file.path(dir_path, "session_info"),
  showWarnings = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    dir_path,
    "session_info",
    "sessionInfo_GSE46602_validation.txt"
  )
)

# =========================================================
# Final message
# =========================================================

cat("\n========================================\n")
cat("External validation completed successfully\n")
cat("========================================\n")