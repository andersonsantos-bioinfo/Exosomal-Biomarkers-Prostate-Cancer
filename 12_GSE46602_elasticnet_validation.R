# =========================================================
# 12_GSE46602_elasticnet_validation.R
# Elastic Net Validation in External Cohort GSE46602
# =========================================================

# =========================================================
# Load required packages
# =========================================================

library(glmnet)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# =========================================================
# Define project directory
# =========================================================

dir_path <- "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer"

# =========================================================
# Create folders
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
# Load external validation matrix
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

# =========================================================
# Load GEO metadata
# =========================================================

library(GEOquery)
library(Biobase)

gse <- getGEO(
  "GSE46602",
  GSEMatrix = TRUE
)

gse <- gse[[1]]

pheno <- pData(gse)

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

sample_info <- data.frame(
  Sample = rownames(pheno),
  Condition = condition,
  stringsAsFactors = FALSE
)

# =========================================================
# Merge metadata
# =========================================================

expr_panel <- merge(
  expr_panel,
  sample_info,
  by = "Sample"
)

# =========================================================
# Convert outcome
# =========================================================

expr_panel$Condition <- factor(
  expr_panel$Condition,
  levels = c("Normal", "Tumor")
)

# =========================================================
# Define genes
# =========================================================

genes_painel <- c(
  "ETV1",
  "HOXC6",
  "SLC45A2",
  "ZIC2"
)

genes_presentes <- genes_painel[
  genes_painel %in% colnames(expr_panel)
]

cat("\nGenes present:\n")
print(genes_presentes)

# =========================================================
# Stop if genes missing
# =========================================================

if(length(genes_presentes) < 2){
  
  stop("Insufficient panel genes in validation matrix.")
  
}

# =========================================================
# Build predictor matrix
# =========================================================

x <- as.matrix(
  expr_panel[, genes_presentes]
)

# =========================================================
# Build response
# =========================================================

y <- ifelse(
  expr_panel$Condition == "Tumor",
  1,
  0
)

# =========================================================
# Train Elastic Net model
# =========================================================

set.seed(123)

cv_fit <- cv.glmnet(
  x,
  y,
  family = "binomial",
  alpha = 0.5,
  nfolds = 5,
  type.measure = "auc"
)

# =========================================================
# Best lambda
# =========================================================

best_lambda <- cv_fit$lambda.min

cat("\nBest lambda:\n")
print(best_lambda)

# =========================================================
# Final model
# =========================================================

final_model <- glmnet(
  x,
  y,
  family = "binomial",
  alpha = 0.5,
  lambda = best_lambda
)

# =========================================================
# Predictions
# =========================================================

pred_prob <- predict(
  final_model,
  newx = x,
  type = "response"
)

pred_prob <- as.numeric(pred_prob)

expr_panel$Predicted_Probability <- pred_prob

# =========================================================
# ROC analysis
# =========================================================

roc_obj <- roc(
  response = expr_panel$Condition,
  predictor = expr_panel$Predicted_Probability,
  levels = c("Normal", "Tumor"),
  direction = "<"
)

auc_value <- auc(roc_obj)

ci_auc <- ci.auc(roc_obj)

cat("\n====================================\n")
cat("Elastic Net External Validation\n")
cat("====================================\n")

cat("\nAUC:\n")
print(auc_value)

cat("\n95% CI:\n")
print(ci_auc)

# =========================================================
# Optimal threshold
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

threshold <- as.numeric(
  best_coords["threshold"]
)

cat("\nOptimal threshold:\n")
print(threshold)

# =========================================================
# Predicted classes
# =========================================================

expr_panel$Predicted_Class <- ifelse(
  expr_panel$Predicted_Probability >= threshold,
  "Tumor",
  "Normal"
)

expr_panel$Predicted_Class <- factor(
  expr_panel$Predicted_Class,
  levels = c("Normal", "Tumor")
)

# =========================================================
# Confusion matrix
# =========================================================

cm <- confusionMatrix(
  expr_panel$Predicted_Class,
  expr_panel$Condition,
  positive = "Tumor"
)

cat("\nConfusion Matrix:\n")
print(cm)

# =========================================================
# Save predictions
# =========================================================

write.table(
  expr_panel,
  file = file.path(
    dir_path,
    "results",
    "GSE46602_ElasticNet_Validation.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================================================
# Save model
# =========================================================

saveRDS(
  final_model,
  file = file.path(
    dir_path,
    "results",
    "GSE46602_ElasticNet_Model.rds"
  )
)

# =========================================================
# Save ROC object
# =========================================================

saveRDS(
  roc_obj,
  file = file.path(
    dir_path,
    "results",
    "GSE46602_ElasticNet_ROC.rds"
  )
)

# =========================================================
# ROC plot PDF
# =========================================================

pdf(
  file.path(
    dir_path,
    "figures",
    "GSE46602_ElasticNet_ROC.pdf"
  ),
  width = 6,
  height = 6
)

plot(
  roc_obj,
  col = "#1f77b4",
  lwd = 3,
  legacy.axes = TRUE,
  main = "Elastic Net Validation - GSE46602",
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
# ROC plot PNG
# =========================================================

png(
  file.path(
    dir_path,
    "figures",
    "GSE46602_ElasticNet_ROC.png"
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
  main = "Elastic Net Validation - GSE46602",
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
# Save coefficients
# =========================================================

coef_df <- as.data.frame(
  as.matrix(
    coef(final_model)
  )
)

write.table(
  coef_df,
  file = file.path(
    dir_path,
    "results",
    "GSE46602_ElasticNet_Coefficients.txt"
  ),
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)

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
    "sessionInfo_GSE46602_elasticnet.txt"
  )
)

# =========================================================
# Final message
# =========================================================

cat("\n========================================\n")
cat("Elastic Net validation completed successfully\n")
cat("========================================\n")