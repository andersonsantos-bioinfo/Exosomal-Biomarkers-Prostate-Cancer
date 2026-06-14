# =========================================================
# 13_GSE46602_PR_analysis.R
# Precision-Recall Analysis
# LASSO Logistic Regression with Cross-Validation
# External Validation - GSE46602
# =========================================================

# =========================================================
# Load required packages
# =========================================================

library(glmnet)
library(PRROC)
library(caret)
library(ggplot2)
library(dplyr)

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
# Load expression matrix
# =========================================================

df <- read.table(
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
# Create class labels
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
# Merge expression + phenotype
# =========================================================

df <- merge(
  df,
  sample_info,
  by = "Sample"
)

# =========================================================
# Convert outcome
# =========================================================

df$Condition <- factor(
  df$Condition,
  levels = c("Normal", "Tumor")
)

# =========================================================
# Define panel genes
# =========================================================

genes_painel <- c(
  "ETV1",
  "HOXC6",
  "SLC45A2",
  "ZIC2"
)

genes_presentes <- genes_painel[
  genes_painel %in% colnames(df)
]

cat("\nGenes present:\n")
print(genes_presentes)

# =========================================================
# Stop if insufficient genes
# =========================================================

if(length(genes_presentes) < 2){
  
  stop("Insufficient genes found in validation matrix.")
  
}

# =========================================================
# Build predictor matrix
# =========================================================

x <- as.matrix(
  df[, genes_presentes]
)

# =========================================================
# Build response vector
# =========================================================

y <- ifelse(
  df$Condition == "Tumor",
  1,
  0
)

# =========================================================
# Class imbalance correction
# =========================================================

class_weights <- ifelse(
  y == 1,
  sum(y == 0) / sum(y == 1),
  1
)

cat("\nClass weights summary:\n")
print(summary(class_weights))

# =========================================================
# LASSO logistic regression with cross-validation
# =========================================================

set.seed(123)

cv_fit <- cv.glmnet(
  x = x,
  y = y,
  family = "binomial",
  alpha = 1,
  weights = class_weights,
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
# Train final LASSO model
# =========================================================

final_model <- glmnet(
  x = x,
  y = y,
  family = "binomial",
  alpha = 1,
  lambda = best_lambda,
  weights = class_weights
)

# =========================================================
# Predicted probabilities
# =========================================================

pred_prob <- predict(
  final_model,
  newx = x,
  type = "response"
)

pred_prob <- as.numeric(pred_prob)

df$Predicted_Probability <- pred_prob

# =========================================================
# Precision-Recall analysis
# =========================================================

pr <- pr.curve(
  scores.class0 = pred_prob[y == 1],
  scores.class1 = pred_prob[y == 0],
  curve = TRUE
)

# =========================================================
# Display PR-AUC
# =========================================================

cat("\n====================================\n")
cat("Precision-Recall Analysis - GSE46602\n")
cat("====================================\n")

cat("\nPR AUC:\n")
print(pr$auc.integral)

# =========================================================
# Prepare PR dataframe
# =========================================================

pr_df <- data.frame(
  Recall = pr$curve[,1],
  Precision = pr$curve[,2],
  Threshold = pr$curve[,3]
)

# =========================================================
# Save PR coordinates
# =========================================================

write.table(
  pr_df,
  file = file.path(
    dir_path,
    "results",
    "GSE46602_PR_coordinates.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================================================
# Save PR object
# =========================================================

saveRDS(
  pr,
  file = file.path(
    dir_path,
    "results",
    "GSE46602_PR_Object.rds"
  )
)

# =========================================================
# Save model
# =========================================================

saveRDS(
  final_model,
  file = file.path(
    dir_path,
    "results",
    "GSE46602_LASSO_PR_Model.rds"
  )
)

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
    "GSE46602_LASSO_Coefficients.txt"
  ),
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)

# =========================================================
# Smooth PR curve
# =========================================================

pr_smooth <- pr_df %>%
  arrange(Recall)

# =========================================================
# Create publication-quality PR plot
# =========================================================

p_pr <- ggplot(
  pr_smooth,
  aes(
    x = Recall,
    y = Precision
  )
) +
  geom_line(
    linewidth = 1.5,
    color = "#1f77b4"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Precision-Recall Curve - GSE46602",
    subtitle = paste0(
      "PR AUC = ",
      round(pr$auc.integral, 3)
    ),
    x = "Recall",
    y = "Precision"
  ) +
  coord_cartesian(
    xlim = c(0,1),
    ylim = c(0,1)
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    panel.grid.minor = element_blank()
  )

# =========================================================
# Show plot
# =========================================================

print(p_pr)

# =========================================================
# Save PDF
# =========================================================

ggsave(
  filename = file.path(
    dir_path,
    "figures",
    "GSE46602_PR_Curve.pdf"
  ),
  plot = p_pr,
  width = 6,
  height = 6
)

# =========================================================
# Save PNG
# =========================================================

ggsave(
  filename = file.path(
    dir_path,
    "figures",
    "GSE46602_PR_Curve.png"
  ),
  plot = p_pr,
  width = 6,
  height = 6,
  dpi = 300
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
    "sessionInfo_GSE46602_PR.txt"
  )
)

# =========================================================
# Final message
# =========================================================

cat("\n========================================\n")
cat("PR analysis completed successfully\n")
cat("========================================\n")