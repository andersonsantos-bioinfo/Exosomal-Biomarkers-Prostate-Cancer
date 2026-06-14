# =========================================================
# 06_ROC_PR_analysis.R
# ROC and Precision-Recall Analysis
# Exosomal Biomarkers in Prostate Cancer
# =========================================================

# =========================================================
# Load libraries
# =========================================================

library(pROC)
library(PRROC)
library(ggplot2)
library(dplyr)

# =========================================================
# Define project directory
# =========================================================

dir_path <- "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer"

# =========================================================
# Create performance directory
# =========================================================

dir.create(
  file.path(dir_path, "performance"),
  showWarnings = FALSE
)

# =========================================================
# Load Elastic Net predictions
# =========================================================

pred_data <- read.table(
  file.path(
    dir_path,
    "models",
    "ElasticNet_Predictions.txt"
  ),
  header = TRUE,
  sep = "\t"
)

# =========================================================
# Inspect data
# =========================================================

head(pred_data)

dim(pred_data)

# =========================================================
# Define outcome
# =========================================================

pred_data$Condition <- as.factor(
  pred_data$Condition
)

table(pred_data$Condition)

# =========================================================
# Create binary labels
# =========================================================

labels <- ifelse(
  pred_data$Condition == levels(pred_data$Condition)[2],
  1,
  0
)

# =========================================================
# Predicted probabilities
# =========================================================

scores <- pred_data$Predicted_Probability

# =========================================================
# ROC analysis
# =========================================================

roc_obj <- roc(
  response = labels,
  predictor = scores
)

# =========================================================
# Calculate AUC
# =========================================================

roc_auc <- auc(roc_obj)

roc_auc

# =========================================================
# Confidence interval
# =========================================================

roc_ci <- ci.auc(roc_obj)

roc_ci

# =========================================================
# Determine optimal cutoff
# =========================================================

best_cutoff <- coords(
  roc_obj,
  "best",
  ret = c(
    "threshold",
    "sensitivity",
    "specificity"
  )
)

best_cutoff

# =========================================================
# ROC curve plot
# =========================================================

roc_plot <- ggroc(
  roc_obj,
  linewidth = 1.2
) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "gray60"
  ) +
  
  labs(
    title = "ROC Curve - Elastic Net Model",
    subtitle = paste(
      "AUC =",
      round(roc_auc, 3)
    ),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )

# =========================================================
# Show ROC curve
# =========================================================

print(roc_plot)

# =========================================================
# Save ROC plot
# =========================================================

ggsave(
  filename = file.path(
    dir_path,
    "performance",
    "ROC_Curve_ElasticNet.pdf"
  ),
  plot = roc_plot,
  width = 7,
  height = 6
)

# =========================================================
# Precision-Recall analysis
# =========================================================

pr_obj <- pr.curve(
  scores.class0 = scores[labels == 1],
  scores.class1 = scores[labels == 0],
  curve = TRUE
)

# =========================================================
# PR AUC
# =========================================================

pr_auc <- pr_obj$auc.integral

pr_auc

# =========================================================
# Prepare PR data
# =========================================================

pr_data <- data.frame(
  Recall = pr_obj$curve[,1],
  Precision = pr_obj$curve[,2]
)

# =========================================================
# Create PR curve
# =========================================================

pr_plot <- ggplot(
  pr_data,
  aes(
    x = Recall,
    y = Precision
  )
) +
  
  geom_line(
    linewidth = 1.2
  ) +
  
  labs(
    title = "Precision-Recall Curve",
    subtitle = paste(
      "AUC-PR =",
      round(pr_auc, 3)
    ),
    x = "Recall",
    y = "Precision"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )

# =========================================================
# Show PR curve
# =========================================================

print(pr_plot)

# =========================================================
# Save PR plot
# =========================================================

ggsave(
  filename = file.path(
    dir_path,
    "performance",
    "PR_Curve_ElasticNet.pdf"
  ),
  plot = pr_plot,
  width = 7,
  height = 6
)

# =========================================================
# Save performance metrics
# =========================================================

performance_metrics <- data.frame(
  Metric = c(
    "ROC_AUC",
    "ROC_CI_Lower",
    "ROC_CI_Upper",
    "PR_AUC",
    "Optimal_Threshold",
    "Sensitivity",
    "Specificity"
  ),
  
  Value = c(
    as.numeric(roc_auc),
    as.numeric(roc_ci[1]),
    as.numeric(roc_ci[3]),
    as.numeric(pr_auc),
    as.numeric(best_cutoff["threshold"]),
    as.numeric(best_cutoff["sensitivity"]),
    as.numeric(best_cutoff["specificity"])
  )
)

write.table(
  performance_metrics,
  file = file.path(
    dir_path,
    "performance",
    "Performance_Metrics.txt"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# =========================================================
# Session information
# =========================================================

sessionInfo()