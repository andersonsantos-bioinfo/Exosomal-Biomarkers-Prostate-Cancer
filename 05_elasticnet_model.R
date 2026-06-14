# =========================================================
# 05_elasticnet_model.R
# Elastic Net Model - TCGA PRAD
# Exosomal Biomarkers in Prostate Cancer
# =========================================================

# =========================================================
# Load libraries
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
# Create model directory
# =========================================================

dir.create(
  file.path(dir_path, "models"),
  showWarnings = FALSE
)

# =========================================================
# Load expression matrix
# =========================================================

expr_data <- read.table(
  file.path(
    dir_path,
    "data_processed",
    "Matriz_Painel_6genes_com_Condition.txt"
  ),
  header = TRUE,
  sep = "\t",
  row.names = 1
)

# =========================================================
# Inspect data
# =========================================================

dim(expr_data)

head(expr_data)

# =========================================================
# Define response variable
# =========================================================

y <- as.factor(expr_data$Condition)

table(y)

# =========================================================
# Create predictor matrix
# =========================================================

x <- expr_data %>%
  dplyr::select(-Condition)

x <- as.matrix(x)

# =========================================================
# Standardize predictors
# =========================================================

x <- scale(x)

# =========================================================
# Cross-validation Elastic Net
# =========================================================

set.seed(123)

cv_model <- cv.glmnet(
  x = x,
  y = y,
  family = "binomial",
  alpha = 0.5,
  nfolds = 10,
  type.measure = "auc"
)

# =========================================================
# Best lambda
# =========================================================

best_lambda <- cv_model$lambda.min

best_lambda

# =========================================================
# Train final Elastic Net model
# =========================================================

elastic_model <- glmnet(
  x = x,
  y = y,
  family = "binomial",
  alpha = 0.5,
  lambda = best_lambda
)

# =========================================================
# Extract coefficients
# =========================================================

coef_table <- as.matrix(
  coef(elastic_model)
)

coef_table

# =========================================================
# Save coefficients
# =========================================================

write.table(
  coef_table,
  file = file.path(
    dir_path,
    "models",
    "ElasticNet_Coefficients.txt"
  ),
  sep = "\t",
  quote = FALSE
)

# =========================================================
# Predict probabilities
# =========================================================

pred_prob <- predict(
  elastic_model,
  newx = x,
  type = "response"
)

pred_prob <- as.numeric(pred_prob)

# =========================================================
# ROC analysis
# =========================================================

roc_obj <- roc(
  response = y,
  predictor = pred_prob
)

# =========================================================
# AUC value
# =========================================================

auc_value <- auc(roc_obj)

auc_value

# =========================================================
# Plot ROC curve
# =========================================================

roc_plot <- ggroc(
  roc_obj,
  linewidth = 1.2
) +
  
  labs(
    title = "Elastic Net ROC Curve",
    subtitle = paste(
      "AUC =",
      round(auc_value, 3)
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
    "models",
    "ElasticNet_ROC_Curve.pdf"
  ),
  plot = roc_plot,
  width = 7,
  height = 6
)

# =========================================================
# Save trained model
# =========================================================

saveRDS(
  elastic_model,
  file = file.path(
    dir_path,
    "models",
    "ElasticNet_Model.rds"
  )
)

# =========================================================
# Save predictions
# =========================================================

prediction_table <- data.frame(
  Sample = rownames(expr_data),
  Condition = y,
  Predicted_Probability = pred_prob
)

write.table(
  prediction_table,
  file = file.path(
    dir_path,
    "models",
    "ElasticNet_Predictions.txt"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# =========================================================
# Session information
# =========================================================

sessionInfo()