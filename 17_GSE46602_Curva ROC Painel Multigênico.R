library(pROC)
library(ggplot2)

# Preparar dados
risk_df$group_bin <- ifelse(risk_df$group == "Tumor", 1, 0)

# ROC
roc_obj <- roc(risk_df$group_bin, risk_df$risk_score)

# Converter para dataframe
roc_df <- data.frame(
  fpr = 1 - roc_obj$specificities,
  tpr = roc_obj$sensitivities
)

# ???? ORDENAR corretamente (remove artefatos)
roc_df <- roc_df[order(roc_df$fpr, roc_df$tpr), ]

# Plot elegante (sem artefatos)
p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
  
  # ???? geom_step evita linhas verticais erradas
  geom_step(color = "#E64B35", linewidth = 1) +
  
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "gray60", linewidth = 0.6) +
  
  labs(
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  
  theme_classic(base_size = 10, base_family = "Arial") +
  
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank()
  )

print(p_roc)