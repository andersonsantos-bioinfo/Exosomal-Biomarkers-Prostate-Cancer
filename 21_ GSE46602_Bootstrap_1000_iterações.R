library(pROC)

set.seed(123)

# Número de iterações
n_boot <- 1000

# Garantir variável binária
risk_df$group_bin <- ifelse(risk_df$group == "Tumor", 1, 0)

# ???? usar o threshold já definido anteriormente
thr_final <- -1.517435

# Armazenar resultados
boot_results <- data.frame(
  AUC = numeric(n_boot),
  Sens = numeric(n_boot),
  Spec = numeric(n_boot)
)

# ---------------------------------------------------------
# LOOP BOOTSTRAP
# ---------------------------------------------------------
for (i in 1:n_boot) {
  
  # 1. Reamostragem com reposição
  idx <- sample(1:nrow(risk_df), replace = TRUE)
  df_boot <- risk_df[idx, ]
  
  # 2. ROC + AUC
  roc_obj <- roc(df_boot$group_bin, df_boot$risk_score)
  boot_results$AUC[i] <- as.numeric(auc(roc_obj))
  
  # 3. Predição com threshold fixo
  pred <- ifelse(df_boot$risk_score >= thr_final, 1, 0)
  
  # 4. Matriz de confusão com níveis fixos
  cm <- table(
    factor(pred, levels = c(0,1)),
    factor(df_boot$group_bin, levels = c(0,1))
  )
  
  TN <- cm[1,1]
  FP <- cm[2,1]
  FN <- cm[1,2]
  TP <- cm[2,2]
  
  # 5. Métricas
  sens <- ifelse((TP+FN)==0, NA, TP/(TP+FN))
  spec <- ifelse((TN+FP)==0, NA, TN/(TN+FP))
  
  boot_results$Sens[i] <- sens
  boot_results$Spec[i] <- spec
}

# ---------------------------------------------------------
# RESUMO (MÉDIA + IC95%)
# ---------------------------------------------------------
summary_boot <- data.frame(
  Metrica = c("AUC", "Sensibilidade", "Especificidade"),
  
  Media = c(
    mean(boot_results$AUC, na.rm = TRUE),
    mean(boot_results$Sens, na.rm = TRUE),
    mean(boot_results$Spec, na.rm = TRUE)
  ),
  
  IC_inf = c(
    quantile(boot_results$AUC, 0.025, na.rm = TRUE),
    quantile(boot_results$Sens, 0.025, na.rm = TRUE),
    quantile(boot_results$Spec, 0.025, na.rm = TRUE)
  ),
  
  IC_sup = c(
    quantile(boot_results$AUC, 0.975, na.rm = TRUE),
    quantile(boot_results$Sens, 0.975, na.rm = TRUE),
    quantile(boot_results$Spec, 0.975, na.rm = TRUE)
  )
)

print(summary_boot)