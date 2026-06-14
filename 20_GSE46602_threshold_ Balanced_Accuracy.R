# Threshold escolhido
thr_final <- -1.517435

# Aplicar classificação
risk_df$pred_final <- ifelse(risk_df$risk_score >= thr_final, 1, 0)

# Matriz de confusão
cm <- table(
  factor(risk_df$pred_final, levels = c(0,1)),
  factor(risk_df$group_bin, levels = c(0,1))
)

cm

# Métricas
TN <- cm[1,1]
FP <- cm[2,1]
FN <- cm[1,2]
TP <- cm[2,2]

sens <- TP/(TP+FN)
spec <- TN/(TN+FP)
bal_acc <- (sens + spec)/2

precision <- TP/(TP+FP)
f1 <- 2 * (precision * sens) / (precision + sens)

list(
  Sensibilidade = sens,
  Especificidade = spec,
  Balanced_Accuracy = bal_acc,
  F1 = f1
)