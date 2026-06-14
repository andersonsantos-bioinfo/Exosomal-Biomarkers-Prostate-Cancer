# ================================================
# Validação Robusta - 5-fold Repeated CV (1000 repetições)
# GSE46602 - Painel de 4 genes
# ================================================

library(caret)
library(pROC)
library(PRROC)
library(tidyverse)

setwd("C:/Users/Administrator/Exossomos_prostata/Fase-Final/Validação_externa/GSE46602")

# Ler os dados
df <- read.csv("dados_modelo_4genes.csv", stringsAsFactors = FALSE, check.names = FALSE)

# Calcular Risk Score
genes <- c("ETV1", "HOXC6", "SLC45A2", "ZIC2")
df$risk_score <- rowSums(scale(df[, genes]))

# Preparar dados para modelagem
df$Class <- factor(df$group, levels = c("Normal", "Tumor"))

cat("Resumo das classes:\n")
print(table(df$Class))
cat("\n")

# ====================== CONFIGURAÇÃO DA VALIDAÇÃO ======================
set.seed(1234)  # para reprodutibilidade

train_control <- trainControl(
  method = "repeatedcv",      # repeated cross-validation
  number = 5,                 # 5-fold
  repeats = 1000,             # 1000 repetições
  classProbs = TRUE,          # calcular probabilidades
  summaryFunction = twoClassSummary,  # para métricas de classificação binária
  savePredictions = "final",
  verboseIter = FALSE
)

# ====================== TREINAMENTO COM REPEATED CV ======================
cat("Iniciando 5-fold Repeated CV (1000 repetições)...\n")

model_cv <- train(
  x = df[, genes],           # apenas os 4 genes
  y = df$Class,
  method = "glmnet",         # Regressão logística com regularização (elastic net)
  family = "binomial",
  trControl = train_control,
  tuneLength = 10,           # testar vários valores de alpha e lambda
  metric = "ROC"
)

cat("Validação cruzada concluída!\n\n")

# ====================== RESULTADOS PRINCIPAIS ======================
cat("=== RESULTADOS DA VALIDAÇÃO 5-FOLD REPETIDA 1000x ===\n")

# Melhor modelo
best_model <- model_cv$results[which.max(model_cv$results$ROC), ]
cat("Melhor AUC médio:", round(best_model$ROC, 4), "\n")
cat("Desvio padrão do AUC:", round(best_model$ROCSD, 4), "\n\n")

# Predições out-of-fold
pred <- model_cv$pred
pred$obs_bin <- ifelse(pred$obs == "Tumor", 1, 0)

# Métricas agregadas
auc_mean <- mean(pred$ROC)   # se quiser média das folds
cat("AUC médio aproximado:", round(auc_mean, 4), "\n")

# ====================== SALVAR RESULTADOS ======================
resultados_cv <- model_cv$results
write.table(resultados_cv, 
            "GSE46602_5foldCV_1000repeticoes_Resultados.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nResultados detalhados salvos em: GSE46602_5foldCV_1000repeticoes_Resultados.txt\n")
cat("Melhor configuração encontrada:\n")
print(best_model)