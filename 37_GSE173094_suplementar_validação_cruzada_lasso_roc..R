# =========================================================
# 37.GSE173094 - Validação Cruzada miRNAs + ROC
# =========================================================

rm(list = ls())

# =========================================================
# PACOTES
# =========================================================

libs <- c(
  "data.table",
  "caret",
  "glmnet",
  "pROC",
  "ggplot2"
)

for(i in libs){
  if(!requireNamespace(i, quietly = TRUE))
    install.packages(i)
  library(i, character.only = TRUE)
}

# =========================================================
# CARREGAR MATRIZ
# =========================================================

df <- fread(
  "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer/GSE173094/GSE173094_normalized_clean.txt",
  fill = TRUE
)

# =========================================================
# COLUNAS DAS AMOSTRAS
# =========================================================

sample_cols <- grep("^UE", colnames(df), value = TRUE)

cat("Total de amostras:", length(sample_cols), "\n")

# =========================================================
# MATRIZ DE EXPRESSÃO
# =========================================================

expr <- as.matrix(df[, ..sample_cols])

rownames(expr) <- df$ID_REF

mode(expr) <- "numeric"

dim(expr)

# =========================================================
# miRNAs DO PAINEL
# =========================================================

selected_miRNAs <- c(
  "hsa-miR-518e-002395",
  "hsa-miR-548d-5p-002237",
  "hsa-miR-142-3p-000464",
  "hsa-miR-23b-000400",
  "hsa-miR-548b-5p-002408",
  "hsa-let-7g-002282"
)

# =========================================================
# GARANTIR QUE EXISTEM
# =========================================================

selected_miRNAs <- intersect(
  selected_miRNAs,
  rownames(expr)
)

cat("miRNAs encontrados:\n")
print(selected_miRNAs)

# =========================================================
# SUBMATRIZ
# =========================================================

expr_sel <- expr[selected_miRNAs, ]

dim(expr_sel)

# =========================================================
# TRANSPOSTA
# amostras x miRNAs
# =========================================================

x <- t(expr_sel)

x <- as.data.frame(x)

dim(x)

head(x)

# =========================================================
# REMOVER COLUNAS COM NA
# =========================================================

x <- x[, colSums(is.na(x)) == 0]

# =========================================================
# DEFINIR CLASSES
# =========================================================

y <- factor(
  c(
    rep("localized", 19),
    rep("metastatic", 23)
  ),
  levels = c("localized", "metastatic")
)

length(y)

# =========================================================
# CHECAGEM CRÍTICA
# =========================================================

cat("nrow(x):", nrow(x), "\n")
cat("length(y):", length(y), "\n")

stopifnot(nrow(x) == length(y))

# =========================================================
# CONTROLE DE CV
# =========================================================

ctrl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

# =========================================================
# MODELO ELASTIC NET
# =========================================================

set.seed(123)

model_cv <- train(
  x = x,
  y = y,
  method = "glmnet",
  metric = "ROC",
  trControl = ctrl,
  preProcess = c("center", "scale"),
  tuneLength = 10
)

# =========================================================
# RESULTADOS
# =========================================================

print(model_cv)

plot(model_cv)

# =========================================================
# PREDIÇÕES
# =========================================================

pred <- model_cv$pred

best <- model_cv$bestTune

pred_best <- pred[
  pred$alpha == best$alpha &
    abs(pred$lambda - best$lambda) < 1e-10,
]

pred_best <- na.omit(pred_best)

head(pred_best)

# =========================================================
# ROC
# =========================================================

roc_cv <- roc(
  response = pred_best$obs,
  predictor = pred_best$metastatic
)

auc_val <- auc(roc_cv)

ci_val <- ci.auc(roc_cv)

cat("\nAUC:\n")
print(auc_val)

cat("\nIC95%:\n")
print(ci_val)

# =========================================================
# CURVA ROC SUAVIZADA
# =========================================================

roc_smooth <- smooth(roc_cv)

roc_df <- data.frame(
  TPR = roc_smooth$sensitivities,
  FPR = 1 - roc_smooth$specificities
)

roc_df <- roc_df[!duplicated(roc_df$FPR), ]

# =========================================================
# PONTO ÓTIMO (YOUDEN)
# =========================================================

best_coords <- coords(
  roc_cv,
  x = "best",
  best.method = "youden",
  ret = c(
    "threshold",
    "sensitivity",
    "specificity"
  ),
  transpose = FALSE
)

youden_point <- data.frame(
  FPR = 1 - as.numeric(best_coords["specificity"]),
  TPR = as.numeric(best_coords["sensitivity"])
)

print(best_coords)

# =========================================================
# PLOT ROC FINAL
# =========================================================

p <- ggplot(
  roc_df,
  aes(x = FPR, y = TPR)
) +
  
  geom_line(
    linewidth = 1,
    color = "#2C7FB8"
  ) +
  
  geom_abline(
    linetype = "dashed",
    color = "grey60"
  ) +
  
  geom_point(
    data = youden_point,
    aes(x = FPR, y = TPR),
    color = "#D7301F",
    size = 3,
    inherit.aes = FALSE
  ) +
  
  annotate(
    "text",
    x = 0.65,
    y = 0.1,
    label = paste0(
      "AUC = ",
      round(auc_val, 3),
      "\n95% CI: ",
      round(ci_val[1], 3),
      "-",
      round(ci_val[3], 3)
    ),
    size = 4
  ) +
  
  labs(
    x = "False Positive Rate",
    y = "True Positive Rate",
    title = "Cross-Validated ROC Curve"
  ) +
  
  theme_classic(base_size = 12)

print(p)

# =========================================================
# SALVAR RESULTADOS
# =========================================================

write.table(
  pred_best,
  "CV_predictions.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ggsave(
  "ROC_CV_GSE173094.png",
  p,
  width = 6,
  height = 5,
  dpi = 600
)

# =========================================================
# FIM
# =========================================================