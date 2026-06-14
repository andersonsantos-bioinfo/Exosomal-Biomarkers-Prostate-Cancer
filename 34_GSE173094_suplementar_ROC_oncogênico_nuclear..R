# =========================================================
# GSE173094 - ROC do programa oncogênico nuclear
# HOXC6 . ETV1 . ZIC2 . SLC45A2
# =========================================================

library(data.table)
library(pROC)
library(ggplot2)

# =========================================================
# 1. CARREGAR MATRIZ DE EXPRESSÃO
# =========================================================

expr_file <- "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer/GSE173094/GSE173094_normalized_clean.txt"

df <- fread(
  expr_file,
  fill = TRUE,
  data.table = FALSE
)

# =========================================================
# 2. MATRIZ DE EXPRESSÃO
# miRNAs x amostras
# =========================================================

expr <- as.matrix(
  df[, grep("^UE", colnames(df))]
)

rownames(expr) <- df$ID_REF

# =========================================================
# 3. DEFINIR miRNAs DO EIXO
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
# 4. GARANTIR QUE EXISTEM
# =========================================================

selected_miRNAs <- intersect(
  selected_miRNAs,
  rownames(expr)
)

cat(
  "miRNAs encontradas:",
  length(selected_miRNAs),
  "\n"
)

print(selected_miRNAs)

# =========================================================
# 5. PADRONIZAÇÃO (Z-SCORE)
# =========================================================

expr_z <- t(
  scale(
    t(expr)
  )
)

# =========================================================
# 6. SCORE ONCOGÊNICO
# média das miRNAs por amostra
# =========================================================

score <- colMeans(
  expr_z[selected_miRNAs, ],
  na.rm = TRUE
)

# =========================================================
# 7. DEFINIR GRUPOS
# =========================================================

group <- factor(
  c(
    rep("localized", 19),
    rep("metastatic", 23)
  ),
  levels = c(
    "localized",
    "metastatic"
  )
)

table(group)

# =========================================================
# 8. ROC
# =========================================================

roc_obj <- roc(
  response = group,
  predictor = score,
  direction = "auto"
)

# =========================================================
# 9. AUC + IC95%
# =========================================================

auc_value <- auc(roc_obj)

ci_auc <- ci.auc(roc_obj)

print(auc_value)

print(ci_auc)

# =========================================================
# 10. PONTO ÓTIMO (YOUDEN)
# =========================================================

best_coords <- coords(
  roc_obj,
  "best",
  ret = c(
    "threshold",
    "sensitivity",
    "specificity"
  )
)

print(best_coords)

# =========================================================
# 11. DATAFRAME ROC
# =========================================================

roc_df <- data.frame(
  specificity = roc_obj$specificities,
  sensitivity = roc_obj$sensitivities
)

roc_df$fpr <- 1 - roc_df$specificity

# remover duplicatas
roc_df <- roc_df[
  !duplicated(roc_df$fpr),
]

# =========================================================
# 12. PONTO DE YOUDEN
# =========================================================

youden_df <- data.frame(
  fpr = 1 - as.numeric(best_coords["specificity"]),
  sensitivity = as.numeric(best_coords["sensitivity"])
)

print(youden_df)

# =========================================================
# 13. CURVA ROC
# =========================================================

p_roc <- ggplot(
  roc_df,
  aes(
    x = fpr,
    y = sensitivity
  )
) +
  geom_line(
    color = "#D62728",
    linewidth = 1
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "gray60",
    linewidth = 0.5
  ) +
  geom_point(
    data = youden_df,
    mapping = aes(
      x = fpr,
      y = sensitivity
    ),
    inherit.aes = FALSE,
    color = "red",
    size = 2.8
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  labs(
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_classic(
    base_family = "Arial",
    base_size = 10
  ) +
  theme(
    axis.text = element_text(
      color = "black",
      face = "plain"
    ),
    axis.title = element_text(
      color = "black",
      face = "plain"
    )
  )

print(p_roc)

# =========================================================
# 14. SALVAR FIGURA
# =========================================================

ggsave(
  filename = "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer/GSE173094/ROC_programa_oncogenico_GSE173094.png",
  plot = p_roc,
  width = 5,
  height = 5,
  dpi = 600
)

# =========================================================
# 15. EXPORTAR RESULTADOS
# =========================================================

resultados <- data.frame(
  Sample = colnames(expr),
  Group = group,
  Score = score
)

write.table(
  resultados,
  file = "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer/GSE173094/oncogenic_program_score.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# =========================================================
# 16. RESUMO FINAL
# =========================================================

cat("\n========================================\n")
cat("AUC =", round(as.numeric(auc_value), 3), "\n")
cat("IC95% =", round(ci_auc[1], 3), "-",
    round(ci_auc[3], 3), "\n")
cat("========================================\n")