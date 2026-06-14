# =========================================================
# GSE173094 - Programa oncogênico nuclear
# HOXC6 . ETV1 . ZIC2 . SLC45A2
# =========================================================

library(data.table)
library(ggplot2)
library(pROC)
library(pheatmap)

# =========================================================
# 1. CARREGAR DADOS
# =========================================================

expr_file <- "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer/GSE173094/GSE173094_normalized_clean.txt"

df <- fread(
  expr_file,
  fill = TRUE,
  data.table = FALSE
)

# =========================================================
# 2. IDENTIFICAR AMOSTRAS
# =========================================================

sample_cols <- grep("^UE", colnames(df), value = TRUE)

cat("Encontradas", length(sample_cols), "amostras\n")

# =========================================================
# 3. MATRIZ DE EXPRESSÃO
# =========================================================

expr <- as.matrix(df[, sample_cols])

rownames(expr) <- df$ID_REF

dim(expr)

# =========================================================
# 4. DEFINIR GRUPOS
# =========================================================

group <- factor(
  c(
    rep("Localized", 19),
    rep("Metastatic", 23)
  ),
  levels = c("Localized", "Metastatic")
)

table(group)

# =========================================================
# 5. miRNAs DO EIXO ONCOGÊNICO
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
# 6. VERIFICAR miRNAs DISPONÍVEIS
# =========================================================

present <- selected_miRNAs %in% rownames(expr)

check_df <- data.frame(
  miRNA = selected_miRNAs,
  Present = present
)

print(check_df)

# =========================================================
# 7. EXTRAIR EXPRESSÃO
# =========================================================

expr_miRNA <- expr[selected_miRNAs[present], ]

expr_miRNA <- as.matrix(expr_miRNA)

dim(expr_miRNA)

# =========================================================
# 8. EXPRESSÃO RELATIVA
# =========================================================

expr_rel <- 40 - expr_miRNA

# =========================================================
# 9. PADRONIZAÇÃO
# =========================================================

expr_scaled <- t(scale(t(expr_rel)))

# =========================================================
# 10. SCORE ONCOGÊNICO
# =========================================================

oncogenic_score <- colMeans(
  expr_scaled,
  na.rm = TRUE
)

# =========================================================
# 11. DATAFRAME FINAL
# =========================================================

result_df <- data.frame(
  Sample = colnames(expr_scaled),
  Group = group,
  Oncogenic_Score = oncogenic_score
)

head(result_df)

# =========================================================
# 12. TESTE ESTATÍSTICO
# =========================================================

ttest_res <- t.test(
  oncogenic_score ~ group
)

print(ttest_res)

# =========================================================
# 13. MODELO LINEAR
# =========================================================

model <- lm(
  oncogenic_score ~ group
)

summary(model)

# =========================================================
# 14. BOXPLOT
# =========================================================

p_box <- ggplot(
  result_df,
  aes(
    x = Group,
    y = Oncogenic_Score,
    fill = Group
  )
) +
  geom_boxplot(
    width = 0.65,
    alpha = 0.85,
    outlier.shape = 16,
    outlier.size = 2
  ) +
  geom_jitter(
    width = 0.08,
    size = 2,
    alpha = 0.8
  ) +
  scale_fill_manual(
    values = c(
      "Localized" = "#5DA5DA",
      "Metastatic" = "#F15854"
    )
  ) +
  labs(
    x = "",
    y = "Oncogenic Axis Score"
  ) +
  theme_classic(
    base_family = "Arial",
    base_size = 10
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(
      color = "black",
      face = "plain"
    ),
    axis.title = element_text(
      color = "black",
      face = "plain"
    )
  )

print(p_box)

# =========================================================
# 15. ROC
# =========================================================

roc_obj <- roc(
  response = group,
  predictor = oncogenic_score,
  levels = c("Localized", "Metastatic")
)

print(auc(roc_obj))

print(ci.auc(roc_obj))

# =========================================================
# 16. ROC SUAVIZADA
# =========================================================

roc_smooth <- smooth(roc_obj)

roc_df <- data.frame(
  TPR = roc_smooth$sensitivities,
  FPR = 1 - roc_smooth$specificities
)

roc_df <- roc_df[
  !duplicated(roc_df$FPR),
]

head(roc_df)

# =========================================================
# 17. PONTO DE YOUDEN
# =========================================================

best_coords <- coords(
  roc_obj,
  x = "best",
  best.method = "youden",
  ret = c(
    "threshold",
    "sensitivity",
    "specificity"
  )
)

print(best_coords)

youden_df <- data.frame(
  FPR = 1 - as.numeric(best_coords["specificity"]),
  TPR = as.numeric(best_coords["sensitivity"])
)

print(youden_df)

# =========================================================
# 18. CURVA ROC FINAL
# =========================================================

p_roc <- ggplot(
  roc_df,
  aes(
    x = FPR,
    y = TPR
  )
) +
  geom_line(
    linewidth = 0.9,
    color = "#2C7FB8",
    lineend = "round"
  ) +
  geom_abline(
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey60"
  ) +
  geom_point(
    data = youden_df,
    mapping = aes(
      x = FPR,
      y = TPR
    ),
    inherit.aes = FALSE,
    size = 2.8,
    color = "red"
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE
  ) +
  labs(
    x = "False Positive Rate",
    y = "True Positive Rate"
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
# 19. HEATMAP
# =========================================================

annotation_col <- data.frame(
  Group = group
)

rownames(annotation_col) <- colnames(expr_scaled)

pheatmap(
  expr_scaled,
  annotation_col = annotation_col,
  show_colnames = FALSE,
  fontsize_row = 10,
  clustering_method = "complete",
  main = ""
)

# =========================================================
# 20. EXPORTAR RESULTADOS
# =========================================================

write.table(
  result_df,
  file = "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer/GSE173094/oncogenic_axis_score.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# =========================================================
# 21. SALVAR FIGURAS
# =========================================================

ggsave(
  filename = "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer/GSE173094/ROC_oncogenic_axis.png",
  plot = p_roc,
  width = 5,
  height = 5,
  dpi = 600
)

ggsave(
  filename = "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer/GSE173094/Boxplot_oncogenic_axis.png",
  plot = p_box,
  width = 4,
  height = 5,
  dpi = 600
)

# =========================================================
# 22. FINAL
# =========================================================

cat("\nAnálise concluída com sucesso.\n")
cat("AUC =", round(as.numeric(auc(roc_obj)), 3), "\n")