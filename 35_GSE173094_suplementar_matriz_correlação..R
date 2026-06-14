# =========================================================
# GSE173094 - MATRIZ DE CORRELAÇÃO
# miRNAs × eixo oncogênico
# =========================================================

library(data.table)
library(Hmisc)
library(pheatmap)
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
# =========================================================

expr <- as.matrix(
  df[, grep("^UE", colnames(df))]
)

rownames(expr) <- df$ID_REF

# =========================================================
# 3. miRNAs DO EIXO ONCOGÊNICO
# =========================================================

selected_miRNAs <- c(
  "hsa-miR-518e-002395",
  "hsa-miR-548d-5p-002237",
  "hsa-miR-142-3p-000464",
  "hsa-miR-23b-000400",
  "hsa-miR-548b-5p-002408",
  "hsa-let-7g-002282"
)

# manter apenas miRNAs presentes
selected_miRNAs <- intersect(
  selected_miRNAs,
  rownames(expr)
)

print(selected_miRNAs)

# =========================================================
# 4. PADRONIZAÇÃO (Z-SCORE)
# =========================================================

expr_z <- t(
  scale(
    t(expr)
  )
)

# =========================================================
# 5. SCORE ONCOGÊNICO
# =========================================================

oncogenic_score <- colMeans(
  expr_z[selected_miRNAs, ],
  na.rm = TRUE
)

head(oncogenic_score)

# =========================================================
# 6. ALINHAR AMOSTRAS
# =========================================================

samples_miRNA <- colnames(expr)
samples_score <- names(oncogenic_score)

common_samples <- intersect(
  samples_miRNA,
  samples_score
)

cat(
  "Amostras em comum:",
  length(common_samples),
  "\n"
)

# alinhar
expr_sub <- expr_z[, common_samples]

score_sub <- oncogenic_score[
  common_samples
]

# =========================================================
# 7. TRANSPOSTA
# amostras × miRNAs
# =========================================================

miRNA_t <- t(expr_sub)

dim(miRNA_t)

# =========================================================
# 8. MATRIZ DE CORRELAÇÃO (SPEARMAN)
# =========================================================

res <- rcorr(
  as.matrix(
    cbind(miRNA_t, score_sub)
  ),
  type = "spearman"
)

# =========================================================
# 9. EXTRAIR CORRELAÇÕES
# =========================================================

cor_vec <- res$r[, ncol(res$r)]

p_vec <- res$P[, ncol(res$P)]

# remover o próprio score
cor_vec <- cor_vec[-length(cor_vec)]
p_vec <- p_vec[-length(p_vec)]

# =========================================================
# 10. DATAFRAME FINAL
# =========================================================

cor_df <- data.frame(
  miRNA = names(cor_vec),
  correlation = cor_vec,
  p_value = p_vec
)

# remover NA
cor_df <- cor_df[
  complete.cases(cor_df),
]

# ordenar
cor_df <- cor_df[
  order(
    -abs(cor_df$correlation)
  ),
]

# =========================================================
# 11. RESULTADOS PRINCIPAIS
# =========================================================

head(cor_df, 20)

# =========================================================
# 12. miRNAs SIGNIFICATIVAS
# =========================================================

cor_sig <- cor_df[
  abs(cor_df$correlation) > 0.4 &
    cor_df$p_value < 0.05,
]

print(cor_sig)

# =========================================================
# 13. EXPORTAR RESULTADOS
# =========================================================

write.table(
  cor_df,
  file = "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer/GSE173094/correlation_miRNA_oncogenic_axis.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  cor_sig,
  file = "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer/GSE173094/correlation_significant_miRNAs.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# =========================================================
# 14. HEATMAP - TOP 20 miRNAs
# =========================================================

top_miRNAs <- cor_df$miRNA[1:20]

heat_data <- cor_df[
  match(top_miRNAs, cor_df$miRNA),
]

rownames(heat_data) <- heat_data$miRNA

heat_matrix <- as.matrix(
  heat_data["correlation"]
)

# =========================================================
# 15. PHEATMAP
# =========================================================

pheatmap(
  heat_matrix,
  color = colorRampPalette(
    c("blue", "white", "red")
  )(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 9,
  fontsize_col = 10,
  border_color = NA,
  main = "miRNAs correlated with oncogenic axis"
)

# =========================================================
# 16. BARPLOT - TOP CORRELAÇÕES
# =========================================================

top_plot <- head(cor_df, 15)

p_bar <- ggplot(
  top_plot,
  aes(
    x = reorder(miRNA, correlation),
    y = correlation,
    fill = correlation
  )
) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  labs(
    x = "",
    y = "Spearman correlation"
  ) +
  theme_classic(
    base_family = "Arial",
    base_size = 10
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(
      color = "black"
    ),
    axis.title = element_text(
      color = "black"
    )
  )

print(p_bar)

# =========================================================
# 17. SALVAR FIGURA
# =========================================================

ggsave(
  filename = "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer/GSE173094/Top_miRNA_correlations.png",
  plot = p_bar,
  width = 5,
  height = 6,
  dpi = 600
)

# =========================================================
# 18. FINAL
# =========================================================

cat("\n=====================================\n")
cat("Análise de correlação concluída.\n")
cat("miRNAs significativas:", nrow(cor_sig), "\n")
cat("=====================================\n")