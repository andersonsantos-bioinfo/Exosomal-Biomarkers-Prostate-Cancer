# =========================================================
# 39.CORRELAÇÃO SPEARMAN
# =========================================================

library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# =========================
# DADOS
# =========================
matriz <- as.matrix(matriz_painel)

# garantir orientação correta
if (nrow(matriz) < ncol(matriz)) {
  # OK (genes x samples)
} else {
  matriz <- t(matriz)
}

# =========================
# VERIFICAR EMAS
# =========================
if (!exists("EMAS")) stop("Objeto EMAS não encontrado")

# alinhar nomes
common_samples <- intersect(colnames(matriz), names(EMAS))

cat("Samples em comum:", length(common_samples), "\n")

matriz <- matriz[, common_samples, drop = FALSE]
emas_vec <- EMAS[common_samples]

# =========================
# GENES
# =========================
genes_sel <- c(
  "SCHLAP1",
  "LINC01475",
  "HOXC6",
  "ETV1",
  "ZIC2",
  "SLC45A2"
)

genes_ok <- intersect(genes_sel, rownames(matriz))

cat("Genes encontrados:", genes_ok, "\n")

expr_sel <- matriz[genes_ok, , drop = FALSE]

# =========================
# DATA FRAME
# =========================
df_corr <- as.data.frame(t(expr_sel))
df_corr$EMAS <- emas_vec

df_corr <- df_corr %>%
  mutate(across(everything(), as.numeric))

# =========================
# CORRELAÇÃO
# =========================
cor_mat <- cor(df_corr, method = "spearman", use = "pairwise.complete.obs")

print(cor_mat)

# =========================
# MOSTRAR NO R (FORÇAR VISUALIZAÇÃO)
# =========================
pheatmap(
  cor_mat,
  color = colorRampPalette(
    rev(brewer.pal(11, "RdBu"))
  )(100),
  display_numbers = TRUE,
  fontsize_number = 10,
  main = "Spearman correlation matrix"
)