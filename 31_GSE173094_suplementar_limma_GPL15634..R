# =========================================================
# 31.GSE173094 - LIMMA + ANOTAÇÃO CORRETA GPL15634
# =========================================================

library(GEOquery)
library(limma)
library(dplyr)
library(readr)

# =========================================================
# 1. CAMINHOS
# =========================================================

base_dir <- "C:/Users/Administrator/Exossomos_prostata/Fase-Final/Validação_suplementar/GSE173094"

expr_file <- file.path(base_dir,
                       "GSE173094_expr_matrix_clean.txt")

# =========================================================
# 2. LER MATRIZ
# =========================================================

df <- read.delim(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

# primeira coluna = IDs numéricos da plataforma
probe_ids <- df[[1]]

# remover primeira coluna da expressão
expr <- df[, -1]

# transformar em matriz
expr <- as.matrix(expr)

# converter para numeric
mode(expr) <- "numeric"

# =========================================================
# 3. OBTER ANOTAÇÃO GPL15634
# =========================================================

gpl <- getGEO("GPL15634", AnnotGPL = TRUE)

annot <- Table(gpl)

# selecionar colunas úteis
annot2 <- annot[, c("ID", "miRNA_ID")]

# remover vazios
annot2 <- annot2[annot2$miRNA_ID != "", ]

# =========================================================
# 4. MAPEAR IDs NUMÉRICOS ??? miRNAs
# =========================================================

# IDs da matriz são 1:381
# IDs da GPL seguem a mesma ordem

annot2$index <- 1:nrow(annot2)

# alinhar
miRNA_names <- annot2$miRNA_ID[
  match(probe_ids, annot2$index)
]

# fallback
miRNA_names[is.na(miRNA_names)] <-
  paste0("miRNA_", probe_ids[is.na(miRNA_names)])

# aplicar rownames
rownames(expr) <- miRNA_names

# =========================================================
# 5. REMOVER DUPLICATAS
# =========================================================

expr <- expr[!duplicated(rownames(expr)), ]

# =========================================================
# 6. DEFINIR GRUPOS
# =========================================================

group <- c(
  rep("Localized", 19),
  rep("Metastatic", 23)
)

group <- factor(
  group,
  levels = c("Localized", "Metastatic")
)

# =========================================================
# 7. DESIGN MATRIX
# =========================================================

design <- model.matrix(~ group)

# =========================================================
# 8. LIMMA
# =========================================================

fit <- lmFit(expr, design)
fit <- eBayes(fit)

# =========================================================
# 9. RESULTADOS
# =========================================================

res <- topTable(
  fit,
  coef = "groupMetastatic",
  number = Inf,
  adjust.method = "BH"
)

# adicionar miRNA explicitamente
res$miRNA <- rownames(res)

# =========================================================
# 10. SCORE MULTICRITÉRIO
# =========================================================

res$abs_logFC <- abs(res$logFC)

res$score <- (
  scale(res$abs_logFC) +
    scale(-log10(res$P.Value)) +
    scale(res$AveExpr)
)

# ordenar
res <- res[order(res$score, decreasing = TRUE), ]

# =========================================================
# 11. TOP BIOMARCADORES
# =========================================================

top_biomarkers <- head(res, 20)

# reorganizar colunas
top_biomarkers <- top_biomarkers[, c(
  "miRNA",
  "logFC",
  "AveExpr",
  "P.Value",
  "adj.P.Val",
  "score"
)]

# =========================================================
# 12. VERIFICAR let-7g
# =========================================================

subset(res, grepl("let-7g", miRNA, ignore.case = TRUE))

# =========================================================
# 13. SALVAR
# =========================================================

write.csv(
  res,
  file.path(base_dir,
            "GSE173094_limma_resultados_completos.csv"),
  row.names = FALSE
)

write.csv(
  top_biomarkers,
  file.path(base_dir,
            "GSE173094_top_biomarcadores.csv"),
  row.names = FALSE
)

cat("\nAnálise concluída com sucesso.\n")