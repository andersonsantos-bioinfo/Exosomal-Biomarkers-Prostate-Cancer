# =========================================================
# GSE173094 - PIPELINE COMPLETO LIMMA + RANQUEAMENTO
# =========================================================

rm(list = ls())

# =========================================================
# PACOTES
# =========================================================

library(limma)
library(dplyr)
library(stringr)

# =========================================================
# CAMINHO BASE
# =========================================================

base_dir <- "C:/Users/Administrator/Exossomos_prostata/Fase-Final/Validação_suplementar/GSE173094"

# =========================================================
# ARQUIVOS
# =========================================================

expr_file <- file.path(
  base_dir,
  "GSE173094_expr_matrix_clean.txt"
)

meta_file <- file.path(
  base_dir,
  "GSE173094_series_matrix.txt"
)

# =========================================================
# IMPORTAR MATRIZ DE EXPRESSÃO
# =========================================================

df <- read.table(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Primeira coluna = genes
rownames(df) <- df[,1]

# Remover primeira coluna
df <- df[,-1]

# Converter para numérico
df[] <- lapply(df, as.numeric)

# Converter em matriz
df <- as.matrix(df)

# =========================================================
# IMPORTAR METADATA
# =========================================================

meta <- readLines(meta_file)

# =========================================================
# EXTRAIR LINHAS DE CARACTERÍSTICAS
# =========================================================

group_lines <- grep(
  "characteristics",
  meta,
  value = TRUE,
  ignore.case = TRUE
)

# =========================================================
# IDENTIFICAR LINHA DE STATUS
# =========================================================

status_line <- grep(
  "localized|metast",
  group_lines,
  value = TRUE,
  ignore.case = TRUE
)[1]

cat("\nLinha identificada:\n")
print(status_line)

# =========================================================
# LIMPAR TEXTO
# =========================================================

status_line <- sub(
  "!Sample_characteristics_ch1\t",
  "",
  status_line
)

vals <- unlist(strsplit(status_line, "\t"))

vals <- gsub("\"", "", vals)

# =========================================================
# CRIAR VETOR DE GRUPOS
# =========================================================

group <- ifelse(
  grepl("localized|local", vals, ignore.case = TRUE),
  "Localized",
  "Metastatic"
)

# Nomear amostras
names(group) <- paste0("UE", 1:length(group))

# =========================================================
# VALIDAR GRUPOS
# =========================================================

cat("\nDistribuição dos grupos:\n")
print(table(group))

# =========================================================
# EXTRAIR MATRIZ EXPRESSÃO
# =========================================================

expr <- df[, grep("^UE", colnames(df))]

# =========================================================
# ALINHAR MATRIZ
# =========================================================

expr <- expr[, names(group)]

# =========================================================
# CHECAGEM
# =========================================================

cat("\nAlinhamento correto:\n")
print(all(colnames(expr) == names(group)))

cat("\nDimensão expressão:\n")
print(dim(expr))

# =========================================================
# REMOVER GENES INVÁLIDOS
# =========================================================

expr <- expr[rowSums(is.na(expr)) == 0, ]

expr <- expr[rowSums(expr) > 0, ]

cat("\nDimensão após limpeza:\n")
print(dim(expr))

# =========================================================
# FATOR ORDENADO
# =========================================================

group <- factor(
  group,
  levels = c("Localized", "Metastatic")
)

# =========================================================
# DESIGN MATRIX
# =========================================================

design <- model.matrix(~ group)

cat("\nDesign matrix:\n")
print(head(design))

# =========================================================
# ANÁLISE DIFERENCIAL - LIMMA
# =========================================================

fit <- lmFit(expr, design)

fit <- eBayes(fit)

# =========================================================
# RESULTADOS
# =========================================================

results <- topTable(
  fit,
  coef = "groupMetastatic",
  number = Inf,
  adjust.method = "BH"
)

# =========================================================
# ADICIONAR MÉTRICAS
# =========================================================

res <- results

res$abs_logFC <- abs(res$logFC)

# =========================================================
# SCORE MULTICRITÉRIO
# =========================================================

res$score <- (
  scale(res$abs_logFC) * 0.4 +
    scale(-log10(res$P.Value)) * 0.4 +
    scale(res$AveExpr) * 0.2
)

# Converter para vetor numérico
res$score <- as.numeric(res$score)

# =========================================================
# RANQUEAMENTO
# =========================================================

res_ranked <- res[
  order(res$score, decreasing = TRUE),
]

# =========================================================
# TOP 20
# =========================================================

cat("\nTOP 20 genes ranqueados:\n")

print(head(res_ranked, 20))

# =========================================================
# FILTRO DE CANDIDATOS
# =========================================================

candidates <- res[
  res$P.Value < 0.01 &
    abs(res$logFC) > 1 &
    res$AveExpr > quantile(res$AveExpr, 0.25),
]

# =========================================================
# UP / DOWN REGULATED
# =========================================================

upregulated <- candidates[
  candidates$logFC > 0,
]

downregulated <- candidates[
  candidates$logFC < 0,
]

# =========================================================
# TOP BIOMARCADORES
# =========================================================

top_biomarkers <- head(res_ranked, 15)

cat("\nTOP BIOMARCADORES:\n")

print(top_biomarkers)

# =========================================================
# SALVAR RESULTADOS
# =========================================================

write.table(
  results,
  file.path(base_dir, "GSE173094_limma_results.txt"),
  sep = "\t",
  quote = FALSE
)

write.table(
  res_ranked,
  file.path(base_dir, "GSE173094_ranked_genes.txt"),
  sep = "\t",
  quote = FALSE
)

write.table(
  candidates,
  file.path(base_dir, "GSE173094_candidate_biomarkers.txt"),
  sep = "\t",
  quote = FALSE
)

write.table(
  upregulated,
  file.path(base_dir, "GSE173094_upregulated.txt"),
  sep = "\t",
  quote = FALSE
)

write.table(
  downregulated,
  file.path(base_dir, "GSE173094_downregulated.txt"),
  sep = "\t",
  quote = FALSE
)

write.table(
  top_biomarkers,
  file.path(base_dir, "GSE173094_top15_biomarkers.txt"),
  sep = "\t",
  quote = FALSE
)

# =========================================================
# SALVAR AMBIENTE R
# =========================================================

today <- Sys.Date()

env_file <- file.path(
  base_dir,
  paste0("R_env_", today, ".RData")
)

save.image(env_file)

cat("\nAmbiente salvo em:\n")
print(env_file)

# =========================================================
# RESUMO FINAL
# =========================================================

cat("\n=============================\n")
cat("ANÁLISE CONCLUÍDA\n")
cat("=============================\n")

cat("\nGenes totais:\n")
print(nrow(results))

cat("\nCandidatos:\n")
print(nrow(candidates))

cat("\nUpregulated:\n")
print(nrow(upregulated))

cat("\nDownregulated:\n")
print(nrow(downregulated))