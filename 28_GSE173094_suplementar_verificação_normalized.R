# =========================================================
# 28.GSE173094_suplementar_¬verificação_normalized_clean.
# =========================================================

path <- "C:/Users/Administrator/Exossomos_prostata/Fase-Final/Validação_suplementar/GSE173094/GSE173094_normalized_clean.txt"

df <- read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("Dimensão do arquivo:\n")
print(dim(df))


# =========================================================
# 2. VERIFICAR ESTRUTURA GERAL
# =========================================================

cat("\nColunas:\n")
print(colnames(df))

cat("\nTipos de dados:\n")
print(sapply(df, class))


# =========================================================
# 3. VERIFICAR MIRNAs E AMOSTRAS
# =========================================================

cat("\nNúmero de miRNAs:\n")
print(nrow(df))

cat("\nNúmero de amostras:\n")
print(ncol(df))


# =========================================================
# 4. CHECAR IDENTIFICADORES
# =========================================================

cat("\nExemplos de miRNAs:\n")
print(head(df$ID_REF, 10))


# =========================================================
# 5. VERIFICAR VALORES NUMÉRICOS (UE SAMPLES)
# =========================================================

expr <- df[, grep("^UE", colnames(df))]

cat("\nResumo dos valores:\n")
print(summary(as.matrix(expr)))


# =========================================================
# 6. VERIFICAR VALORES PROBLEMÁTICOS
# =========================================================

cat("\nQuantidade de NA:\n")
print(sum(is.na(expr)))

cat("\nQuantidade de valores iguais a 40 (provável limite de Ct):\n")
print(sum(expr == 40, na.rm = TRUE))


# =========================================================
# 7. VERIFICAR VARIAÇÃO ENTRE AMOSTRAS
# =========================================================

cat("\nDesvio padrão por amostra:\n")
print(summary(apply(expr, 2, sd, na.rm = TRUE)))


# =========================================================
# 8. DETECTAR POSSÍVEL PROBLEMA DE NORMALIZAÇÃO
# =========================================================

cat("\nAmplitude geral (min/max):\n")
print(range(expr, na.rm = TRUE))


# =========================================================
# 9. CHECAR DISTRIBUIÇÃO GLOBAL
# =========================================================

hist(as.vector(as.matrix(expr)),
     main = "Distribuição dos valores (UE samples)",
     xlab = "Expressão",
     breaks = 50)