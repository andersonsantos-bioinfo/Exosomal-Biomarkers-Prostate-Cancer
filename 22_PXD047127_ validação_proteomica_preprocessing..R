# =========================================================
# 22_PXD047127_validacao_proteomica_preprocessing.R
# Proteomic Validation Cohort - PXD047127
# Preprocessing Pipeline
# Exosomal Biomarkers in Prostate Cancer
# =========================================================

# =========================================================
# PACOTES
# =========================================================
library(tidyverse)
library(data.table)
library(janitor)

# =========================================================
# DIRETÓRIOS
# =========================================================
base_dir <- "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer"

proteomics_dir <- file.path(
  base_dir,
  "Validação_Proteômica"
)

processed_dir <- file.path(
  proteomics_dir,
  "processed"
)

results_dir <- file.path(
  proteomics_dir,
  "results"
)

figures_dir <- file.path(
  proteomics_dir,
  "figures"
)

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# ARQUIVOS DISPONÍVEIS
# =========================================================
cat("\n=====================================\n")
cat("ARQUIVOS ENCONTRADOS\n")
cat("=====================================\n")

print(list.files(proteomics_dir))

# =========================================================
# MATRIZ PRINCIPAL
# =========================================================
# Mais indicada para validação proteômica

input_file <- file.path(
  proteomics_dir,
  "report.unique_genes_matrix.tsv"
)

# =========================================================
# VERIFICAR ARQUIVO
# =========================================================
if (!file.exists(input_file)) {
  
  stop(
    paste(
      "\nArquivo não encontrado:\n",
      input_file
    )
  )
}

cat("\nArquivo selecionado:\n")
cat(input_file, "\n")

# =========================================================
# IMPORTAÇÃO
# =========================================================
cat("\n=====================================\n")
cat("IMPORTANDO MATRIZ PROTEÔMICA\n")
cat("=====================================\n")

expr_raw <- fread(
  input_file,
  sep = "\t",
  data.table = FALSE
)

# =========================================================
# LIMPEZA NOMES COLUNAS
# =========================================================
expr_raw <- expr_raw %>%
  clean_names()

# =========================================================
# DIMENSÕES INICIAIS
# =========================================================
cat("\n=====================================\n")
cat("DIMENSÕES INICIAIS\n")
cat("=====================================\n")

print(dim(expr_raw))

# =========================================================
# COLUNAS DISPONÍVEIS
# =========================================================
cat("\n=====================================\n")
cat("COLUNAS DISPONÍVEIS\n")
cat("=====================================\n")

print(colnames(expr_raw))

# =========================================================
# IDENTIFICAR COLUNA DE GENE
# =========================================================
possible_gene_cols <- c(
  "gene",
  "gene_symbol",
  "gene_name",
  "genes",
  "symbol"
)

gene_col <- intersect(
  possible_gene_cols,
  colnames(expr_raw)
)[1]

if (is.na(gene_col)) {
  
  gene_col <- colnames(expr_raw)[1]
  
  cat(
    "\nNenhuma coluna padrão encontrada.\n",
    "Usando primeira coluna como identificador.\n"
  )
}

cat("\nColuna identificadora:", gene_col, "\n")

# =========================================================
# DEFINIR ROW NAMES
# =========================================================
rownames(expr_raw) <- expr_raw[[gene_col]]

# =========================================================
# REMOVER COLUNA DE GENE
# =========================================================
expr <- expr_raw %>%
  dplyr::select(-all_of(gene_col))

# =========================================================
# CONVERTER PARA NUMÉRICO
# =========================================================
expr <- expr %>%
  mutate(
    across(
      everything(),
      as.numeric
    )
  )

expr <- as.data.frame(expr)

# =========================================================
# REMOVER PROTEÍNAS COM MUITOS NAs
# =========================================================
na_threshold <- 0.50

keep <- rowMeans(is.na(expr)) < na_threshold

expr <- expr[keep, ]

cat("\nProteínas mantidas:", nrow(expr), "\n")

# =========================================================
# IMPUTAÇÃO POR MEDIANA
# =========================================================
expr <- expr %>%
  mutate(
    across(
      everything(),
      ~ifelse(
        is.na(.),
        median(., na.rm = TRUE),
        .
      )
    )
  )

# =========================================================
# LOG2 TRANSFORMATION
# =========================================================
expr_log <- log2(expr + 1)

# =========================================================
# Z-SCORE NORMALIZATION
# =========================================================
expr_scaled <- t(
  scale(
    t(expr_log)
  )
)

expr_scaled <- as.data.frame(expr_scaled)

# =========================================================
# DIMENSÕES FINAIS
# =========================================================
cat("\n=====================================\n")
cat("DIMENSÕES FINAIS\n")
cat("=====================================\n")

print(dim(expr_scaled))

# =========================================================
# EXPORTAR MATRIZES
# =========================================================
write.csv(
  expr,
  file.path(
    processed_dir,
    "PXD047127_expression_processed.csv"
  )
)

write.csv(
  expr_log,
  file.path(
    processed_dir,
    "PXD047127_expression_log2.csv"
  )
)

write.csv(
  expr_scaled,
  file.path(
    processed_dir,
    "PXD047127_expression_scaled.csv"
  )
)

# =========================================================
# ESTATÍSTICAS GERAIS
# =========================================================
summary_df <- data.frame(
  Metric = c(
    "Proteins",
    "Samples"
  ),
  Value = c(
    nrow(expr_scaled),
    ncol(expr_scaled)
  )
)

write.csv(
  summary_df,
  file.path(
    results_dir,
    "PXD047127_preprocessing_summary.csv"
  ),
  row.names = FALSE
)

# =========================================================
# DISTRIBUIÇÃO DAS AMOSTRAS
# =========================================================
plot_df <- expr_log %>%
  pivot_longer(
    cols = everything(),
    names_to = "Sample",
    values_to = "Expression"
  )

# =========================================================
# QC BOXPLOT
# =========================================================
qc_plot <- ggplot(
  plot_df,
  aes(
    x = Sample,
    y = Expression
  )
) +
  
  geom_boxplot(
    fill = "#4A86B8",
    alpha = 0.85,
    outlier.size = 0.2
  ) +
  
  labs(
    title = "Proteomic Expression Distribution - PXD047127",
    x = "",
    y = "log2(Expression + 1)"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    
    plot.title = element_text(
      face = "bold",
      size = 14
    )
  )

# =========================================================
# MOSTRAR PLOT
# =========================================================
print(qc_plot)

# =========================================================
# SALVAR FIGURA
# =========================================================
ggsave(
  filename = file.path(
    figures_dir,
    "PXD047127_Boxplot_QC.png"
  ),
  plot = qc_plot,
  width = 11,
  height = 6,
  dpi = 600
)

# =========================================================
# SALVAR OBJETOS IMPORTANTES
# =========================================================
save(
  expr,
  expr_log,
  expr_scaled,
  expr_raw,
  file = file.path(
    processed_dir,
    "PXD047127_preprocessed_objects.RData"
  )
)

# =========================================================
# SESSION INFO
# =========================================================
writeLines(
  capture.output(sessionInfo()),
  file.path(
    results_dir,
    "sessionInfo_PXD047127.txt"
  )
)

# =========================================================
# FINALIZAÇÃO
# =========================================================
cat("\n=====================================\n")
cat("PREPROCESSAMENTO CONCLUÍDO\n")
cat("=====================================\n")

cat("\nArquivos exportados para:\n")
cat(processed_dir, "\n")

cat("\nFigura QC salva em:\n")
cat(figures_dir, "\n")