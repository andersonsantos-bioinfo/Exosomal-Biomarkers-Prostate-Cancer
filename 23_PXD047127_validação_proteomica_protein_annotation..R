# =========================================================
# 23_PXD047127_protein_annotation.R
# Proteomic Validation Cohort - PXD047127
# Protein Annotation and Target Mapping
# Exosomal Biomarkers in Prostate Cancer
# =========================================================

# =========================================================
# PACOTES
# =========================================================
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
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

dir.create(results_dir, showWarnings = FALSE)
dir.create(figures_dir, showWarnings = FALSE)

# =========================================================
# CARREGAR MATRIZ PROCESSADA
# =========================================================
input_file <- file.path(
  processed_dir,
  "PXD047127_expression_log2.csv"
)

if (!file.exists(input_file)) {
  
  stop(
    paste(
      "\nArquivo não encontrado:\n",
      input_file
    )
  )
}

cat("\n=====================================\n")
cat("IMPORTANDO MATRIZ PROCESSADA\n")
cat("=====================================\n")

expr <- read.csv(
  input_file,
  row.names = 1,
  check.names = FALSE
)

# =========================================================
# DIMENSÕES
# =========================================================
cat("\n=====================================\n")
cat("DIMENSÕES DA MATRIZ\n")
cat("=====================================\n")

print(dim(expr))

# =========================================================
# GENES DO PAINEL MULTIGÊNICO
# =========================================================
targets <- c(
  "ETV1",
  "HOXC6",
  "SLC45A2",
  "ZIC2",
  "SCHLAP1",
  "LINC01475"
)

cat("\n=====================================\n")
cat("GENES ALVO\n")
cat("=====================================\n")

print(targets)

# =========================================================
# IDENTIFICAR GENES ENCONTRADOS
# =========================================================
all_genes <- rownames(expr)

found_genes <- targets[
  targets %in% all_genes
]

missing_genes <- targets[
  !targets %in% all_genes
]

cat("\n=====================================\n")
cat("GENES ENCONTRADOS\n")
cat("=====================================\n")

print(found_genes)

cat("\n=====================================\n")
cat("GENES AUSENTES\n")
cat("=====================================\n")

print(missing_genes)

# =========================================================
# EXPORTAR GENES ENCONTRADOS
# =========================================================
write.csv(
  data.frame(
    Found_Genes = found_genes
  ),
  file.path(
    results_dir,
    "PXD047127_found_genes.csv"
  ),
  row.names = FALSE
)

# =========================================================
# EXPORTAR GENES AUSENTES
# =========================================================
write.csv(
  data.frame(
    Missing_Genes = missing_genes
  ),
  file.path(
    results_dir,
    "PXD047127_missing_genes.csv"
  ),
  row.names = FALSE
)

# =========================================================
# BUSCA PARCIAL (FAMÍLIAS)
# =========================================================
cat("\n=====================================\n")
cat("BUSCA PARCIAL DE FAMÍLIAS\n")
cat("=====================================\n")

family_ETV <- grep(
  "ETV",
  all_genes,
  value = TRUE
)

family_ZIC <- grep(
  "ZIC",
  all_genes,
  value = TRUE
)

family_SLC45 <- grep(
  "SLC45",
  all_genes,
  value = TRUE
)

cat("\nGenes ETV encontrados:\n")
print(family_ETV)

cat("\nGenes ZIC encontrados:\n")
print(family_ZIC)

cat("\nGenes SLC45 encontrados:\n")
print(family_SLC45)

# =========================================================
# EXPORTAR BUSCAS PARCIAIS
# =========================================================
write.csv(
  data.frame(ETV_Family = family_ETV),
  file.path(
    results_dir,
    "PXD047127_ETV_family.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(ZIC_Family = family_ZIC),
  file.path(
    results_dir,
    "PXD047127_ZIC_family.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(SLC45_Family = family_SLC45),
  file.path(
    results_dir,
    "PXD047127_SLC45_family.csv"
  ),
  row.names = FALSE
)

# =========================================================
# CONVERSÃO PARA ENTREZ IDs
# =========================================================
cat("\n=====================================\n")
cat("CONVERSÃO ENTREZ IDs\n")
cat("=====================================\n")

gene_df <- bitr(
  all_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

cat("\nGenes convertidos:\n")
print(head(gene_df))

# =========================================================
# EXPORTAR ENTREZ IDs
# =========================================================
write.csv(
  gene_df,
  file.path(
    results_dir,
    "PXD047127_Entrez_IDs.csv"
  ),
  row.names = FALSE
)

# =========================================================
# LISTA FINAL ENTREZ
# =========================================================
entrez_ids <- unique(
  gene_df$ENTREZID
)

# =========================================================
# RESUMO FINAL
# =========================================================
summary_df <- data.frame(
  Metric = c(
    "Total proteins",
    "Detected target genes",
    "Missing target genes",
    "Genes converted to Entrez"
  ),
  
  Value = c(
    length(all_genes),
    length(found_genes),
    length(missing_genes),
    length(entrez_ids)
  )
)

write.csv(
  summary_df,
  file.path(
    results_dir,
    "PXD047127_annotation_summary.csv"
  ),
  row.names = FALSE
)

# =========================================================
# TABELA CONSOLIDADA DOS TARGETS
# =========================================================
target_table <- data.frame(
  Target = targets,
  Detected = targets %in% found_genes
)

write.csv(
  target_table,
  file.path(
    results_dir,
    "PXD047127_target_detection_table.csv"
  ),
  row.names = FALSE
)

# =========================================================
# SESSION INFO
# =========================================================
writeLines(
  capture.output(sessionInfo()),
  file.path(
    results_dir,
    "sessionInfo_PXD047127_annotation.txt"
  )
)

# =========================================================
# FINALIZAÇÃO
# =========================================================
cat("\n=====================================\n")
cat("ANNOTATION CONCLUÍDA\n")
cat("=====================================\n")

cat("\nGenes detectados:\n")
print(found_genes)

cat("\nArquivos exportados para:\n")
cat(results_dir, "\n")