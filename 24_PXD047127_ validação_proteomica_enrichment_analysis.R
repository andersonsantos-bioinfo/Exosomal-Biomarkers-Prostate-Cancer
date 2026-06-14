# =========================================================
# ENRIQUECIMENTO FUNCIONAL (GO + REACTOME)
# VERSÃO CORRIGIDA E ROBUSTA
# PXD047127 - Vesículas Extracelulares
# =========================================================

rm(list = ls())

# =========================================================
# PACOTES
# =========================================================

packages <- c(
  "data.table",
  "clusterProfiler",
  "ReactomePA",
  "org.Hs.eg.db",
  "enrichplot",
  "ggplot2",
  "dplyr",
  "stringr"
)

installed <- packages %in% installed.packages()

if(any(!installed)){
  install.packages(packages[!installed])
}

# Bioconductor
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_pkgs <- c(
  "clusterProfiler",
  "ReactomePA",
  "org.Hs.eg.db",
  "enrichplot"
)

for(p in bioc_pkgs){
  if(!requireNamespace(p, quietly = TRUE)){
    BiocManager::install(p)
  }
}

# carregar
lapply(packages, library, character.only = TRUE)

# =========================================================
# CAMINHO
# =========================================================

path <- "C:/Users/Administrator/Exossomos_prostata/Fase-Final/Validação Proteômica"

# =========================================================
# LEITURA DA MATRIZ PROTEICA
# =========================================================

pg <- fread(file.path(path, "report.pg_matrix.tsv"))

# conferir
dim(pg)

# =========================================================
# EXTRAIR GENES
# =========================================================

genes_raw <- pg$Genes

# remover NAs
genes_raw <- genes_raw[!is.na(genes_raw)]

# separar genes compostos por ;
genes_split <- strsplit(genes_raw, ";")

# vetor único
all_genes <- unique(unlist(genes_split))

# limpar espaços
all_genes <- trimws(all_genes)

# remover vazios
all_genes <- all_genes[all_genes != ""]

# conferir
cat("Total de genes detectados:", length(all_genes), "\n")

head(all_genes)

# =========================================================
# CONVERSÃO SYMBOL -> ENTREZ
# =========================================================

gene_df <- bitr(
  all_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# genes convertidos
entrez_ids <- unique(gene_df$ENTREZID)

cat("Genes convertidos para ENTREZ:", length(entrez_ids), "\n")

# =========================================================
# GO ENRICHMENT
# =========================================================

ego <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# conferir
ego_df <- as.data.frame(ego)

cat("GO terms enriquecidos:", nrow(ego_df), "\n")

# salvar tabela
write.csv(
  ego_df,
  file.path(path, "GO_Enrichment_Results.csv"),
  row.names = FALSE
)

# =========================================================
# REACTOME ENRICHMENT
# =========================================================

reactome <- enrichPathway(
  gene          = entrez_ids,
  organism      = "human",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# dataframe
reactome_df <- as.data.frame(reactome)

cat("Reactome pathways enriquecidas:", nrow(reactome_df), "\n")

# salvar tabela
write.csv(
  reactome_df,
  file.path(path, "Reactome_Enrichment_Results.csv"),
  row.names = FALSE
)

# =========================================================
# TOP RESULTADOS
# =========================================================

cat("\n====================================\n")
cat("TOP GO TERMS\n")
cat("====================================\n")

print(
  ego_df[, c("Description", "p.adjust")] %>%
    head(10)
)

cat("\n====================================\n")
cat("TOP REACTOME PATHWAYS\n")
cat("====================================\n")

print(
  reactome_df[, c("Description", "p.adjust")] %>%
    head(10)
)

# =========================================================
# VISUALIZAÇÃO GO
# =========================================================

pdf(
  file.path(path, "GO_dotplot.pdf"),
  width = 10,
  height = 8
)

print(
  dotplot(
    ego,
    showCategory = 15,
    font.size = 10
  ) +
    ggtitle("GO Biological Process Enrichment")
)

dev.off()

# =========================================================
# VISUALIZAÇÃO REACTOME
# =========================================================

pdf(
  file.path(path, "Reactome_dotplot.pdf"),
  width = 10,
  height = 8
)

print(
  dotplot(
    reactome,
    showCategory = 15,
    font.size = 10
  ) +
    ggtitle("Reactome Pathway Enrichment")
)

dev.off()

# =========================================================
# GENES DE INTERESSE (ASSINATURA FUNCIONAL)
# =========================================================

signature_genes <- c(
  "FN1",
  "ITGB1",
  "ACTB",
  "VIM",
  "ACTN1",
  "ACTN4",
  "CDC42",
  "RHOA"
)

found_signature <- intersect(signature_genes, all_genes)

cat("\n====================================\n")
cat("ASSINATURA FUNCIONAL DETECTADA\n")
cat("====================================\n")

print(found_signature)

# =========================================================
# SALVAR AMBIENTE
# =========================================================

save.image(
  file = file.path(
    path,
    paste0("Enrichment_Proteomics_", Sys.Date(), ".RData")
  )
)

cat("\n====================================\n")
cat("ANÁLISE CONCLUÍDA COM SUCESSO\n")
cat("====================================\n")