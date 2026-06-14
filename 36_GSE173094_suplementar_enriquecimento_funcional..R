# =========================================================
# 36_GSE173094_suplementar_enriquecimento_funcional
# Interseção entre genes alvo dos miRNAs
# + enriquecimento funcional
# + GO
# + KEGG
# + network plot compatível
# =========================================================

# =========================================================
# 1. PACOTES
# =========================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages <- c(
  "multiMiR",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "ggplot2"
)

for(pkg in packages){
  
  if(!requireNamespace(pkg, quietly = TRUE)){
    BiocManager::install(pkg)
  }
  
  library(pkg, character.only = TRUE)
}

# =========================================================
# 2. VERIFICAR OBJETO cor_sig
# =========================================================

if(!exists("cor_sig")){
  stop("Objeto 'cor_sig' não encontrado no ambiente.")
}

# =========================================================
# 3. EXTRAIR miRNAs SIGNIFICATIVOS
# =========================================================

miRNAs <- unique(cor_sig$miRNA)

cat("\nNúmero de miRNAs:\n")
print(length(miRNAs))

head(miRNAs)

# =========================================================
# 4. LIMPAR NOMES DOS miRNAs
# =========================================================

# remover IDs finais da plataforma
miRNAs_clean <- sub("-\\d+$", "", miRNAs)

# garantir nomenclatura madura
miRNAs_clean <- ifelse(
  grepl("-3p|-5p$", miRNAs_clean),
  miRNAs_clean,
  paste0(miRNAs_clean, "-5p")
)

# remover duplicados
miRNAs_clean <- unique(miRNAs_clean)

cat("\nmiRNAs após limpeza:\n")
print(miRNAs_clean)

# =========================================================
# 5. BUSCAR ALVOS VALIDADOS
# =========================================================

targets <- get_multimir(
  mirna = miRNAs_clean,
  table = "validated"
)

targets_df <- targets@data

cat("\nDimensão da tabela de alvos:\n")
print(dim(targets_df))

# =========================================================
# 6. EXTRAIR GENES ALVO
# =========================================================

target_genes <- unique(
  targets_df$target_symbol
)

target_genes <- na.omit(target_genes)

cat("\nNúmero total de genes alvo:\n")
print(length(target_genes))

head(target_genes)

# =========================================================
# 7. DEFINIR EIXO ONCOGÊNICO
# =========================================================

genes_eixo <- c(
  "HOXC6",
  "ETV1",
  "ZIC2",
  "SLC45A2"
)

cat("\nGenes do eixo oncogênico:\n")
print(genes_eixo)

# =========================================================
# 8. INTERSEÇÃO
# =========================================================

intersect_genes <- intersect(
  target_genes,
  genes_eixo
)

cat("\nInterseção encontrada:\n")
print(intersect_genes)

# =========================================================
# 9. TESTE HIPERGEOMÉTRICO
# =========================================================

universo <- keys(
  org.Hs.eg.db,
  keytype = "SYMBOL"
)

# número de interseções
k <- length(intersect_genes)

# genes alvo
m <- length(target_genes)

# genes eixo
n <- length(genes_eixo)

# universo total
N <- length(universo)

# teste hipergeométrico
p_value <- phyper(
  q = k - 1,
  m = m,
  n = N - m,
  k = n,
  lower.tail = FALSE
)

cat("\n=====================================\n")
cat("TESTE HIPERGEOMÉTRICO\n")
cat("=====================================\n")

cat("Interseções:", k, "\n")
cat("Genes alvo:", m, "\n")
cat("Genes eixo:", n, "\n")
cat("Universo:", N, "\n")
cat("P-value:", p_value, "\n")

# =========================================================
# 10. CONVERTER PARA ENTREZ
# =========================================================

genes_entrez <- bitr(
  target_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

genes_entrez <- genes_entrez[
  !duplicated(genes_entrez$ENTREZID),
]

cat("\nGenes convertidos para ENTREZ:\n")
print(head(genes_entrez))

# =========================================================
# 11. ENRIQUECIMENTO KEGG
# =========================================================

kegg <- enrichKEGG(
  gene = genes_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

cat("\nTop KEGG pathways:\n")
print(head(kegg@result))

# =========================================================
# 12. ENRIQUECIMENTO GO
# =========================================================

go <- enrichGO(
  gene = genes_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

cat("\nTop GO Biological Processes:\n")
print(head(go@result))

# =========================================================
# 13. PROCESSOS BIOLÓGICOS RELEVANTES
# =========================================================

genes_interessantes <- subset(
  go@result,
  grepl(
    "transcription|differentiation|development|stem|migration|epithelial|plasticity",
    Description,
    ignore.case = TRUE
  )
)

cat("\nProcessos relacionados ao programa oncogênico:\n")
print(head(genes_interessantes, 20))

# =========================================================
# 14. DOTPLOT GO
# =========================================================

p_go <- dotplot(
  go,
  showCategory = 15,
  font.size = 10
) +
  ggtitle("GO Biological Processes") +
  theme_bw(base_size = 12)

print(p_go)

# =========================================================
# 15. DOTPLOT KEGG
# =========================================================

p_kegg <- dotplot(
  kegg,
  showCategory = 15,
  font.size = 10
) +
  ggtitle("KEGG Pathways") +
  theme_bw(base_size = 12)

print(p_kegg)

# =========================================================
# 16. NETWORK PLOT COMPATÍVEL
# =========================================================

cat("\nGerando network plot...\n")

try({
  
  p_net <- cnetplot(
    go,
    showCategory = 8,
    node_label = "all"
  )
  
  print(p_net)
  
}, silent = TRUE)

# =========================================================
# 17. TABELA FINAL DE RESULTADOS
# =========================================================

final_table <- data.frame(
  Metric = c(
    "miRNAs significativos",
    "Genes alvo validados",
    "Genes eixo oncogênico",
    "Interseções observadas",
    "P-value hipergeométrico"
  ),
  
  Value = c(
    length(miRNAs_clean),
    length(target_genes),
    length(genes_eixo),
    k,
    signif(p_value, 3)
  )
)

print(final_table)

# =========================================================
# 18. EXPORTAR RESULTADOS
# =========================================================

write.table(
  target_genes,
  file = "GSE173094_validated_targets.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  intersect_genes,
  file = "GSE173094_oncogenic_intersection.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  go@result,
  file = "GSE173094_GO_results.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  kegg@result,
  file = "GSE173094_KEGG_results.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  genes_interessantes,
  file = "GSE173094_GO_oncogenic_programs.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  final_table,
  file = "GSE173094_summary_statistics.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================================================
# 19. SALVAR AMBIENTE
# =========================================================

save(
  miRNAs,
  miRNAs_clean,
  target_genes,
  intersect_genes,
  genes_interessantes,
  kegg,
  go,
  final_table,
  file = "GSE173094_oncogenic_program_analysis.RData"
)

# =========================================================
# 20. FINAL
# =========================================================

cat("\n========================================\n")
cat("ANÁLISE FINALIZADA COM SUCESSO\n")
cat("========================================\n")
cat("Arquivos exportados:\n")
cat("- GSE173094_validated_targets.txt\n")
cat("- GSE173094_oncogenic_intersection.txt\n")
cat("- GSE173094_GO_results.txt\n")
cat("- GSE173094_KEGG_results.txt\n")
cat("- GSE173094_GO_oncogenic_programs.txt\n")
cat("- GSE173094_summary_statistics.txt\n")
cat("- GSE173094_oncogenic_program_analysis.RData\n")
cat("========================================\n")