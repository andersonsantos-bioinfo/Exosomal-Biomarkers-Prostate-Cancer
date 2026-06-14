# =========================================================
# ATIVAÇÃO DO EIXO ONCOGÊNICO
# TCGA-PRAD
# SCHLAP1 / LINC01475 vs
# HOXC6 - ETV1 - ZIC2 - SLC45A2
# =========================================================

# =========================
# PACOTES
# =========================
library(tidyverse)
library(ggpubr)
library(rstatix)
library(pheatmap)
library(RColorBrewer)

# =========================
# MATRIZ DE EXPRESSÃO
# =========================
matriz_expr <- as.matrix(matriz_painel)

# verificar genes
genes_total <- c(
  "SCHLAP1",
  "LINC01475",
  "HOXC6",
  "ETV1",
  "ZIC2",
  "SLC45A2"
)

genes_total %in% rownames(matriz_expr)

# =========================
# DIRETÓRIO
# =========================
outdir <- "C:/Users/Administrator/Exossomos_prostata/Fase-Final/Validação_suplementar/GSE173094/Figuras/Eixo_Oncogenico"

dir.create(outdir,
           recursive = TRUE,
           showWarnings = FALSE)

# =========================
# EXTRAÇÃO DA EXPRESSÃO
# =========================
expr <- matriz_expr[
  rownames(matriz_expr) %in% genes_total,
]

expr <- as.data.frame(t(expr))

expr$Sample <- rownames(expr)

# =========================
# FUNÇÃO TERCIS
# =========================
criar_grupo_tercil <- function(x){
  
  q1 <- quantile(x, 0.33, na.rm = TRUE)
  q2 <- quantile(x, 0.67, na.rm = TRUE)
  
  case_when(
    x <= q1 ~ "Low",
    x >= q2 ~ "High",
    TRUE ~ NA_character_
  )
}

# =========================
# CRIAR GRUPOS
# =========================
expr$SCHLAP1_group <- criar_grupo_tercil(expr$SCHLAP1)

expr$LINC01475_group <- criar_grupo_tercil(expr$LINC01475)

# remover intermediários
dados_sch <- expr %>%
  filter(!is.na(SCHLAP1_group))

dados_linc <- expr %>%
  filter(!is.na(LINC01475_group))

# =========================
# GENES DO EIXO
# =========================
genes_eixo <- c(
  "HOXC6",
  "ETV1",
  "ZIC2",
  "SLC45A2"
)

# =========================================================
# FUNÇÃO DE ANÁLISE
# =========================================================
analisar_gene <- function(data,
                          grupo_var,
                          gene_alvo){
  
  form <- as.formula(
    paste0(gene_alvo, " ~ ", grupo_var)
  )
  
  teste <- wilcox_test(
    data,
    formula = form
  ) %>%
    add_significance()
  
  p <- ggplot(
    data,
    aes_string(
      x = grupo_var,
      y = gene_alvo,
      fill = grupo_var
    )
  ) +
    
    geom_boxplot(
      width = 0.65,
      alpha = 0.9,
      outlier.shape = NA
    ) +
    
    geom_jitter(
      width = 0.12,
      size = 1.4,
      alpha = 0.6
    ) +
    
    stat_compare_means(
      method = "wilcox.test",
      label = "p.format",
      size = 4
    ) +
    
    scale_fill_manual(
      values = c(
        "Low" = "#3C5488",
        "High" = "#B2182B"
      )
    ) +
    
    labs(
      x = "",
      y = paste0(gene_alvo, " expression")
    ) +
    
    theme_classic(base_family = "Arial") +
    
    theme(
      legend.position = "none",
      axis.title = element_text(
        face = "bold",
        size = 12
      ),
      axis.text = element_text(
        color = "black",
        size = 10
      )
    )
  
  return(list(
    plot = p,
    stats = teste
  ))
}

# =========================================================
# SCHLAP1-HIGH
# =========================================================

stats_sch <- list()

for(gene in genes_eixo){
  
  res <- analisar_gene(
    data = dados_sch,
    grupo_var = "SCHLAP1_group",
    gene_alvo = gene
  )
  
  stats_sch[[gene]] <- res$stats
  
  ggsave(
    filename = file.path(
      outdir,
      paste0(
        "SCHLAP1_",
        gene,
        "_boxplot.png"
      )
    ),
    plot = res$plot,
    width = 4,
    height = 4,
    dpi = 600
  )
}

# =========================================================
# LINC01475-HIGH
# =========================================================

stats_linc <- list()

for(gene in genes_eixo){
  
  res <- analisar_gene(
    data = dados_linc,
    grupo_var = "LINC01475_group",
    gene_alvo = gene
  )
  
  stats_linc[[gene]] <- res$stats
  
  ggsave(
    filename = file.path(
      outdir,
      paste0(
        "LINC01475_",
        gene,
        "_boxplot.png"
      )
    ),
    plot = res$plot,
    width = 4,
    height = 4,
    dpi = 600
  )
}

# =========================================================
# TABELAS ESTATÍSTICAS
# =========================================================

tabela_sch <- bind_rows(
  stats_sch,
  .id = "Gene"
)

tabela_linc <- bind_rows(
  stats_linc,
  .id = "Gene"
)

write.csv(
  tabela_sch,
  file.path(
    outdir,
    "Tabela_Wilcoxon_SCHLAP1.csv"
  ),
  row.names = FALSE
)

write.csv(
  tabela_linc,
  file.path(
    outdir,
    "Tabela_Wilcoxon_LINC01475.csv"
  ),
  row.names = FALSE
)

# =========================================================
# HEATMAP INTEGRATIVO
# =========================================================

heatmap_genes <- c(
  "SCHLAP1",
  "LINC01475",
  genes_eixo
)

mat_heat <- expr[, heatmap_genes]

mat_heat <- scale(mat_heat)

annotation <- data.frame(
  SCHLAP1 = expr$SCHLAP1_group,
  LINC01475 = expr$LINC01475_group
)

rownames(annotation) <- rownames(expr)

png(
  file.path(
    outdir,
    "Heatmap_Eixo_Oncogenico.png"
  ),
  width = 2600,
  height = 2200,
  res = 320
)

pheatmap(
  t(mat_heat),
  annotation_col = annotation,
  show_colnames = FALSE,
  fontsize = 11,
  clustering_method = "ward.D2",
  border_color = NA,
  color = colorRampPalette(
    rev(brewer.pal(11, "RdBu"))
  )(100)
)

dev.off()

# =========================================================
# FINAL
# =========================================================

cat(
  "\nAnálises do eixo oncogênico concluídas.\n"
)