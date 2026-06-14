# =========================================================
# 25_PXD047127_validacao_proteomica_visualization.R
# Proteomic Validation Cohort - PXD047127
# GO + Reactome Visualization
# CORRIGIDO PARA OBJETOS clusterProfiler
# =========================================================

# =========================================================
# PACOTES
# =========================================================
library(tidyverse)
library(ggplot2)
library(patchwork)

# =========================================================
# DIRETÓRIO PRINCIPAL
# =========================================================
base_dir <- "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer"

# =========================================================
# DIRETÓRIO PROTEÔMICA
# =========================================================
proteomics_dir <- file.path(
  base_dir,
  "Validação_Proteômica"
)

# =========================================================
# RESULTS E FIGURES
# =========================================================
results_dir <- file.path(
  proteomics_dir,
  "results"
)

figures_dir <- file.path(
  proteomics_dir,
  "figures"
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# =========================================================
# ARQUIVOS DE ENRICHMENT
# =========================================================
go_file <- file.path(
  results_dir,
  "GO_enrichment_results.csv"
)

reactome_file <- file.path(
  results_dir,
  "Reactome_enrichment_results.csv"
)

# =========================================================
# VERIFICAR EXISTÊNCIA
# =========================================================
if (!file.exists(go_file)) {
  
  stop(
    paste(
      "\nArquivo GO não encontrado:\n",
      go_file
    )
  )
}

if (!file.exists(reactome_file)) {
  
  stop(
    paste(
      "\nArquivo Reactome não encontrado:\n",
      reactome_file
    )
  )
}

# =========================================================
# IMPORTAÇÃO
# =========================================================
cat("\n=====================================\n")
cat("IMPORTANDO ENRICHMENT RESULTS\n")
cat("=====================================\n")

go_df <- read.csv(
  go_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

reactome_df <- read.csv(
  reactome_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# =========================================================
# CONVERTER PARA DATA.FRAME PADRÃO
# =========================================================
go_df <- as.data.frame(go_df)

reactome_df <- as.data.frame(reactome_df)

# =========================================================
# VERIFICAR COLUNAS
# =========================================================
cat("\n=====================================\n")
cat("COLUNAS GO\n")
cat("=====================================\n")

print(colnames(go_df))

cat("\n=====================================\n")
cat("COLUNAS REACTOME\n")
cat("=====================================\n")

print(colnames(reactome_df))

# =========================================================
# GARANTIR COLUNAS NUMÉRICAS
# =========================================================
go_df$p.adjust <- as.numeric(go_df$p.adjust)
go_df$Count <- as.numeric(go_df$Count)

reactome_df$p.adjust <- as.numeric(reactome_df$p.adjust)
reactome_df$Count <- as.numeric(reactome_df$Count)

# =========================================================
# REMOVER NAs
# =========================================================
go_df <- go_df %>%
  filter(!is.na(p.adjust))

reactome_df <- reactome_df %>%
  filter(!is.na(p.adjust))

# =========================================================
# CONVERTER GeneRatio
# =========================================================
convert_ratio <- function(x) {
  
  sapply(x, function(y) {
    
    vals <- unlist(strsplit(y, "/"))
    
    as.numeric(vals[1]) /
      as.numeric(vals[2])
  })
}

go_df$GeneRatio_numeric <- convert_ratio(
  go_df$GeneRatio
)

reactome_df$GeneRatio_numeric <- convert_ratio(
  reactome_df$GeneRatio
)

# =========================================================
# TOP TERMS
# =========================================================
go_top <- go_df %>%
  arrange(p.adjust) %>%
  head(10)

reactome_top <- reactome_df %>%
  arrange(p.adjust) %>%
  head(10)

# =========================================================
# GO DOTPLOT
# =========================================================
go_plot <- ggplot(
  go_top,
  aes(
    x = GeneRatio_numeric,
    y = reorder(
      Description,
      GeneRatio_numeric
    )
  )
) +
  
  geom_point(
    aes(
      size = Count,
      color = -log10(p.adjust)
    ),
    alpha = 0.9
  ) +
  
  scale_color_gradient(
    low = "#9ecae1",
    high = "#08519c"
  ) +
  
  labs(
    title = "GO Biological Process Enrichment",
    subtitle = "PXD047127 Proteomic Cohort",
    x = "Gene Ratio",
    y = "",
    color = expression(-log[10](adjusted~p)),
    size = "Count"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    
    axis.title = element_text(
      face = "bold"
    )
  )

# =========================================================
# REACTOME DOTPLOT
# =========================================================
reactome_plot <- ggplot(
  reactome_top,
  aes(
    x = GeneRatio_numeric,
    y = reorder(
      Description,
      GeneRatio_numeric
    )
  )
) +
  
  geom_point(
    aes(
      size = Count,
      color = -log10(p.adjust)
    ),
    alpha = 0.9
  ) +
  
  scale_color_gradient(
    low = "#bdd7e7",
    high = "#2171b5"
  ) +
  
  labs(
    title = "Reactome Pathway Enrichment",
    subtitle = "PXD047127 Proteomic Cohort",
    x = "Gene Ratio",
    y = "",
    color = expression(-log[10](adjusted~p)),
    size = "Count"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    
    axis.title = element_text(
      face = "bold"
    )
  )

# =========================================================
# FIGURA COMBINADA
# =========================================================
combined_plot <- go_plot / reactome_plot

print(combined_plot)

# =========================================================
# EXPORTAR FIGURA COMBINADA
# =========================================================
ggsave(
  filename = file.path(
    figures_dir,
    "PXD047127_GO_Reactome_Combined.png"
  ),
  plot = combined_plot,
  width = 10,
  height = 12,
  dpi = 600
)

# =========================================================
# GO BARPLOT
# =========================================================
go_barplot <- ggplot(
  go_top,
  aes(
    x = reorder(
      Description,
      Count
    ),
    y = Count
  )
) +
  
  geom_col(
    fill = "#4A86B8",
    alpha = 0.9
  ) +
  
  coord_flip() +
  
  labs(
    title = "Top GO Biological Processes",
    x = "",
    y = "Protein Count"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    
    axis.title = element_text(
      face = "bold"
    )
  )

# =========================================================
# EXPORTAR GO BARPLOT
# =========================================================
ggsave(
  filename = file.path(
    figures_dir,
    "PXD047127_GO_Barplot.png"
  ),
  plot = go_barplot,
  width = 8,
  height = 6,
  dpi = 600
)

# =========================================================
# REACTOME BARPLOT
# =========================================================
reactome_barplot <- ggplot(
  reactome_top,
  aes(
    x = reorder(
      Description,
      Count
    ),
    y = Count
  )
) +
  
  geom_col(
    fill = "#3B6EA8",
    alpha = 0.9
  ) +
  
  coord_flip() +
  
  labs(
    title = "Top Reactome Pathways",
    x = "",
    y = "Protein Count"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    
    axis.title = element_text(
      face = "bold"
    )
  )

# =========================================================
# EXPORTAR REACTOME BARPLOT
# =========================================================
ggsave(
  filename = file.path(
    figures_dir,
    "PXD047127_Reactome_Barplot.png"
  ),
  plot = reactome_barplot,
  width = 8,
  height = 6,
  dpi = 600
)

# =========================================================
# TABELA RESUMO
# =========================================================
summary_terms <- bind_rows(
  
  go_top %>%
    dplyr::select(
      Description,
      p.adjust
    ) %>%
    mutate(Database = "GO"),
  
  reactome_top %>%
    dplyr::select(
      Description,
      p.adjust
    ) %>%
    mutate(Database = "Reactome")
)

write.csv(
  summary_terms,
  file.path(
    results_dir,
    "PXD047127_top_enrichment_terms.csv"
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
    "sessionInfo_PXD047127_visualization.txt"
  )
)

# =========================================================
# FINALIZAÇÃO
# =========================================================
cat("\n=====================================\n")
cat("VISUALIZAÇÃO CONCLUÍDA\n")
cat("=====================================\n")

cat("\nFiguras exportadas para:\n")
cat(figures_dir, "\n")