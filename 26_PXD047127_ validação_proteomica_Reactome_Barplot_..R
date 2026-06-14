# =========================================================
# PXD047127_Reactome_Barplot_PUBLICATION.R
# Publication-ready Reactome Barplot
# =========================================================

library(tidyverse)
library(ggplot2)
library(stringr)

# =========================================================
# DIRETÓRIOS
# =========================================================
base_dir <- "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer"

results_dir <- file.path(
  base_dir,
  "Validação_Proteômica",
  "results"
)

figures_dir <- file.path(
  base_dir,
  "Validação_Proteômica",
  "figures"
)

dir.create(figures_dir,
           recursive = TRUE,
           showWarnings = FALSE)

# =========================================================
# IMPORTAR DADOS
# =========================================================
reactome_file <- file.path(
  results_dir,
  "Reactome_enrichment_results.csv"
)

reactome_df <- read.csv(
  reactome_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# =========================================================
# LIMPEZA
# =========================================================
reactome_df <- reactome_df %>%
  
  mutate(
    Count = as.numeric(Count),
    p.adjust = as.numeric(p.adjust)
  ) %>%
  
  filter(
    !is.na(Count),
    !is.na(p.adjust)
  ) %>%
  
  arrange(p.adjust)

# =========================================================
# TOP 10
# =========================================================
reactome_top <- head(
  reactome_df,
  10
)

# =========================================================
# QUEBRAR NOMES LONGOS
# =========================================================
reactome_top$Description <- str_wrap(
  reactome_top$Description,
  width = 45
)

# =========================================================
# REORDENAR
# =========================================================
reactome_top$Description <- factor(
  reactome_top$Description,
  levels = rev(reactome_top$Description)
)

# =========================================================
# BARPLOT
# =========================================================
p <- ggplot(
  reactome_top,
  aes(
    x = Description,
    y = Count
  )
) +
  
  geom_col(
    fill = "#4A86B8",
    width = 0.75,
    alpha = 0.9
  ) +
  
  coord_flip() +
  
  labs(
    title = "Top Reactome Pathways",
    subtitle = "PXD047127 Proteomic Cohort",
    x = NULL,
    y = "Protein Count"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    
    plot.subtitle = element_text(
      size = 11
    ),
    
    axis.title.x = element_text(
      face = "bold",
      size = 13
    ),
    
    axis.text.y = element_text(
      size = 10
    ),
    
    axis.text.x = element_text(
      size = 10
    ),
    
    plot.margin = margin(
      t = 10,
      r = 20,
      b = 10,
      l = 20
    )
  )

# =========================================================
# MOSTRAR
# =========================================================
print(p)

# =========================================================
# EXPORTAR
# =========================================================
ggsave(
  filename = file.path(
    figures_dir,
    "PXD047127_Reactome_Barplot_Publication.png"
  ),
  plot = p,
  width = 11,
  height = 7,
  dpi = 600,
  bg = "white"
)

cat("\nFigura publication-ready salva com sucesso.\n")