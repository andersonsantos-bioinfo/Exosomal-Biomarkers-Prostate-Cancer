# =========================================================
# 40.ASSOCIAÇÃO ENTRE
# SCHLAP1 / LINC01475 / EMAS
# E GLEASON + ISUP
# TCGA-PRAD
# =========================================================

# =========================================================
# PACOTES
# =========================================================

library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)

# =========================================================
# MATRIZ DE EXPRESSÃO
# =========================================================

matriz <- as.matrix(matriz_painel)

# garantir genes nas linhas
if(nrow(matriz) > ncol(matriz)){
  matriz <- t(matriz)
}

# =========================================================
# PADRONIZAR IDS TCGA
# =========================================================

samples_expr <- substr(
  colnames(matriz),
  1,
  12
)

# =========================================================
# DATAFRAME EXPRESSÃO
# =========================================================

df_expr <- data.frame(
  
  Sample = samples_expr,
  
  SCHLAP1 = as.numeric(
    matriz["SCHLAP1", ]
  ),
  
  LINC01475 = as.numeric(
    matriz["LINC01475", ]
  ),
  
  EMAS = as.numeric(EMAS)
)

# remover duplicados
df_expr <- df_expr %>%
  distinct(Sample, .keep_all = TRUE)

# =========================================================
# FUNÇÃO PARA TERCIS
# =========================================================

criar_grupo_tercil <- function(x){
  
  q1 <- quantile(
    x,
    0.33,
    na.rm = TRUE
  )
  
  q2 <- quantile(
    x,
    0.67,
    na.rm = TRUE
  )
  
  ifelse(
    x <= q1,
    "Low",
    
    ifelse(
      x >= q2,
      "High",
      "Mid"
    )
  )
}

# =========================================================
# CRIAR GRUPOS
# =========================================================

df_expr$grupo_sch <- criar_grupo_tercil(
  df_expr$SCHLAP1
)

df_expr$grupo_linc <- criar_grupo_tercil(
  df_expr$LINC01475
)

df_expr$grupo_emas <- criar_grupo_tercil(
  df_expr$EMAS
)

# =========================================================
# DADOS CLÍNICOS
# =========================================================

df_clin <- clinical %>%
  
  dplyr::select(
    
    submitter_id,
    
    gleason_score,
    tumor_grade,
    ajcc_pathologic_t
  ) %>%
  
  rename(
    Sample = submitter_id
  )

# =========================================================
# MERGE
# =========================================================

dados <- left_join(
  
  df_expr,
  df_clin,
  
  by = "Sample"
)

# =========================================================
# LIMPEZA
# =========================================================

dados <- dados %>%
  
  filter(
    !is.na(gleason_score)
  )

dados$gleason_score <- as.numeric(
  dados$gleason_score
)

# =========================================================
# REMOVER MID
# =========================================================

dados_sch <- dados %>%
  filter(grupo_sch != "Mid")

dados_linc <- dados %>%
  filter(grupo_linc != "Mid")

dados_emas <- dados %>%
  filter(grupo_emas != "Mid")

# =========================================================
# TESTES WILCOXON
# =========================================================

teste_sch <- wilcox.test(
  
  gleason_score ~ grupo_sch,
  
  data = dados_sch
)

teste_linc <- wilcox.test(
  
  gleason_score ~ grupo_linc,
  
  data = dados_linc
)

teste_emas <- wilcox.test(
  
  gleason_score ~ grupo_emas,
  
  data = dados_emas
)

# =========================================================
# RESULTADOS
# =========================================================

teste_sch
teste_linc
teste_emas

# =========================================================
# TEMA PUBLICAÇÃO
# =========================================================

tema_pub <- theme_classic(
  base_family = "Arial"
) +
  
  theme(
    
    text = element_text(size = 10),
    
    legend.position = "none",
    
    plot.title = element_blank(),
    
    axis.title = element_text(face = "plain"),
    
    axis.text = element_text(face = "plain")
  )

# =========================================================
# BOXPLOT SCHLAP1
# =========================================================

p1 <- ggboxplot(
  
  dados_sch,
  
  x = "grupo_sch",
  y = "gleason_score",
  
  fill = "grupo_sch",
  
  palette = c(
    "#4575B4",
    "#D73027"
  ),
  
  add = "jitter",
  size = 0.4
  
) +
  
  stat_compare_means(
    method = "wilcox.test"
  ) +
  
  xlab("SCHLAP1") +
  ylab("Gleason score") +
  
  tema_pub

# =========================================================
# BOXPLOT LINC01475
# =========================================================

p2 <- ggboxplot(
  
  dados_linc,
  
  x = "grupo_linc",
  y = "gleason_score",
  
  fill = "grupo_linc",
  
  palette = c(
    "#4575B4",
    "#D73027"
  ),
  
  add = "jitter",
  size = 0.4
  
) +
  
  stat_compare_means(
    method = "wilcox.test"
  ) +
  
  xlab("LINC01475") +
  ylab("Gleason score") +
  
  tema_pub

# =========================================================
# BOXPLOT EMAS
# =========================================================

p3 <- ggboxplot(
  
  dados_emas,
  
  x = "grupo_emas",
  y = "gleason_score",
  
  fill = "grupo_emas",
  
  palette = c(
    "#4575B4",
    "#D73027"
  ),
  
  add = "jitter",
  size = 0.4
  
) +
  
  stat_compare_means(
    method = "wilcox.test"
  ) +
  
  xlab("EMAS") +
  ylab("Gleason score") +
  
  tema_pub

# =========================================================
# FIGURA FINAL
# =========================================================

fig_final <- (
  p1 | p2 | p3
)

# visualizar
fig_final

# =========================================================
# SALVAR FIGURA
# =========================================================

ggsave(
  
  filename = paste0(
    
    "C:/Users/Administrator/Exossomos_prostata/",
    "Fase-Final/Validação_suplementar/",
    "GSE173094/Figuras/",
    "Gleason_Association_TCGA.png"
  ),
  
  plot = fig_final,
  
  width = 10,
  height = 4,
  dpi = 300
)

# =========================================================
# TABELA RESUMO
# =========================================================

tabela_resultados <- data.frame(
  
  Assinatura = c(
    "SCHLAP1",
    "LINC01475",
    "EMAS"
  ),
  
  W = c(
    teste_sch$statistic,
    teste_linc$statistic,
    teste_emas$statistic
  ),
  
  p_valor = c(
    teste_sch$p.value,
    teste_linc$p.value,
    teste_emas$p.value
  )
)

tabela_resultados