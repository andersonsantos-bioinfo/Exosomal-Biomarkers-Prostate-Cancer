# =========================================================
# 41.ASSOCIAÇÃO ENTRE
# SCHLAP1 / LINC01475 / EMAS
# E ISUP + T STAGE
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
library(scales)

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

dados$gleason_score <- as.numeric(
  dados$gleason_score
)

dados <- dados %>%
  
  filter(
    !is.na(gleason_score)
  )

# =========================================================
# CRIAR ISUP
# =========================================================

dados$ISUP <- case_when(
  
  dados$gleason_score <= 6 ~ "ISUP 1",
  
  dados$gleason_score == 7 ~ "ISUP 2-3",
  
  dados$gleason_score >= 8 ~ "ISUP 4-5",
  
  TRUE ~ NA_character_
)

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

dados$grupo_sch <- criar_grupo_tercil(
  dados$SCHLAP1
)

dados$grupo_linc <- criar_grupo_tercil(
  dados$LINC01475
)

dados$grupo_emas <- criar_grupo_tercil(
  dados$EMAS
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
# VERIFICAÇÃO
# =========================================================

table(dados$ISUP)

table(dados$grupo_sch)

table(dados$grupo_linc)

table(dados$grupo_emas)

# =========================================================
# TESTES QUI-QUADRADO - ISUP
# =========================================================

tab_isup_sch <- table(
  dados_sch$grupo_sch,
  dados_sch$ISUP
)

tab_isup_linc <- table(
  dados_linc$grupo_linc,
  dados_linc$ISUP
)

tab_isup_emas <- table(
  dados_emas$grupo_emas,
  dados_emas$ISUP
)

teste_isup_sch <- chisq.test(
  tab_isup_sch
)

teste_isup_linc <- chisq.test(
  tab_isup_linc
)

teste_isup_emas <- chisq.test(
  tab_isup_emas
)

# =========================================================
# TESTES QUI-QUADRADO - T STAGE
# =========================================================

dados_t_sch <- dados_sch %>%
  filter(!is.na(ajcc_pathologic_t))

dados_t_linc <- dados_linc %>%
  filter(!is.na(ajcc_pathologic_t))

dados_t_emas <- dados_emas %>%
  filter(!is.na(ajcc_pathologic_t))

tab_t_sch <- table(
  dados_t_sch$grupo_sch,
  dados_t_sch$ajcc_pathologic_t
)

tab_t_linc <- table(
  dados_t_linc$grupo_linc,
  dados_t_linc$ajcc_pathologic_t
)

tab_t_emas <- table(
  dados_t_emas$grupo_emas,
  dados_t_emas$ajcc_pathologic_t
)

teste_t_sch <- chisq.test(
  tab_t_sch
)

teste_t_linc <- chisq.test(
  tab_t_linc
)

teste_t_emas <- chisq.test(
  tab_t_emas
)

# =========================================================
# RESULTADOS
# =========================================================

teste_isup_sch
teste_isup_linc
teste_isup_emas

teste_t_sch
teste_t_linc
teste_t_emas

# =========================================================
# TEMA PUBLICAÇÃO
# =========================================================

tema_pub <- theme_classic(
  base_family = "Arial"
) +
  
  theme(
    
    text = element_text(size = 10),
    
    legend.position = "right",
    
    plot.title = element_blank(),
    
    axis.title = element_text(face = "plain"),
    
    axis.text = element_text(face = "plain")
  )

# =========================================================
# BARPLOT ISUP - SCHLAP1
# =========================================================

p1 <- ggplot(
  
  dados_sch,
  
  aes(
    x = grupo_sch,
    fill = ISUP
  )
  
) +
  
  geom_bar(
    position = "fill"
  ) +
  
  scale_y_continuous(
    labels = percent
  ) +
  
  xlab("SCHLAP1") +
  ylab("Proportion") +
  
  tema_pub

# =========================================================
# BARPLOT ISUP - LINC01475
# =========================================================

p2 <- ggplot(
  
  dados_linc,
  
  aes(
    x = grupo_linc,
    fill = ISUP
  )
  
) +
  
  geom_bar(
    position = "fill"
  ) +
  
  scale_y_continuous(
    labels = percent
  ) +
  
  xlab("LINC01475") +
  ylab("Proportion") +
  
  tema_pub

# =========================================================
# BARPLOT ISUP - EMAS
# =========================================================

p3 <- ggplot(
  
  dados_emas,
  
  aes(
    x = grupo_emas,
    fill = ISUP
  )
  
) +
  
  geom_bar(
    position = "fill"
  ) +
  
  scale_y_continuous(
    labels = percent
  ) +
  
  xlab("EMAS") +
  ylab("Proportion") +
  
  tema_pub

# =========================================================
# FIGURA FINAL ISUP
# =========================================================

fig_isup <- (
  p1 | p2 | p3
)

# visualizar
fig_isup

# =========================================================
# SALVAR FIGURA ISUP
# =========================================================

ggsave(
  
  filename = paste0(
    
    "C:/Users/Administrator/Exossomos_prostata/",
    "Fase-Final/Validação_suplementar/",
    "GSE173094/Figuras/",
    "ISUP_Association_TCGA.png"
  ),
  
  plot = fig_isup,
  
  width = 10,
  height = 4,
  dpi = 300
)

# =========================================================
# BARPLOTS T STAGE
# =========================================================

p4 <- ggplot(
  
  dados_t_sch,
  
  aes(
    x = grupo_sch,
    fill = ajcc_pathologic_t
  )
  
) +
  
  geom_bar(
    position = "fill"
  ) +
  
  scale_y_continuous(
    labels = percent
  ) +
  
  xlab("SCHLAP1") +
  ylab("Proportion") +
  
  labs(fill = "T stage") +
  
  tema_pub

# ---------------------------------------------------------

p5 <- ggplot(
  
  dados_t_linc,
  
  aes(
    x = grupo_linc,
    fill = ajcc_pathologic_t
  )
  
) +
  
  geom_bar(
    position = "fill"
  ) +
  
  scale_y_continuous(
    labels = percent
  ) +
  
  xlab("LINC01475") +
  ylab("Proportion") +
  
  labs(fill = "T stage") +
  
  tema_pub

# ---------------------------------------------------------

p6 <- ggplot(
  
  dados_t_emas,
  
  aes(
    x = grupo_emas,
    fill = ajcc_pathologic_t
  )
  
) +
  
  geom_bar(
    position = "fill"
  ) +
  
  scale_y_continuous(
    labels = percent
  ) +
  
  xlab("EMAS") +
  ylab("Proportion") +
  
  labs(fill = "T stage") +
  
  tema_pub

# =========================================================
# FIGURA FINAL T STAGE
# =========================================================

fig_tstage <- (
  p4 | p5 | p6
)

# visualizar
fig_tstage

# =========================================================
# SALVAR FIGURA T STAGE
# =========================================================

ggsave(
  
  filename = paste0(
    
    "C:/Users/Administrator/Exossomos_prostata/",
    "Fase-Final/Validação_suplementar/",
    "GSE173094/Figuras/",
    "TStage_Association_TCGA.png"
  ),
  
  plot = fig_tstage,
  
  width = 10,
  height = 4,
  dpi = 300
)

# =========================================================
# TABELA RESUMO - ISUP
# =========================================================

tabela_isup <- data.frame(
  
  Assinatura = c(
    "SCHLAP1",
    "LINC01475",
    "EMAS"
  ),
  
  Qui_quadrado = c(
    teste_isup_sch$statistic,
    teste_isup_linc$statistic,
    teste_isup_emas$statistic
  ),
  
  gl = c(
    teste_isup_sch$parameter,
    teste_isup_linc$parameter,
    teste_isup_emas$parameter
  ),
  
  p_valor = c(
    teste_isup_sch$p.value,
    teste_isup_linc$p.value,
    teste_isup_emas$p.value
  )
)

tabela_isup

# =========================================================
# TABELA RESUMO - T STAGE
# =========================================================

tabela_tstage <- data.frame(
  
  Assinatura = c(
    "SCHLAP1",
    "LINC01475",
    "EMAS"
  ),
  
  Qui_quadrado = c(
    teste_t_sch$statistic,
    teste_t_linc$statistic,
    teste_t_emas$statistic
  ),
  
  gl = c(
    teste_t_sch$parameter,
    teste_t_linc$parameter,
    teste_t_emas$parameter
  ),
  
  p_valor = c(
    teste_t_sch$p.value,
    teste_t_linc$p.value,
    teste_t_emas$p.value
  )
)

tabela_tstage