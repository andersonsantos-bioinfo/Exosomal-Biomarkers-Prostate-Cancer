# ================================================
# Distribuição dos AUCs - 5-fold CV repetida 1000x
# Versão limpa - Arial 10, sem anotação de média
# ================================================

library(ggplot2)
library(tidyverse)

setwd("C:/Users/Administrator/Exossomos_prostata/Fase-Final/Validação_externa/GSE46602")

# Carregar o objeto model_cv (se não estiver carregado)
# load("model_cv_GSE46602.RData")   # descomente se necessário

# Extrair os valores de AUC de todas as repetições
auc_values <- model_cv$results$ROC

cat("Número de AUCs:", length(auc_values), "\n")
cat("AUC médio:", round(mean(auc_values), 4), "\n")
cat("Desvio padrão:", round(sd(auc_values), 4), "\n\n")

# Gráfico da distribuição
ggplot(data.frame(AUC = auc_values), aes(x = AUC)) +
  geom_histogram(binwidth = 0.008, 
                 fill = "#1F77B4", 
                 color = "white", 
                 alpha = 0.85) +
  labs(x = "AUC-ROC por repetição",
       y = "Frequência") +
  theme_minimal(base_size = 10) +
  theme(
    text = element_text(family = "Arial"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10)
  )

# Salvar
ggsave("C:/Users/Administrator/Exossomos_prostata/Fase-Final/Validação_externa/GSE46602/Distribuicao_AUCs_1000repeticoes.tiff",
       width = 9, 
       height = 6.5, 
       dpi = 300, 
       compression = "lzw")

cat("Gráfico da distribuição dos AUCs salvo com sucesso (Arial 10, sem legenda).\n")