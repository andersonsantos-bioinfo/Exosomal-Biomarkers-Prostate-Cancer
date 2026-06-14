library(ggplot2)

# Garantir ordem dos grupos
risk_df$group <- factor(risk_df$group, levels = c("Normal", "Tumor"))

# Gráfico refinado
p <- ggplot(risk_df, aes(x = group, y = risk_score, fill = group)) +
  
  # Boxplot (estrutura principal)
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,   # evita duplicação visual dos outliers
    color = "black"
  ) +
  
  # Pontos individuais (mais elegantes)
  geom_jitter(
    width = 0.10,
    size = 1,
    alpha = 0.45,
    color = "black"
  ) +
  
  # Cores suaves e profissionais
  scale_fill_manual(values = c(
    "Normal" = "#4DBBD5",
    "Tumor"  = "#E64B35"
  )) +
  
  # Labels (sem título)
  labs(
    x = NULL,
    y = "Risk Score"
  ) +
  
  # Tema publication-ready
  theme_classic(base_size = 10, base_family = "Arial") +
  
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.title.y = element_text(size = 10),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

# Exibir
print(p)