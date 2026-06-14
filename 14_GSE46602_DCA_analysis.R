library(tidyverse)

#-----------------------------------
# Preparação
#-----------------------------------
df <- risk_df %>%
  dplyr::select(group_bin, pred_final) %>%
  na.omit()

# thresholds adequados ao seu caso
thresholds <- seq(0.30, 0.95, by = 0.01)

# função Net Benefit
calc_nb <- function(thresh, y, p){
  pred_class <- ifelse(p >= thresh, 1, 0)
  
  TP <- sum(pred_class == 1 & y == 1)
  FP <- sum(pred_class == 1 & y == 0)
  N  <- length(y)
  
  (TP/N) - (FP/N) * (thresh / (1 - thresh))
}

#-----------------------------------
# Bootstrap
#-----------------------------------
set.seed(123)

B <- 500  # pode usar 1000 para artigo final

boot_results <- matrix(NA, nrow = B, ncol = length(thresholds))

for(i in 1:B){
  
  idx <- sample(1:nrow(df), replace = TRUE)
  df_boot <- df[idx, ]
  
  boot_results[i, ] <- sapply(thresholds, calc_nb,
                              y = df_boot$group_bin,
                              p = df_boot$pred_final)
}

#-----------------------------------
# Estatísticas
#-----------------------------------
nb_mean  <- apply(boot_results, 2, mean)
nb_lower <- apply(boot_results, 2, quantile, probs = 0.025)
nb_upper <- apply(boot_results, 2, quantile, probs = 0.975)

# treat-all
prev <- mean(df$group_bin)
treat_all <- prev - (1 - prev) * (thresholds / (1 - thresholds))

#-----------------------------------
# Dataframe final
#-----------------------------------
dca_ci <- data.frame(
  threshold = thresholds,
  mean = nb_mean,
  lower = nb_lower,
  upper = nb_upper,
  treat_all = treat_all,
  treat_none = 0
)

#-----------------------------------
# Plot PUBLICÁVEL
#-----------------------------------
library(ggplot2)

ggplot(dca_ci, aes(x = threshold)) +
  
  # IC sombreado
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "#1F77B4",
              alpha = 0.2) +
  
  # curva do modelo
  geom_line(aes(y = mean),
            color = "#1F77B4",
            linewidth = 1) +
  
  # treat-all
  geom_line(aes(y = treat_all),
            color = "#D62728",
            linewidth = 0.8) +
  
  # treat-none
  geom_line(aes(y = treat_none),
            color = "black",
            linewidth = 0.8) +
  
  labs(
    x = "Threshold probability",
    y = "Net benefit"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )