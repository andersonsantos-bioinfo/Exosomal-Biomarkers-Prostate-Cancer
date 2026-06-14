# =========================================================
# 27_GSE173094_suplementar_preprocessing.R
# Supplementary Validation Cohort
# Exosomal Biomarkers in Prostate Cancer
# =========================================================

# =========================================================
# PACOTES
# =========================================================
packages <- c(
  "tidyverse",
  "data.table",
  "limma",
  "caret",
  "ggplot2"
)

installed <- rownames(installed.packages())

for(pkg in packages){
  
  if(!pkg %in% installed){
    
    install.packages(pkg)
  }
  
  library(pkg,
          character.only = TRUE)
}

# =========================================================
# DIRETÓRIOS
# =========================================================
base_path <- "C:/Users/Administrator/Exossomos_prostata/Fase-Final/Validação_suplementar/GSE173094"

output_base <- "C:/Users/Administrator/Documents/Exosomal-Biomarkers-Prostate-Cancer/GSE173094"

processed_dir <- file.path(output_base, "processed")
figures_dir   <- file.path(output_base, "figures")
results_dir   <- file.path(output_base, "results")

dir.create(processed_dir,
           recursive = TRUE,
           showWarnings = FALSE)

dir.create(figures_dir,
           recursive = TRUE,
           showWarnings = FALSE)

dir.create(results_dir,
           recursive = TRUE,
           showWarnings = FALSE)

# =========================================================
# IMPORTAR MATRIZ LIMPA
# =========================================================
expr <- read.delim(
  file.path(
    base_path,
    "GSE173094_expr_matrix_clean.txt"
  ),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# =========================================================
# INSPEÇÃO
# =========================================================
cat("\n=====================================\n")
cat("DIMENSÃO ORIGINAL\n")
cat("=====================================\n")

print(dim(expr))

# =========================================================
# PRIMEIRA COLUNA = IDs
# =========================================================
rownames(expr) <- expr[,1]

expr <- expr[,-1]

# =========================================================
# CONVERTER NUMÉRICO
# =========================================================
expr <- data.frame(
  lapply(expr, function(x){
    
    suppressWarnings(
      as.numeric(as.character(x))
    )
    
  }),
  check.names = FALSE
)

# =========================================================
# MATRIZ
# =========================================================
expr_matrix <- as.matrix(expr)

mode(expr_matrix) <- "numeric"

# =========================================================
# REMOVER FEATURES TOTALMENTE NA
# =========================================================
expr_matrix <- expr_matrix[
  rowSums(is.na(expr_matrix)) < ncol(expr_matrix),
]

# =========================================================
# IMPUTAÇÃO
# =========================================================
expr_matrix[is.na(expr_matrix)] <- median(
  expr_matrix,
  na.rm = TRUE
)

# =========================================================
# VERIFICAR QUANTIS
# =========================================================
qx <- quantile(
  expr_matrix,
  probs = c(
    0,
    0.25,
    0.5,
    0.75,
    0.99,
    1
  ),
  na.rm = TRUE
)

print(qx)

# =========================================================
# LOG2
# =========================================================
LogC <- FALSE

if(
  !any(is.na(qx))
){
  
  LogC <- (
    qx[5] > 100 ||
      (qx[6] - qx[1] > 50)
  )
}

if(isTRUE(LogC)){
  
  cat("\nAplicando log2...\n")
  
  expr_matrix[expr_matrix <= 0] <- NA
  
  expr_matrix <- log2(expr_matrix)
  
  expr_matrix[is.na(expr_matrix)] <- median(
    expr_matrix,
    na.rm = TRUE
  )
}

# =========================================================
# NORMALIZAÇÃO QUANTILE
# =========================================================
expr_norm <- normalizeBetweenArrays(
  expr_matrix,
  method = "quantile"
)

# =========================================================
# REMOVER VARIÂNCIA ZERO
# =========================================================
nzv <- nearZeroVar(
  t(expr_norm)
)

if(length(nzv) > 0){
  
  expr_norm <- expr_norm[-nzv, ]
}

# =========================================================
# DIMENSÃO FINAL
# =========================================================
cat("\n=====================================\n")
cat("DIMENSÃO FINAL\n")
cat("=====================================\n")

print(dim(expr_norm))

# =========================================================
# EXPORTAR MATRIZ
# =========================================================
write.csv(
  expr_norm,
  file.path(
    processed_dir,
    "GSE173094_expression_processed.csv"
  )
)

# =========================================================
# PCA
# =========================================================
pca <- prcomp(
  t(expr_norm),
  scale. = TRUE
)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)

# =========================================================
# PCA PLOT
# =========================================================
p <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2
  )
) +
  
  geom_point(
    color = "#4A86B8",
    size = 3,
    alpha = 0.85
  ) +
  
  theme_classic(base_size = 13) +
  
  labs(
    title = "PCA - GSE173094",
    subtitle = "Supplementary Validation Cohort"
  )

print(p)

# =========================================================
# EXPORTAR FIGURA
# =========================================================
ggsave(
  file.path(
    figures_dir,
    "GSE173094_PCA_QC.png"
  ),
  p,
  width = 6,
  height = 5,
  dpi = 600
)

# =========================================================
# SALVAR WORKSPACE
# =========================================================
save(
  expr_norm,
  pca,
  pca_df,
  file = file.path(
    processed_dir,
    "GSE173094_preprocessed.RData"
  )
)

# =========================================================
# SESSION INFO
# =========================================================
writeLines(
  capture.output(sessionInfo()),
  file.path(
    results_dir,
    "sessionInfo_preprocessing.txt"
  )
)

# =========================================================
# FINALIZAÇÃO
# =========================================================
cat("\n=====================================\n")
cat("PREPROCESSAMENTO CONCLUÍDO\n")
cat("=====================================\n")