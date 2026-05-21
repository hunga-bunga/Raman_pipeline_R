raman_pca_plot <- function(dataset) {
  message("A calcular o PCA e a gerar o gráfico de clusters...")
  
  # 1. Executar o PCA na matriz de dados
  # scale. = TRUE é importante se os dados não estiverem perfeitamente normalizados
  pca_res <- prcomp(dataset$data, center = TRUE, scale. = FALSE)
  
  # 2. Calcular a variância explicada pelas duas primeiras componentes (PC1 e PC2)
  var_explicada <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
  lbl_pc1 <- paste0("PC1 (", round(var_explicada[1], 1), "%)")
  lbl_pc2 <- paste0("PC2 (", round(var_explicada[2], 1), "%)")
  
  # 3. Criar um data frame com as coordenadas dos pontos e as etiquetas originais
  df_pca <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    Classe = as.factor(dataset$y_train)
  )
  
  # 4. Desenhar o gráfico de dispersão (Scatter Plot) nativo do R
  # Definir cores bonitas para as tuas 5 classes
  cores_classes <- c("#e63946", "#457b9d", "#1d3557", "#2a9d8f", "#f4a261")
  
  plot(df_pca$PC1, df_pca$PC2, 
       col = cores_classes[df_pca$Classe], 
       pch = 19, cex = 1.2,
       xlab = lbl_pc1, ylab = lbl_pc2,
       main = "Análise de Componentes Principais (PCA) - Clusters Raman")
  
  # Adicionar uma grelha de fundo
  grid(col = "gray90", lty = "dotted")
  
  # Adicionar a legenda com os nomes reais das tuas classes
  legend("topright", legend = levels(df_pca$Classe), 
         col = cores_classes, pch = 19, bty = "n", cex = 0.9)
  
  return(pca_res)
}