
raman_radar_class_averages <- function(dataset) {
  message("A calcular médias por classe e a gerar projeções polares (Radar Plot)...")
  
  matriz_dados <- dataset$data
  wavenumbers <- dataset$wavenumbers
  classes <- as.factor(dataset$y_train)
  unique_classes <- levels(classes)
  N_wavenumbers <- length(wavenumbers)
  
  # Ângulos circulares para o mapeamento completo do espetro (0 a 2*pi)
  theta <- seq(0, 2 * pi, length.out = N_wavenumbers)
  
  # Configurar o layout gráfico para mostrar todas as classes lado a lado numa linha
  # Exemplo: se tens 5 classes, divide a janela em 1 linha e 5 colunas
  par(mfrow = c(1, length(unique_classes)))
  
  # Cores vibrantes para os contornos de cada classe
  cores_linhas <- c("#e63946", "#457b9d", "#1d3557", "#2a9d8f", "#f4a261")
  
  for (i in 1:length(unique_classes)) {
    classe_atual <- unique_classes[i]
    
    # Isolar as amostras pertencentes a esta classe e calcular a média por coluna
    indices_classe <- which(classes == classe_atual)
    espetro_medio <- colMeans(matriz_dados[indices_classe, , drop = FALSE])
    
    # Normalizar o raio entre 0 e 1 baseado nos limites globais do espetro
    r_norm <- (espetro_medio - min(matriz_dados)) / (max(matriz_dados) - min(matriz_dados))
    
    x_polar <- r_norm * cos(theta)
    y_polar <- r_norm * sin(theta)
    
    # Desenhar o gráfico polar quadrado
    plot(NULL, xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3), asp = 1,
         xlab = "", ylab = "", axes = FALSE,
         main = paste("Classe:", classe_atual))
    
    # Grelha interna de anéis circulares
    for(raio in c(0.2, 0.4, 0.6, 0.8, 1.0)) {
      pts_circulo <- seq(0, 2*pi, length.out = 100)
      lines(raio * cos(pts_circulo), raio * sin(pts_circulo), col = "gray90", lty = 1)
    }
    lines(1.05 * cos(pts_circulo), 1.05 * sin(pts_circulo), col = "black", lwd = 1)
    
    # Polígono preenchido com fundo escuro e glow transparente
    polygon(x_polar, y_polar, col = "#0d1b2a", border = NA)
    polygon(x_polar, y_polar, col = rgb(0.9, 0.1, 0.4, 0.2), border = NA)
    
    # Linha de contorno com a cor específica daquela classe
    lines(x_polar, y_polar, col = cores_linhas[i], lwd = 2)
    
    # Centro
    points(0, 0, col = "#00b4d8", pch = 16, cex = 1)
  }
  
  # Resetar o layout de ecrã do RStudio para o padrão 1x1
  par(mfrow = c(1, 1))
}