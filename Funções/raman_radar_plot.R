#Radar Gráfico/Visualização

raman_radar_plot <- function(dataset, sample_idx = 1) {
  message("A aplicar Transformação Radar e a gerar gráfico (Coordenadas Polares)...")
  
  matriz_dados <- dataset$data
  wavenumbers <- dataset$wavenumbers
  N_wavenumbers <- length(wavenumbers)
  
  # 1. Mapeia os wavenumbers uniformemente em ângulos de 0 a 2*pi
  theta <- seq(0, 2 * pi, length.out = N_wavenumbers)
  
  # 2. Converte as intensidades de todas as amostras para X e Y polares
  radar_xy_list <- list()
  for(idx in 1:nrow(matriz_dados)) {
    r <- matriz_dados[idx, ]
    
    # Normaliza o raio entre 0 e 1 para manter a proporção do círculo
    r_norm <- (r - min(r)) / (max(r) - min(r))
    if(any(is.na(r_norm))) r_norm <- rep(0, length(r))
    
    x_polar <- r_norm * cos(theta)
    y_polar <- r_norm * sin(theta)
    
    radar_xy_list[[idx]] <- data.frame(X_polar = x_polar, Y_polar = y_polar)
  }
  
  # 3. CONSTRUÇÃO AUTOMÁTICA DO GRÁFICO (Estilo da imagem)
  # Extrai os dados da amostra escolhida para o plot
  dados_plot <- radar_xy_list[[sample_idx]]
  
  # Configura a janela gráfica com proporção perfeitamente quadrada (asp = 1)
  plot(NULL, xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3), asp = 1,
       xlab = "", ylab = "", axes = FALSE,
       main = paste("Class Control - Average 3D Raman Spectra Projection\n(Amostra", sample_idx, ")"))
  
  # Desenha as linhas de grelha circulares (anéis concêntricos)
  valores_raio <- c(0.2, 0.4, 0.6, 0.8, 1.0)
  for(raio in valores_raio) {
    # Desenha o círculo gerando pontos de 0 a 2*pi
    pts_circulo <- seq(0, 2*pi, length.out = 200)
    lines(raio * cos(pts_circulo), raio * sin(pts_circulo), col = "gray85", lty = 1)
  }
  
  # Desenha o círculo exterior principal (a borda preta fina)
  lines(1.05 * cos(pts_circulo), 1.05 * sin(pts_circulo), col = "black", lwd = 1.2)
  
  # Coloca os marcadores de Wavenumbers ao redor do círculo (como na imagem)
  # Seleciona 8 pontos cardeais/colaterais para legendar
  idx_labels <- round(seq(1, N_wavenumbers, length.out = 9))[1:8]
  for(i in idx_labels) {
    # Calcula a posição do texto ligeiramente fora do círculo (raio = 1.15)
    ang <- theta[i]
    txt_x <- 1.18 * cos(ang)
    txt_y <- 1.18 * sin(ang)
    
    # Arredonda o valor do wavenumber para mostrar na borda
    text(txt_x, txt_y, labels = round(wavenumbers[i]), cex = 0.8, col = "black")
  }
  
  # Preenche a área do espetro com a cor sólida de fundo (Azul Escuro da imagem)
  polygon(dados_plot$X_polar, dados_plot$Y_polar, col = "#0d1b2a", border = NA)
  
  # Desenha o contorno vibrante por cima (Efeito degradé/fogo da imagem)
  # Fazemos múltiplas camadas transparentes para dar o efeito de brilho (glow)
  polygon(dados_plot$X_polar, dados_plot$Y_polar, col = rgb(0.9, 0.1, 0.4, 0.4), border = NA)
  polygon(dados_plot$X_polar * 0.95, dados_plot$Y_polar * 0.95, col = rgb(1, 0.6, 0, 0.3), border = NA)
  
  # Desenha a linha de contorno final externa
  lines(dados_plot$X_polar, dados_plot$Y_polar, col = "#e63946", lwd = 1.5)
  
  # Adiciona um ponto azul brilhante no centro exato (origem)
  points(0, 0, col = "#00b4d8", pch = 16, cex = 1.2)
  
  # Salva o resultado no formato do pipeline e retorna
  dataset$data_polar <- radar_xy_list
  return(dataset)
}