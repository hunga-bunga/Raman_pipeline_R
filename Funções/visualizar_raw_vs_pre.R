visualizar_raw_vs_pre <- function(df_raw, df_processado, metodo_norm = "L2") {
  # Configura a janela gráfica para 2 gráficos verticais
  par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
  
  # --- 1. GRÁFICO: DADOS ORIGINAIS (RAW) ---
  matriz_raw <- df_raw$data
  wavenumbers_raw <- df_raw$wavenumbers
  
  # Calcular métricas estatísticas dos dados brutos
  media_raw <- colMeans(matriz_raw, na.rm = TRUE)
  min_raw   <- apply(matriz_raw, 2, min, na.rm = TRUE)
  max_raw   <- apply(matriz_raw, 2, max, na.rm = TRUE)
  
  # Gráfico base vazio com os limites corretos do Raw
  plot(NULL, xlim = range(wavenumbers_raw), ylim = range(matriz_raw, na.rm = TRUE),
       xlab = "Wavenumber (cm-1)", ylab = "Intensidade (U.A.)",
       main = paste("Dataset ORIGINAL (Média Geral + Banda Total, N =", nrow(matriz_raw), ")"))
  grid(col = "gray85")
  
  # Desenhar a área de variabilidade (sombra cinzenta)
  polygon(c(wavenumbers_raw, rev(wavenumbers_raw)), c(max_raw, rev(min_raw)),
          col = rgb(0.8, 0.8, 0.8, alpha = 0.4), border = NA)
  
  # Desenhar a linha da média original a vermelho
  lines(wavenumbers_raw, media_raw, col = "red", lwd = 2.5)
  
  
  # --- 2. GRÁFICO: DADOS PRÉ-PROCESSADOS (DINÂMICO) ---
  matriz_proc <- df_processado$data
  wavenumbers_proc <- df_processado$wavenumbers
  
  # Calcular métricas estatísticas dos dados limpos
  media_proc <- colMeans(matriz_proc, na.rm = TRUE)
  min_proc   <- apply(matriz_proc, 2, min, na.rm = TRUE)
  max_proc   <- apply(matriz_proc, 2, max, na.rm = TRUE)
  
  # Definir cores e rótulos baseados no método escolhido
  if (toupper(metodo_norm) == "L2") {
    cor_linha  <- "darkblue"
    cor_sombra <- rgb(0.7, 0.8, 0.9, alpha = 0.5) # Sombra azulada para L2
    label_y    <- "Intensidade Normalizada L2"
    titulo_pos <- paste("Dataset PRÉ-PROCESSADO (Média Geral + Banda L2, N =", nrow(matriz_proc), ")")
  } else if (toupper(metodo_norm) == "MINMAX" || toupper(metodo_norm) == "MIN-MAX") {
    cor_linha  <- "darkgreen"
    cor_sombra <- rgb(0.7, 0.9, 0.7, alpha = 0.5) # Sombra esverdeada para MinMax
    label_y    <- "Intensidade Normalizada Min-Max [0, 1]"
    titulo_pos <- paste("Dataset PRÉ-PROCESSADO (Média Geral + Banda Min-Max, N =", nrow(matriz_proc), ")")
  } else {
    # Fallback genérico caso passes outro nome
    cor_linha  <- "black"
    cor_sombra <- rgb(0.6, 0.6, 0.6, alpha = 0.5)
    label_y    <- paste("Intensidade Normalizada", metodo_norm)
    titulo_pos <- paste("Dataset PRÉ-PROCESSADO (Média Geral + Banda", metodo_norm, ", N =", nrow(matriz_proc), ")")
  }
  
  # Desenhar o gráfico processado com as variáveis dinâmicas
  plot(NULL, xlim = range(wavenumbers_proc), ylim = range(matriz_proc, na.rm = TRUE),
       xlab = "Wavenumber (cm-1)", ylab = label_y, main = titulo_pos)
  grid(col = "gray85")
  
  # Desenhar a área de variabilidade pós-processamento
  polygon(c(wavenumbers_proc, rev(wavenumbers_proc)), c(max_proc, rev(min_proc)),
          col = cor_sombra, border = NA)
  
  # Desenhar a linha da média processada
  lines(wavenumbers_proc, media_proc, col = cor_linha, lwd = 2.5)
  
  # Restaurar definições de ecrã padrão do R
  par(mfrow = c(1, 1))
}