# --- FUNÇÃO DE NORMALIZAÇÃO EXPANDIDA ---
raman_norm <- function(dataset, method = "min-max", target_peak = 520, search_window = 20) {
  
  # Extrair a matriz de dados (amostras nas linhas, wavenumbers nas colunas)
  matriz_dados <- dataset$data
  wavenumber <- dataset$wavenumbers
  
  if (method == "min-max") {
    message("A aplicar Normalização Min-Max...")
    # (x - min) / (max - min) para cada linha
    norm_data <- t(apply(matriz_dados, 1, function(x) {
      (x - min(x)) / (max(x) - min(x))
    }))
    
  } else if (method == "L2") {
    message("A aplicar Normalização L2 (Vetorial)...")
    # x / sqrt(sum(x^2)) para cada linha
    norm_data <- t(apply(matriz_dados, 1, function(x) {
      norma_l2 <- sqrt(sum(x^2))
      if(norma_l2 == 0) return(x) # Evitar divisão por zero
      x / norma_l2
    }))
    
  } else if (method == "silicon") {
    message("A aplicar Normalização pelo Pico de Silício (Eixo Y)...")
    
    # 1. Define os limites da janela de busca ao redor do pico alvo
    lower_bound <- target_peak - search_window
    upper_bound <- target_peak + search_window
    
    # 2. Cria a máscara lógica usando os wavenumbers do dataset
    mask <- (wavenumber > lower_bound) & (wavenumber < upper_bound)
    
    if (!any(mask)) {
      stop("Erro: Nenhum dado encontrado na janela de busca especificada para o Silício.")
    }
    
    # 3. Processa cada amostra (linha) para normalizar pelo valor máximo do pico de Si
    norm_data <- t(apply(matriz_dados, 1, function(x) {
      # Isola as intensidades da janela daquela amostra específica
      intensity_si <- x[mask]
      # Encontra o valor máximo de intensidade no pico do Silício
      silica_peak_value <- max(intensity_si, na.rm = TRUE)
      
      if (silica_peak_value == 0) return(x) # Evitar divisão por zero
      
      # Divide o espectro inteiro pelo valor desse pico
      x / silica_peak_value
    }))
    
  } else {
    stop("Método inválido. Escolha entre 'min-max', 'L2' ou 'silicon'.")
  }
  
  # Atualizar o dataset com os novos dados e retornar
  dataset$data <- norm_data
  return(dataset)
}