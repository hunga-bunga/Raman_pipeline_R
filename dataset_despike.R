#Função Despiking

dataset_despike <- function(dataset, ma = 10, threshold = 7) {
  
  # 1. Definição da lógica interna
  clean_single_spectrum <- function(y) {
    n <- length(y)
    med_y <- median(y)
    mad_y <- median(abs(y - med_y))
    
    # Proteção contra divisão por zero
    if (mad_y == 0) return(y)
    
    # Cálculo do Modified Z-Score
    y_z <- 0.6745 * (y - med_y) / mad_y
    spikes <- abs(y_z) > threshold
    spike_idx <- which(spikes)
    
    y_out <- y
    
    # Substituição dos spikes pela média da vizinhança
    for (i in spike_idx) {
      w <- seq(max(1, i - ma), min(n, i + ma))
      we <- w[!spikes[w]] # Usar apenas pontos que não são spikes para a média
      if (length(we) > 0) {
        y_out[i] <- mean(y[we])
      }
    }
    return(y_out)
  }
  
  message("A executar Despiking no dataset...")
  
  # 2. Aplicação à matriz inteira
  # t(apply(...)) garante que a estrutura [amostras x wavenumbers] se mantém
  dataset$data <- t(apply(dataset$data, 1, clean_single_spectrum))
  
  # Registo da operação para histórico
  cat("Concluído: ", sum(dataset$data != dataset$data), " potenciais spikes tratados.\n")
  
  return(dataset)
}