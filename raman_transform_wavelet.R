raman_transform_wavelet <- function(dataset, level = 1) {
  message("A aplicar Wavelet Transform (Haar DWT)...")
  matriz_dados <- dataset$data
  
  # Função interna para aplicar DWT Haar a um vetor de raiz
  dwt_haar_1d <- function(x) {
    n <- length(x)
    # Se o comprimento for ímpar, ajustamos o último elemento
    if (n %% 2 != 0) x <- c(x, x[n])
    
    # Filtro Haar de raiz
    approx <- (x[seq(1, n, by=2)] + x[seq(2, n, by=2)]) / sqrt(2)
    detail <- (x[seq(1, n, by=2)] - x[seq(2, n, by=2)]) / sqrt(2)
    
    return(c(approx, detail))
  }
  
  # Aplica o nível de decomposição desejado
  transform_matrix <- t(apply(matriz_dados, 1, function(row) {
    res <- row
    for(l in 1:level) {
      # Decompõe apenas a parte de aproximação a cada nível
      sub_n <- length(res) / (2^(l-1))
      res[1:sub_n] <- dwt_haar_1d(res[1:sub_n])
    }
    return(res)
  }))
  
  dataset$data <- transform_matrix
  return(dataset)
}