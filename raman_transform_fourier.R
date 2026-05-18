raman_transform_fourier <- function(dataset) {
  message("A aplicar Fourier Transform (Power Spectrum)...")
  matriz_dados <- dataset$data
  
  power_spectrum <- t(apply(matriz_dados, 1, function(x) {
    # Aplica a Transformada Rápida de Fourier nativa do R
    fft_res <- fft(x)
    # Calcula a magnitude ao quadrado de cada componente de frequência (Power)
    # Retorna apenas a primeira metade devido à simetria do sinal real
    N <- length(x)
    mag <- (Mod(fft_res)^2) / N
    return(mag[1:floor(N/2)])
  }))
  
  dataset$data <- power_spectrum
  # Ajusta o vetor de frequências correspondente
  dataset$wavenumbers <- seq(0, 0.5, length.out = ncol(power_spectrum))
  
  return(dataset)
}