#Alinhamento

align_silicon_peak <- function(wavenumber, intensity, target_peak = 520, search_window = 20) {
  
  # 1. Define os limites da janela de busca ao redor do pico alvo
  lower_bound <- target_peak - search_window
  upper_bound <- target_peak + search_window
  
  # 2. Cria a máscara lógica para filtrar a região do Silício
  mask <- (wavenumber > lower_bound) & (wavenumber < upper_bound)
  
  # Se nenhum ponto for encontrado na janela, interrompe com um erro
  if (!any(mask)) {
    stop("Nenhum dado encontrado na janela de busca especificada.")
  }
  
  wavenumber_si <- wavenumber[mask]
  intensity_si <- intensity[mask]
  
  # 3. Encontra a posição (wavenumber) onde ocorre a maior intensidade
  # which.max() equivale ao np.argmax() do Python
  peak_x <- wavenumber_si[which.max(intensity_si)]
  
  # 4. Computa o desvio (shift)
  shift <- target_peak - peak_x
  
  # 5. Corrige o eixo X inteiro somando o desvio
  aligned_wavenumber <- wavenumber + shift
  
  # Retorna uma lista contendo o novo eixo X corrigido e o valor do shift calculado
  return(list(
    aligned_wavenumber = aligned_wavenumber,
    shift = shift
  ))
}
